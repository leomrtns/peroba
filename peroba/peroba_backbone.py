import logging, ete3, argparse, treeswift
import matplotlib
import pkg_resources 
import numpy as np, pandas as pd

from utils import *
import common
import regression as ml 

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_backbone %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

current_working_dir = os.getcwd()

class PerobaBackbone:
    g_csv = None
    g_seq = None # dictionary
    l_csv = None
    l_seq = None # dictionary 
    trees = None ## these are treeswift trees
    # subset of columns from that may be useful (drop others)
    cols = ["sequence_name", "central_sample_id", "submission_org_code", "submission_org", "adm2",
            "collection_date", "country" , "cov_id", "sequencing_org", "sequencing_org_code", "sequencing_submission_date",
            "lineage", "lineage_support", "special_lineage","uk_lineage", "phylotype",
            "peroba_freq_acgt", "peroba_freq_n", "peroba_seq_uid", "source_age", "source_sex", ## until here  from global, below is local
            "adm2_private", "Repeat Sample ID", "icu_admission", "PCR Ct value", "No. Reads", "Mapped Reads", 
            "No. Bases (Mb)", "Coverage (X)", "Average read length", "Basic QC", "High Quality QC", "Missing bases (N)",
            "Consensus SNPs"] 

    def __init__ (self, peroba_db): # not split between global and local yet;
        self.g_csv = peroba_db[0]  # formatted metadata
        self.g_seq = {x.id:x for x in peroba_db[1]} # dictionary
        self.trees = [peroba_db[2]]  ## trees is a list; peroba_db[2] is treeswift already (assuming no dup names)
        cols = [x for x in self.cols if x in peroba_db[0].columns]
        self.g_csv = self.g_csv[cols] # remove other columns
        logger.info("Imported %s rows from database", str(self.g_csv.shape[0]))
    
    def add_local_data_and_sequences (self, csv=None, sequence=None, replace = False):
        if not sequence:
            logger.warning("Nothing to merge without sequences")
            self.split_data_sequences()
            return

        name_dict = {x:y for x,y in zip(self.g_csv["central_sample_id"], self.g_csv["sequence_name"])}
        seqs_in_global = [] # list of local unaligned sequences already in global

        logger.info("Updating sequence names if they are on global database")
        s_short = []; s_long = []; s_new = [] ## all will store original seq name
        for seq in sequence:
            if seq.id not in self.g_csv.index: # seqname is not a global (long) name
                if seq.id in name_dict.keys(): # sequence is present, but with short name
                    s_short.append (seq.id) # appends short name, while seqs_in_global[] has long name
                    seq.id = name_dict[seq.id] # receive coguk long name
                    seqs_in_global.append (seq.id) 
                else:
                    logger.warning("Sequence %s was not found in global database, will be added by hand", seq.id)
                    s_new.append (seq.id)
                    self.g_csv.loc[str(seq.id)] = pd.Series({
                        'sequence_name':seq.iq, 
                        'central_sample_id':seq.iq, 
                        'submission_org_code':"NORW",
                        'submission_org':"Norwich"})
            else: # sequence has long, official name
                if "NORW" in seq.id: ## we only consider replacing local seqs, otherwise database migth have newer  
                    s_long.append(seq.id)
                    seqs_in_global.append(seq.id)
                else:
                    seq.id = None
        logger.info("%s sequences found with long name and %s with short name on database. %s new sequences (not in database).",
                str(len(s_long)), str(len(s_short)), str(len(s_new)))

        sequence = {x.id:x for x in sequence if x.id is not None}
        if replace: # align all, which will replace existing ones in g_seq
            seq_list = [x for x in sequences.values()]
        else: # alignment will contain only those not found in global
            seq_list = [x for x in sequences.values() if x.id not in seqs_in_global]
        sqn = [x.id for x in seq_list]
        logger.info("Will update frequency info and align %s sequences", str(len(sqn)))
        self.g_csv.loc[ self.g_csv["sequence_name"].isin(sqn),"peroba_freq_acgt" ] = np.nan # recalculate
        self.g_csv.loc[ self.g_csv["sequence_name"].isin(sqn),"peroba_freq_n" ] = np.nan
        self.g_csv, sequence = add_sequence_counts_to_metadata (self.g_csv, sequence, from_scratch=False) # seqname must be in index

        # merge sequences (align local first, since global are already aligned)
        ref_seq = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas")  
        aln = minimap2_align_seqs(seq_list, reference_path=ref_seq)
        logger.info("Finished aligning %s local sequences", str(len(aln)))
        self.g_seq.update({x.id:x for x in aln})

        if csv:
            cols = [x for x in self.cols if x in csv.columns]
            csv = csv[cols] # remove other columns
            # global sequences with a match 
            csv = csv[ csv["central_sample_id"].isin(s_short + s_new) ] 
            logger.info("Local data contains %s rows matching sequences (excluding sequences matching long names in global database", 
                    str(csv.shape[0]))
            # merge metadata, to then spli into local and global datasets
            csv["peroba_seq_uid"] = csv["sequence_name"]
            csv.reset_index (drop=False, inplace=True) ## drop=False makes index (central_sample_id) become a column
            csv.set_index (self.g_csv.index.names, drop = True, inplace = True) # drop to avoid an extra 'peroba_seq_uid' column
            self.g_csv = common.df_merge_metadata_by_index (csv, self.g_csv) 
        self.split_data_sequences()
        return

    def split_data_sequences (self):
        self.l_csv = self.g_csv[  self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
        self.g_csv = self.g_csv[ ~self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
        self.l_seq = {x:y for x,y in self.g_seq.items() if x in self.l_csv["sequence_name"]}
        self.g_seq = {x:y for x,y in self.g_seq.items() if x in self.g_csv["sequence_name"]}
        logger.info ("split data into %s global and %s local sequences", str(self.g_csv.shape[0]), str(self.l_csv.shape[0]))

def read_peroba_database (f_prefix): 
    if f_prefix[-1] == ".": f_prefix = f_prefix[:-1] ## both `perobaDB.0621` and `perobaDB.0621.` are valid
    fname = f_prefix + common.suffix["metadata"]
    logger.info(f"Reading database metadata from \'{fname}\'")
    metadata = pd.read_csv (fname, compression="infer", index_col="peroba_seq_uid") 

    fname = f_prefix + common.suffix["tree"]
    logger.info(f"Reading database tree from \'{fname}\'")
    treestring = open(fname).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = treeswift.read_tree_newick (treestring) 

    fname = f_prefix + common.suffix["alignment"]
    logger.info(f"Reading database sequences from \'{fname}\'")
    sequences = common.read_fasta (fname, check_name = False)

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same \
            number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequences, tree]

def generate_backbone_dataset (database, csv, sequences, trees, replace, prefix):
    ## treeswift has node.set_label(label) to change name when we find NORW official name
    bb = PerobaBackbone (database)
    bb.add_local_data_and_sequences (csv, sequences, replace)
#    bb.add_trees (trees)

    
class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="""
    peroba_backbone is the script that generates a global backbone data set (COGUK+GISAID) given a local one (NORW).
    It depends on the prefix for a perobaDB set of files (from `peroba_database`), like "perobaDB.0519".
    It's recommended that you also local sequences, even without CSV metadata. You can furthermore add a newick file with extra 
    trees (the tree from previous run is a good choice).
    """, 
    usage='''peroba_backbone <perobaDB> [options]''')

    parser.add_argument('perobaDB')
    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
            help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-i', '--input', action="store", help="Directory where perobaDB files are. Default: working directory")
    parser.add_argument('-c', '--csv', metavar='csv', help="csv table with metadata from NORW")
    parser.add_argument('-s', '--sequences', metavar='fasta.bz2', help="extra sequences from NORW")
    parser.add_argument('-t', '--trees', metavar='', help="file with trees in newick format to help produce backbone")
    parser.add_argument('-o', '--output', action="store", help="Output database directory. Default: working directory")
    parser.add_argument('-r', '--replace', default=False, action='store_true', help="replace database sequence with local version")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: 
        output_d = os.path.join(current_working_dir, args.output)
        pathlib.Path(output_d).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    else: 
        output_d = current_working_dir
    prefix = os.path.join(output_d, "peroba_backbone." + datetime.datetime.now().strftime("%m%d_%H%M") + ".")

    if args.input: input_d = os.path.join(current_working_dir, args.input)
    else: input_d = current_working_dir

    logger.info("Reading metadata, sequences, and tree from peroba_database")
    database = read_peroba_database (os.path.join(input_d, args.perobaDB)) # something like "my_folder/perobaDB.0515"

    csv = None
    if (args.csv):
        fname = os.path.join(current_working_dir, args.csv)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.csv)
        if not os.path.exists(fname):
            logger.warning ("Could not find local CSV file {args.csv}; Will proceed without it")
        else:
            logger.info("Reading CSV file with metadata from NORW")
            csv = df_read_genome_metadata (fname, index_name = "central_sample_id")

    sequences = None
    if (args.sequences):
        fname = os.path.join(current_working_dir, args.sequences)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.sequences)
        if not os.path.exists(fname):
            logger.warning ("Could not find sequence file {args.sequences}; Will proceed without it")
        else:
            logger.info("Reading fasta file with sequences from NORW")
            sequences = common.read_fasta (fname, check_name = False)

    trees = None
    if (args.trees):
        fname = os.path.join(current_working_dir, args.trees)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.trees)
        if not os.path.exists(fname):
            logger.warning ("Could not find tree file {args.trees}; Will proceed without it")
        else:
            logger.info("Reading file with current trees and checking for duplicate names")
            treestring = [x.rstrip().replace("\'","").replace("\"","").replace("[&R]","") for x in open(fname)]
            trees = []
            for i,trs in enumerate (treestring): ## I use ete3 to remove duplicate leaves  
                tre = ete3.Tree(trs)
                tree_length = len([leaf.name for leaf in tre.iter_leaves()])
                tree_leaves = {str(leaf.name):leaf for leaf in tre.iter_leaves()} # dup leaves will simply overwrite node information
                if (tree_length > len(tree_leaves)):
                    tre.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
                    logger.warning(f"Found duplicated leaf names in input treefile {i}, will keep one at random")
                logger.info("%s leaves in treefile %s", len(tree_leaves), str(i))
                trees.append(tre)

    generate_backbone_dataset (database, csv, sequences, trees, args.replace, prefix)


if __name__ == '__main__':
    main()
