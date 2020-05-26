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

# TODO: we can create a "priority" column based on geography, to sort/preference

prefsort = [   # order of preference for samples (True => smaller values preferred)
    ['peroba_freq_acgt', False],
    [ 'peroba_freq_n', True],
    [ "PCR Ct value", True],
    [ "Coverage (X)", False],
    [ "submission_org", True],  # just to add preference for any value instead of nan (i.e. COGUK is preferred)
    [ 'lineage_support', False],
    ['collection_datetime', False] 
]

class PerobaBackbone:
    g_csv = None
    g_seq = None # dictionary
    l_csv = None
    l_seq = None # dictionary 
    trees = None ## these are treeswift trees
    # subset of columns from that may be useful (drop others)
    cols = ["sequence_name", "central_sample_id", "submission_org_code", "submission_org", "adm2", "collection_datetime", 
            "collection_date", "country" , "cov_id", "sequencing_org", "sequencing_org_code", "sequencing_submission_date",
            "lineage", "lineage_support", "special_lineage","uk_lineage", "phylotype",
            "peroba_freq_acgt", "peroba_freq_n", "peroba_seq_uid", "source_age", "source_sex", ## until here  from global, below is local
            "adm2_private", "Repeat Sample ID", "icu_admission", "PCR Ct value", "No. Reads", "Mapped Reads", 
            "No. Bases (Mb)", "Coverage (X)", "Average read length", "Basic QC", "High Quality QC", "Missing bases (N)",
            "Consensus SNPs"] 
    sort_cols = None 
    sort_ascend = None

    def __init__ (self, peroba_db): # not split between global and local yet;
        self.g_csv = peroba_db[0]  # formatted metadata
        self.g_seq = {x.id:x for x in peroba_db[1]} # dictionary
        self.trees = [peroba_db[2]]  ## trees is a list; peroba_db[2] is treeswift already (assuming no dup names)
        cols = [x for x in self.cols if x in peroba_db[0].columns]
        self.g_csv = self.g_csv[cols] # remove other columns
        
        pref = [x for x in prefsort if x[0] in self.g_csv.columns] # sort only existing columns
        self.sort_cols = [x[0] for x in pref] # follow this order (thus can't use dict...)
        self.sort_ascend = [x[1] for x in pref] # same order as sort_cols
        
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
                        'sequence_name':seq.id, 
                        'central_sample_id':seq.id, 
                        'submission_org_code':"NORW",
                        'submission_org':"Norwich"})
            else: # sequence has long, official name
                if "NORW" in seq.id: ## we only consider replacing local seqs, otherwise database migth have newer  
                    s_long.append(seq.id)
                    seqs_in_global.append(seq.id)
                else:
                    seq.id = None
        logger.info("%s long-named sequences and %s short-named sequences found on database. %s new sequences (not in database).",
                str(len(s_long)), str(len(s_short)), str(len(s_new)))
        logger.info("Now I'll check for obvious duplicate sequences, with identical names (only first will be used)")
        s_long = list(set(s_long)); s_short= list(set(s_short)); s_new = list(set(s_new))
        logger.info("%s long-named sequences and %s short-named sequences found on database. %s new sequences (not in database).",
                str(len(s_long)), str(len(s_short)), str(len(s_new)))

        # check for duplicates where we have both long and short names (short has precedence)
        logger.info("Checking for non-obvious duplicate sequences, where both short- and long-named versions appear")
        duplicate_long = []
        for seq in sequence:
            if seq.id in name_dict.keys() and name_dict[seq.id] in s_long:
                logger.warning("Duplicate: both %s and %s were found in local sequences (will keep first, with short name)",
                        seq.id, name_dict[seq.id]);
                duplicate_long.append(name_dict[seq.id])

        sequence = {x.id:x for x in sequence if x.id is not None and x.id not in duplicate_long}
        if replace: # align all, which will replace existing ones in g_seq
            seq_list = [x for x in sequence.values()]
        else: # alignment will contain only those not found in global
            seq_list = [x for x in sequence.values() if x.id not in seqs_in_global]
        sqn = [x.id for x in seq_list]
        logger.info("Will update frequency info and align %s sequences", str(len(sqn)))
        self.g_csv.loc[self.g_csv["sequence_name"].isin(sqn), "peroba_freq_acgt"] = np.nan # recalculate
        self.g_csv.loc[self.g_csv["sequence_name"].isin(sqn), "peroba_freq_n"] = np.nan
        logger.info("DEBUG")

        self.g_csv, sequence = common.add_sequence_counts_to_metadata (self.g_csv, sequence, from_scratch=False) # seqname must be in index

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

    def remove_seq_tree_based_on_metadata (self, seqnames = None):
        # only global sequences are removed 
        logger.info("removing global sequences after global metadata update")
        if seqnames is None: # use this variable carefully since it destroys the correspondence 
            seqnames = self.g_csv["sequence_name"].tolist()
        newseqs = {x:y for x,y in self.g_seq.items() if x in seqnames}
        # add all local seqnames to avoid being pruned 
        seqnames += self.l_csv["sequence_name"].tolist()
        newtrees = []
        for i,t in enumerate(self.trees):
            remove_leaf = []
            for l in t.traverse_leaves():
                lab = l.get_label()
                if lab not in seqnames:
                    remove_leaf.append(lab)
            if len(remove_leaf) > 0:
                logger.warning("In tree id %s, %s leaves were pruned (absent from global metadata)", str(i), str(len(remove_leaf)))
                newtrees.append(t.extract_tree_without(remove_leaf))
            else:
                newtrees.append(t)
        return  newseqs, newtrees

    def add_trees (self, trees = None):
        if trees is None or len(trees) < 1: return
        # leaf names will be "sequence_name", and we must translate them both on global and local
        name_d = {x:y for x,y in zip(self.g_csv["central_sample_id"], self.g_csv["sequence_name"])}
        name_d.update({x:y for x,y in zip(self.g_csv["sequence_name"], self.g_csv["sequence_name"])})
        name_d.update({x:y for x,y in zip(self.l_csv["central_sample_id"], self.l_csv["sequence_name"])})
        name_d.update({x:y for x,y in zip(self.l_csv["sequence_name"], self.l_csv["sequence_name"])})

        logger.info("Storing %s trees with leaves mapped to sequence names",str(len(trees)))
        for i,et in enumerate(trees):
            treestring = et.write() 
            remove_leaf = []
            t = treeswift.read_tree_newick (treestring)
            for l in t.traverse_leaves():
                lab = l.get_label()
                if lab in name_d.keys():
                    l.set_label(name_d[lab])
                else: # leaf is not in our set
                    remove_leaf.append(lab)
            if len(remove_leaf) > 0:
                logger.warning("In tree id %s the following leaves were not found and will be pruned:\n%s", str(i), "\n".join(remove_leaf))
                self.trees.append(t.extract_tree_without(remove_leaf))
            else:
                self.trees.append(t)

    def remove_duplicates (self, blocks = 4, leaf_size = 500, radius=0.00001):
        logger.info("Removing duplicates (identical sequences)")
        clusters = ml.list_duplicates (self.g_seq, blocks, leaf_size, radius)

        df = pd.DataFrame([[y,i] for i,x in enumerate(clusters) for y in x], columns=["sequence_name","peroba_tmp"])
        df = df.sort_values(by=["peroba_tmp"], ascending=True) # same sequence belongs to several clusters
        df = df.groupby("sequence_name").aggregate("first") # now each sequence belongs only to lowest cluster number
        logger.info("Number of unique sequences: %s", str(len(df["peroba_tmp"].unique())) )

        self.g_csv.reset_index (drop=False, inplace=True) ## merge will destroy index...
        df = self.g_csv.merge(df, on="sequence_name")
        df = df.sort_values(by=self.sort_cols, ascending=self.sort_ascend)
        df = df.groupby("peroba_tmp").aggregate("first") # only one sequence from each cluster, following self.order_col preference
        df.set_index ("peroba_seq_uid", drop = True, inplace = True)
        self.g_csv = df
        self.g_seqs, self.g_trees = self.remove_seq_tree_based_on_metadata()
    
    def find_neighbours (self, blocks = 2000, leaf_size = 500, dist_blocks = 1, nn = 20):
        logger.info(f"Finding neighbours to local sequences, using a distance of {dist_blocks} to {blocks} segments")
        neighbours = ml.list_r_neighbours (self.g_seq, self.l_seq, blocks, leaf_size, dist_blocks)
       
        blocks = 2 * blocks; leaf_size = leaf_size/2
        logger.info("Found %s neighbours; now will find their %s closest neighbours on %s segments", 
                len(neighbours), str(nn), str(blocks))
        aln_d = {x:self.g_seq[x] for x in neighbours}
        neighbours = ml.list_n_neighbours (self.g_seq, aln_d, blocks, leaf_size, nn)
        logger.info("Found %s neighbours", len(neighbours))
        return neighbours

def save_global_from_seqnames (bb, seqnames, prefix):
    seqs, trees = bb.remove_seq_tree_based_on_metadata(seqnames)
    save_sequences (seqs, prefix)
    fname = prefix + ".trees.nhx" 
    with open(fname,"w") as fw:
        for t in trees:
            fw.write(str(t))
    logger.info(f"Finished saving trees to file {fname}")

def save_sequences (seqs, prefix):
    fname = prefix + ".aln.xz" 
    logger.info(f"Saving sequences to file {fname}")
    mode = "wb"
    if   "bz2" in fname[-5:]: this_open = bz2.open
    elif "gz"  in fname[-5:]: this_open = gzip.open
    elif "xz"  in fname[-5:]: this_open = lzma.open
    else:  
        this_open = open
        mode = "w"
    with this_open(fname,mode) as fw: 
        for name, rec in seqs.items():
            if rec:  ## missing/query sequences
                seq = str(rec.seq)
                fw.write(str(f">{name}\n{seq}\n").encode())
                rec.id = name ## make sure alignment will have same names
    logger.info(f"Finished saving alignment")

def read_peroba_database (f_prefix): 
    if f_prefix[-1] == ".": f_prefix = f_prefix[:-1] ## both `perobaDB.0621` and `perobaDB.0621.` are valid
    fname = f_prefix + common.suffix["metadata"]
    logger.info(f"Reading database metadata from \'{fname}\'")
    metadata = pd.read_csv (fname, compression="infer", index_col="peroba_seq_uid") 
    metadata = common.df_finalise_metadata (metadata) 

    fname = f_prefix + common.suffix["tree"]
    logger.info(f"Reading database tree from \'{fname}\'")
    treestring = open(fname).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = treeswift.read_tree_newick (treestring) 

    fname = f_prefix + common.suffix["alignment"]
    logger.info(f"Reading database sequences from \'{fname}\'")
    sequences = common.read_fasta (fname, check_name = False)

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequences, tree]

def generate_backbone_dataset (database, csv, sequences, trees, replace, prefix):
    # initialisation
    bb = PerobaBackbone (database)
    bb.add_local_data_and_sequences (csv, sequences, replace)
    bb.add_trees (trees)
    bb.remove_duplicates()
    # more methods come here
    neighbours = bb.find_neighbours()
    save_global_from_seqnames (bb, neighbours, prefix + "coguk_NN")
    # finally, save all NORW sequences
    save_sequences (bb.l_seq, prefix + "norw")
    
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
            csv = common.df_read_genome_metadata (fname, index_name = "central_sample_id")
            csv = common.df_finalise_metadata (csv)

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
