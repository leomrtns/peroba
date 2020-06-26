#!/usr/bin/env python

import logging, ete3, argparse

from peroba.utils import *
from peroba import common 

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='perobaDB %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

current_working_dir = os.getcwd()

class PerobaDatabase:   
    """ Basic structure with all available data: sequence, tree, and epi-metadata, with file location information
    All the merge/append between data of the same type should be done elsewhere: this class will drop data that does not
    have correspondence between data, seq, and tree.
    This class does check for duplicates, with arbitrary (aka not clever) removal.
    This is a "global" DataSequenceTree, which is not modified by peroba, and should be run whenever the upstream
    databases are updated (COGUK and GISAID). 
    """
    metadata = None
    raw = None  ## also metadata, but with _all_ rows, even duplicates or missing sequences
    tree = None
    no_pruning = True
    tree_leaves = None    # dictionary with {name:node}
    sequences = None      # dictionary with sename:SeqRecord(); peroba.utils and other lowlevel use a lot of lists of SeqRecord()
    prefix = None

    def __init__ (self, tree, metadata, sequence, prefix, pruning = False): 
        self.prefix = prefix
        if pruning is True: self.no_pruning = False 
        else: self.no_pruning = True
        logger.info(f"Starting database \'{prefix}\'")
        if (len(tree)):     self.add_treefile (tree)
        if (len(metadata)): self.add_metadata (metadata)
        if (len(sequence)): self.add_sequences (sequence)

    def save_to_db_files (self, alignment):
        logger.info("Full merge between sequences, metadata, and tree before saving database")
        self.merge_data_sequence ()
        self.merge_sequence_tree ()
        self.merge_data_tree ()

        if self.metadata is not None: 
            fname = self.prefix + common.suffix["metadata"]
            logger.info(f"Saving metadata to file {fname}")
            self.metadata.to_csv (fname)
            desc = "perobaDB merged metadata information {}".format(datetime.datetime.now().strftime("%Y-%m-%d"))
            fname = self.prefix + ".html"
            logger.info(f"Generating an HTML report of metadata variables in {fname}")
            common.metadata_to_html (self.metadata, fname, desc) 

        else:
            logger.warning("Metadata missing: this won't be a valid database")

        if self.raw is not None: 
            fname = self.prefix + common.suffix["raw"]
            logger.info(f"Saving raw (metadata) table to file {fname}")
            self.raw.to_csv (fname)
        else:
            logger.warning("Metadata missing: this won't be a valid database")

        if self.tree is not None: 
            fname = self.prefix + common.suffix["tree"]
            logger.info(f"Saving tree to file {fname}")
            self.tree.write(format=1, outfile=fname)
        else:
            logger.warning("Tree missing: this won't be a valid database")

        if self.sequences is not None: 
            fname = self.prefix + common.suffix["sequences"]
            logger.info(f"Saving unaligned sequences to file {fname}")
            mode = "wb"
            if   "bz2" in fname[-5:]: this_open = bz2.open
            elif "gz"  in fname[-5:]: this_open = gzip.open
            elif "xz"  in fname[-5:]: this_open = lzma.open
            else:  
                this_open = open
                mode = "w"
            with this_open(fname, mode) as fw: 
                for name, rec in self.sequences.items():
                    if rec:  ## missing/query sequences
                        seq = str(rec.seq)
                        fw.write(str(f">{name}\n{seq}\n").encode())
                        rec.id = name ## make sure alignment will have same names

            aln = self.update_alignment (alignment, seqs_per_block = 1000)
            fname = self.prefix + common.suffix["alignment"]
            logger.info(f"Saving alignment to file {fname}")
            with this_open(fname,"wb") as fw: 
                for sequence in aln: # list of sequences (not dictionary as above) 
                    fw.write(str(f">{sequence.id}\n{sequence.seq}\n").encode())
            logger.info(f"Finished saving alignment")
        else:
            logger.warning("Sequences missing: this won't be a valid database")

    def add_treefile (self, tree): 
        self.tree = common.read_ete_treefile (tree[0], multi = False)
        self.tree_leaves = {str(leaf.name):leaf for leaf in self.tree.iter_leaves()} # dup leaves will simply overwrite node information

        self.merge_data_tree ()
        self.merge_sequence_tree ()

    def add_metadata (self, metadata):
        """ Given list of CSV/TSV files, merge them as possible. 
        """
        df1 = None
        for f in metadata:
            if "tsv" in f[-8:]: sep = "\t" # tsv or csv will be determined by suffix
            else:               sep = ","
            logger.info(f"Reading metadata file {f}")
            df2 = common.df_read_genome_metadata (f, sep = sep, index_name = "peroba_seq_uid")
            if df1 is None: df1 = df2
            else: df1 = common.df_merge_metadata_by_index (df1, df2) 
        df1 = common.df_finalise_metadata (df1)
        self.metadata = df1
        self.raw = df1.copy() ## not merged with sequence
        logger.info("Finished reading metadata files. Merged (raw) dataframe has shape %s (rows x cols)", str(self.metadata.shape))
        
        self.merge_data_tree ()
        self.merge_data_sequence ()

## TODO: allow for merging of NORW sequences maybe? or add merged file on report 

    def add_sequences (self, sequence): 
        logger.info(f"Reading fasta sequence files (if not set individually below)") 
        self.sequences = dict(); 
        for f in sequence:
            logger.info(f"Reading sequence file {f}")
            seqs = common.read_fasta (f, check_name = True) # list
            self.sequences.update({x.id:x for x in seqs})  # dictionary of SeqRecord() (so duplicates are simply overwritten)

        logger.info("Database now has %s valid sequences", str(len(self.sequences)))
        self.merge_data_sequence ()
        self.merge_sequence_tree ()

    def update_alignment (self, alignment, seqs_per_block):
        ref_seq = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas")  
        if alignment is None or len(alignment) < 1:
            logger.info(f"Aligning all sequences with mafft (no alignment file found)")
            aln = common.align_sequences_in_blocks ( [x for x in self.sequences.values()], reference_file = ref_seq, seqs_per_block=seqs_per_block)
            return aln

        aln = dict() # aln is dict but prealign is list
        for f in alignment:
            logger.info(f"Reading Alignment file {f}")
            seqs = common.read_fasta (f, check_name = True) # list
            aln.update({x.id:x for x in seqs})  #  duplicates are overwritten
        
        prealign = [x for x in aln.values() if x.id in self.sequences.keys()] 
        prealn_names = [x.id for x in prealign]
        remain = [x for x in self.sequences.values() if x.id not in prealn_names]
        logger.info(f"From %s sequences, %s were found in alignment (originally with %s sequences)", 
                len(self.sequences), len(prealign), len(aln))
        logger.info(f"Aligning remaining sequences with mafft")
        aln_list = common.align_sequences_in_blocks (remain, reference_file = ref_seq, seqs_per_block=seqs_per_block)
        return prealign + aln_list

    def merge_data_sequence (self): # keep only intersection
        if self.sequences is None or self.metadata is None: return
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_seqs = [x for x in self.sequences.keys()] # copy needed since I'll delete elements (modifying dict)

        ## remove sequences absent from the metadata
        idx_sequence_only = set (id_seqs) - set (id_meta) # names unique to sequence, not found in metadata
        idx_matched = set(id_seqs) & set(id_meta)
        df_matched = self.metadata.loc[idx_matched,] ## these rows have well-behaved sequences, with same name 

        # sequences found will be new SeqRecords, old ones can be removed
        for seqname in idx_sequence_only:
            del self.sequences[seqname] # delete this dict element
        self.metadata = df_matched
        logger.info("DataSeq:: New dataframe size, after removing missing sequences: %s (r x c)", str(self.metadata.shape))
        # remove table rows without sequence information < REDUNDANT with df_match? >
        samples_found_in_sequences = self.metadata.index.isin([x for x in self.sequences.keys()])
        self.metadata = self.metadata[samples_found_in_sequences]
        logger.debug("DataSeq:: Number of dataframe rows with sequence information = %s (should be same as above)", str(sum(samples_found_in_sequences)))
        self.metadata, self.sequences = common.add_sequence_counts_to_metadata (self.metadata, self.sequences, from_scratch = False) ## add column if missing 
        self.raw = self.raw.combine_first (self.metadata[["peroba_freq_n", "peroba_freq_acgt"]]) ## add columns by index

    def merge_sequence_tree (self): # keep only intersection
        if self.no_pruning or self.tree is None or self.sequences is None: return
        id_seqs = [x for x in self.sequences.keys()]   # copy needed since we may delete elements 
        id_tree = [x for x in self.tree_leaves.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove leaves without a sequence
        leaves_to_delete = set (id_tree) - set (id_seqs)
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        logger.info("SeqTree:: Number of tree leaves to be deleted (no equivalent sequence): %s", str(len(leaves_to_delete)))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## sequences not in the tree
        logger.debug("SeqTree:: Tree will not be used to remove sequences")

    def merge_data_tree (self): # keep only intersection
        if self.no_pruning or self.tree is None or self.metadata is None: return
        ## remove leaves not in metadata
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_tree = [x for x in self.tree_leaves.keys()]   # copy needed since I'll delete elements (modifying dict)
        leaves_to_delete = set(id_tree) - set(id_meta) # warning: loops can be quite slow
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        logger.info("DataTree:: Number of tree leaves to be deleted (no associated metadata): %s", str(len(leaves_to_delete)))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## remove table rows absent from tree or just mark tree as incomplete
        logger.debug("DataTree:: Tree will not be used to remove rows from metadata")

class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def fill_file_location (filelist, directory):
    if isinstance (filelist, str):
        filelist = [filelist]
    fullname = []
    for f in filelist: ## FIXME: priority is directory not working_dir
        fname = os.path.join(current_working_dir, f)
        if not os.path.exists(fname):
            fname = os.path.join(directory, f)
        if not os.path.exists(fname):
            logger.warning (f"Could not find file {f}; Will proceed anyway")
        else:
            logger.info(f"Found file {fname}")
            fullname.append(fname)
    return fullname 

def main():
    parser = ParserWithErrorHelp(
    description="peroba_database is the script that generates from scratch the integrated data from COGUK and GISAID",
    usage='''peroba_database [options] ''')

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
        help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
#    parser.add_argument('-m', '--metadata', metavar='csv', type=argparse.FileType('r'), ## this opens the file !
    parser.add_argument('-m', '--metadata', metavar='csv', nargs='+', required=True, 
        help="csv/tsv files with metadata information (from COGUK or GISAID)")
    parser.add_argument('-t', '--tree', metavar='treefile', required=True,  help="single treefile in newick format")
    parser.add_argument('-s', '--fasta', metavar='fas', nargs='+', required=True, 
        help="fasta files with unaligned genomes from COGUK, GISAID")
    parser.add_argument('-a', '--alignment', metavar='aln', nargs="+", help="aligned sequences from previous iteration (speeds up calculations)")
    parser.add_argument('-i', '--input', action="store", help="Directory where input files are. Default: working directory")
    parser.add_argument('-o', '--output', action="store", help="Output database directory. Default: working directory")
    parser.add_argument('-f', '--force-pruning', dest = "force", default=False, action='store_true', 
        help="Prune tree, removing leaves absent from metadata and sequences. Slow, few benefits?")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: 
        output_d = os.path.join(current_working_dir, args.output)
        pathlib.Path(output_d).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    else: 
        output_d = current_working_dir
    prefix = os.path.join(output_d, common.prefix["database"] + datetime.datetime.now().strftime("%m%d"))

    if args.input: input_d = os.path.join(current_working_dir, args.input)
    else: input_d = current_working_dir

    tree = fill_file_location (args.tree, input_d)
    meta = fill_file_location (args.metadata, input_d)
    fasta = fill_file_location (args.fasta, input_d)
    align = fill_file_location (args.alignment, input_d)

    prb_db = PerobaDatabase (tree = tree, metadata = meta, sequence = fasta, prefix = prefix, pruning = args.force) 
    prb_db.save_to_db_files (align) ## alignment added only here 

if __name__ == '__main__':
    main()
