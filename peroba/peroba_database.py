#!/usr/bin/env python

## TODO: "UNKNOWN SOURCE" and "UNKNOWN" are the same adm2 (in cog)

import logging, ete3, argparse
#import numpy as np, pandas as pd, seaborn as sns
#from Bio import Seq, SeqIO, Align, AlignIO, Phylo, Alphabet, pairwise2
#from Bio.SeqRecord import SeqRecord
#import ete3 
#import datetime, time, codecs, sys, gzip, bz2, re, glob, pickle, collections, subprocess, os, errno, random, itertools, pathlib

from utils import *

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='perobaDB %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

cwd = os.getcwd()
#outdir = os.path.join(cwd, args.outdir)

class PerobaDatabase:   
    """ Basic structure with all available data: sequence, tree, and epi-metadata, with file location information
    All the merge/append between data of the same type should be done elsewhere: this class will drop data that does not
    have correspondence between data, seq, and tree.
    This class does check for duplicates, with arbitrary (aka not clever) removal.
    This is a "global" DataSequenceTree, which is not modified by peroba, and should be run whenever the upstream
    databases are updated (COGUK and GISAID). 
    """
    db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
    timestamp = ""
    metadata = None
    tree     = None
    tree_leaves = None    # dictionary with {name:node}
    sequences = None      # dictionary with sename:SeqRecord(); peroba.utils and other lowlevel use a lot of lists of SeqRecord()
    f_prefix = "perobaDB." # should never change, this file name is protected by the confraternity of bioinformaticians
    f_suffix = {
            "metadata": ".metadata.csv.gz",
            "tree": ".tree.nhx",
            "sequences": ".sequences.fasta.bz2"
            }

    def __init__ (self, from_db = None, timestamp = None, treefile=None, raw_metadata_filelist=None, 
            sequence_filelist = None, input_path = None, peroba_db_path = None): 
        if from_db is not None: 
            self.copy_from_PerobaDatabase (original = from_db)
            return
        if peroba_db_path: self.db_path = peroba_db_path
        else: self.db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
        self.reset_timestamp(timestamp) ## uses f_prefix
        logger.info(f"Starting database \'{self.timestamp}\'")
        if (treefile): 
            self.add_treefile (treefile, use_db_path=input_path)
        if (raw_metadata_filelist): 
            self.add_metadata_from_raw_csv_files (raw_metadata_filelist, use_db_path=input_path)
        if (sequence_filelist): 
            self.add_sequences (sequence_filelist, use_db_path=input_path)

    def reset_timestamp(self, timestamp = None):
        if (timestamp):  self.timestamp = timestamp
        else:            self.timestamp = datetime.datetime.now().strftime("%m%d") 

    def copy_from_PerobaDatabase (self, original = None, prefix = None): # untested
        if original is None: return
        if prefix is None and self.f_prefix is None: self.f_prefix = original.f_prefix 
        if prefix is not None: self.f_prefix = prefix
        self.db_path   = original.db_path
        self.timestamp = original.timestamp 
        logger.info(f"Copying database \'{self.timestamp}\'")
        self.metadata  = original.metadata.copy(deep=True)
        self.tree = ete3.Tree (original.tree)
        self.tree_leaves = copy.deepcopy (original.tree_leaves)
        self.sequences = copy.deepcopy (original.sequences) # not deep enough, we need to create new SeqRecord to be truly independent
 
    def read_from_db_files (self, timestamp=None, directory=None, subfolder=None): # subfolder is untested
        if directory is None: directory = self.db_path + "/"
        if subfolder is None: directory = directory
        else:                 directory = f"{directory}/{subfolder}/"
        if timestamp is None: # must find the most recent files; assumes default naming
            names = [x.replace(self.f_suffix["metadata"],"")[-4:] 
                    for x in glob.glob(f"{directory}/{self.f_prefix}*{self.f_suffix['metadata']}")] # careful with single and double quotes
            timestamp = str(sorted(names, key=int)[-1]) # highest timestamp => most recent?
            logger.debug("List of possible timestamps: %s", " ".join(names))

        self.timestamp = timestamp
        path = directory + self.f_prefix + timestamp # directory + file prefix
        logger.info(f"Reading database from \'{self.f_prefix}{timestamp}\' files")
        
        self.add_treefile (treefile=path + self.f_suffix["tree"])
        self.add_metadata_from_clean_dataframe (metadata = path + self.f_suffix["metadata"], use_db_path = False)
        self.add_sequences (sequence_filelist = path + self.f_suffix["sequences"], use_db_path = False)

    def save_to_db_files (self, timestamp=None, directory=None, subfolder=None):
        """ save your initial big table etc, since we may modify them
        """
        if timestamp is None: self.reset_timestamp()
        else:                 self.timestamp = str (timestamp)
        if directory is None: directory = self.db_path + "/"
        if subfolder is None: directory = directory
        else:                 directory = f"{directory}/{subfolder}/"

        prefix = self.f_prefix + self.timestamp
        pathlib.Path(directory).mkdir(parents=True, exist_ok=True) # python 3.5+
        path = directory + prefix # directory + file prefix("perobaDB") + timestamp
        logger.info(f"Saving database to files \'{prefix}\' on directory \'{directory}\'")
        self.metadata.to_csv (path + self.f_suffix["metadata"])
        self.tree.write(format=1, outfile=path + self.f_suffix["tree"])
        with bz2.open(path + self.f_suffix["sequences"],"wb") as fw: 
            for name, rec in self.sequences.items():
                if rec:  ## missing/query sequences
                    seq = str(rec.seq)
                    fw.write(str(f">{name}\n{seq}\n").encode())

    def add_treefile (self, treefile=None, use_db_path=False): 
        if (treefile is None or self.tree is not None): return
        if isinstance (treefile, str):
            if use_db_path is True: treefile = f"{self.db_path}/{treefile}"
            elif use_db_path:       treefile = f"{use_db_path}/{treefile}"
            else:                   treefile = treefile
        treestring = open(treefile).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
        self.tree = ete3.Tree(treestring)
        # check duplicated names (ete3 accepts without checking, but complains later)
        tree_length = len([leaf.name for leaf in self.tree.iter_leaves()])
        self.tree_leaves = {str(leaf.name):leaf for leaf in self.tree.iter_leaves()} # dup leaves will simply overwrite node information
        if (tree_length > len(self.tree_leaves)):
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.info("Read file %s with %s leaves", treefile, str(len(self.tree_leaves)))

        self.merge_data_tree ()
        self.merge_sequence_tree ()

    def add_metadata_from_clean_dataframe (self, metadata = None, use_db_path = False):
        """ Assume all information is already formatted (from previous run, or curated by hand)
        use_db_path can be a custom directory, or "True" if default DB is to be used
        """
        if (metadata is None or self.metadata is not None): return
        if isinstance (metadata, str):
            if use_db_path is True:
                metadatafile = f"{self.db_path}/{metadata}"
            elif use_db_path:
                metadatafile = f"{use_db_path}/{metadata}"
            else: 
                metadatafile = metadata
            metadata = df_read_genome_metadata (metadatafile, index_name = "peroba_seq_uid")
            # we can assume df.index.name is now "peroba_seq_uid", without duplicates
        self.metadata = metadata
        logger.info("Clean dataframe has shape %s (rows x cols)", str(self.metadata.shape))

        self.merge_data_tree ()
        self.merge_data_sequence ()

    def add_metadata_from_raw_csv_files (self, metadata_filelist = None, use_db_path = False):
        """ Given list of CSV/TSV files, merge them as possible. 
        If all files are in same directory structure, you can set `use_db_path`; o.w. please give relative paths for each file
        """
        if (self.metadata is not None): return
        if isinstance (metadata_filelist, str):
            metadata_filelist = [metadata_filelist] ## must be a list even if a single file 
        if use_db_path is True:
            dbpath = f"{self.db_path}/"
        elif use_db_path:
            dbpath = f"{use_db_path}/"
        else:
            dbpath = ""
        df1 = None
        logger.info(f"Will read metadata files from directory {dbpath} (if not set idividually, see below)")
        for f in metadata_filelist:
            filepath = dbpath + f
            if "tsv" in f: sep = "\t"
            else:          sep = ","
            logger.info(f"Reading metadata file {f}")
            df2 = df_read_genome_metadata (filepath, sep = sep, index_name = "peroba_seq_uid")
            if df1 is None: df1 = df2
            else: df1 = df_merge_metadata_by_index (df1, df2) 
        df1 = df_finalise_metadata (df1)
        logger.info("Finished reading input metadata files")
        self.add_metadata_from_clean_dataframe (df1)

    def add_sequences (self, sequence_filelist = None, use_db_path = False):
        if (self.sequences is not None): return
        if isinstance (sequence_filelist, str):
            sequence_filelist = [sequence_filelist]
        if use_db_path is True:
            dbpath = f"{self.db_path}/"
        elif use_db_path:
            dbpath = f"{use_db_path}/"
        else:
            dbpath = ""
        
        logger.info(f"Reading fasta sequence files from {dbpath} (if not set individually below)") 
        if self.sequences is None:  self.sequences = dict(); 
        for f in sequence_filelist:
            filepath = dbpath + f
            logger.info(f"Reading sequence file {f}")
            seqs = read_fasta (filepath, check_name = True) # list
            self.sequences.update({x.id:x for x in seqs})  # dictionary of SeqRecord() (so duplicates are simply overwritten)

        logger.info("Database now has %s valid sequences", str(len(self.sequences)))
        self.merge_data_sequence ()
        self.merge_sequence_tree ()

    def merge_DST (self, use_tree_to_remove_samples = False):
        logger.info("Full merge between sequences, metadata, and tree requested")
        self.merge_data_sequence ()
        self.merge_sequence_tree (use_tree_to_remove_samples = use_tree_to_remove_samples)
        self.merge_data_tree (use_tree_to_remove_samples = use_tree_to_remove_samples)

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
        logger.info("New dataframe size, after removing missing sequences: %s (r x c)", str(self.metadata.shape))
        # remove table rows without sequence information < REDUNDANT with df_match? >
        samples_found_in_sequences = self.metadata.index.isin([x for x in self.sequences.keys()])
        self.metadata = self.metadata[samples_found_in_sequences]
        logger.debug("Number of dataframe rows with sequence information = %s (should be same as above)", str(sum(samples_found_in_sequences)))

        #self.add_sequence_counts_to_metadata (from_scratch = False) ## add column if missing ## TODO: needed here or  downstream?

    def merge_sequence_tree (self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.sequences is None:
            return
        id_seqs = [x for x in self.sequences.keys()]   # copy needed since we may delete elements 
        id_tree = [x for x in self.tree_leaves.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove leaves without a sequence
        leaves_to_delete = set (id_tree) - set (id_seqs)
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        logger.info("Number of tree leaves to be deleted (no equivalent sequence): %s", str(len(leaves_to_delete)))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## sequences not in the tree
        if use_tree_to_remove_samples: 
            id_tree = [x for x in self.tree_leaves.keys()] 
            sequences_to_remove = set (id_seqs) - set (id_tree)
            if len(sequences_to_remove): # or we remove sequences or we infer new tree (later)
                logger.info("Number of sequences removed removed (not in the tree): %s", str(len(sequences_to_remove)))
                for seq in sequences_to_remove:
                    del self.sequences[seq]
        else:
            logger.info("Tree will not be used to remove sequences")

    def merge_data_tree (self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.metadata is None: return
        ## remove leaves not in metadata
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_tree = [x for x in self.tree_leaves.keys()]   # copy needed since I'll delete elements (modifying dict)
        leaves_to_delete = set(id_tree) - set(id_meta) # warning: loops can be quite slow
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        logger.info("Number of tree leaves to be deleted (no associated metadata): %s", str(len(leaves_to_delete)))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## remove table rows absent from tree or just mark tree as incomplete
        if use_tree_to_remove_samples:  
            samples_found_in_tree = self.metadata.index.isin([x for x in self.tree_leaves.keys()])
            if not all(samples_found_in_tree): # or we purge metadata rows or we infer new tree (later)
                logger.info("Metadata size after removing rows not found in the tree: %s", str(len(samples_found_in_tree)))
                self.metadata = self.metadata[samples_found_in_tree]
        else:   
            logger.info("Tree will not be used to remove rows from metadata")
    
    def add_sequence_counts_to_metadata (self, from_scratch = False):
        """ counts the proportions of indel/uncertain bases (i.e `N` and `-`) as well as number of ACGT.
        Notice that they do not sum up to one since we have partial uncertainty (`W`,`S`, etc.) that can be used
        TODO: This should not be done here, but in downstream analysis 
        """
        def calc_freq_N (index):
            if index not in self.sequences: return 1.
            if self.sequences[index] is None: return 1. ## missing sequences
            genome = str(self.sequences[index].seq); l = len(genome)
            if (l):
                number_Ns = sum([genome.upper().count(nuc) for nuc in ["N", "-"]])
                return number_Ns / l
            else: return 1.
        def calc_freq_ACGT (index):
            if index not in self.sequences: return 0.
            if self.sequences[index] is None: return 0. ## missing sequences
            genome = str(self.sequences[index].seq); l = len(genome)
            if (l):
                number_ACGTs = sum([genome.upper().count(nuc) for nuc in ["A", "C", "G", "T"]])
                return number_ACGTs / len(genome)
            else: return 0.

        if from_scratch or "peroba_freq_n" not in self.metadata.columns: # map() sends the index to lambda function
            self.metadata["peroba_freq_n"] = self.metadata.index.map(lambda x: calc_freq_N (x))  
        else: # only update null values 
            nilvalues = self.metadata[self.metadata["peroba_freq_n"].isnull()].index.map(lambda x: calc_freq_N (x))  
            if len(nilvalues) > 0:
                self.metadata[self.metadata["peroba_freq_n"].isnull()] = nilvalues 
        
        if from_scratch or "peroba_freq_acgt" not in self.metadata.columns: # map() sends the index to lambda function
            self.metadata["peroba_freq_acgt"] = self.metadata.index.map(lambda x: calc_freq_ACGT (x))  
        else: # only update null values 
            nilvalues = self.metadata[self.metadata["peroba_freq_acgt"].isnull()].index.map(lambda x: calc_freq_ACGT (x))  
            if len(nilvalues) > 0:
                self.metadata[self.metadata["peroba_freq_acgt"].isnull()] = nilvalues 
        
class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="peroba_database is the script that generates from scratch the integrated data from COGUK and GISAID",
    usage='''peroba_database [options] metadata (csv) files''')

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
            help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
#    parser.add_argument('-m', '--metadata', metavar='csv', type=argparse.FileType('r'), ## this opens the file !
    parser.add_argument('-m', '--metadata', metavar='csv', nargs='+', required=True, 
            help="csv/tsv files with metadata information (from COGUK or GISAID)")
    parser.add_argument('-t', '--tree', metavar='treefile', required=True,  help="single treefile in newick format")
    parser.add_argument('-s', '--fasta', metavar='fas', nargs='+', required=True, 
            help="fasta files with unaligned genomes from COGUK, GISAID")
    parser.add_argument('-o','--output', action="store", help="Output database directory. Default: working directory")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: output = os.path.join(cwd, args.output)
    else: output = cwd

    prb_db = PerobaDatabase (treefile=args.tree, raw_metadata_filelist=args.metadata, 
            sequence_filelist = args.fasta, input_path = None, peroba_db_path = output) 
    prb_db.save_to_db_files () 

if __name__ == '__main__':
    main()
