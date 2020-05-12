#!/usr/bin/env python

#import utils 
from utils import *

class DataSeqTree:   ## python3 doesn't use base class object anymore ;)
    """ Basic structure with all available data: sequence, tree, and epi-metadata, with file location information
    All the merge/append between data of the same type should be done elsewhere, since this class only drops data that is not
    at intersection. This class does check for duplicates, with arbitrary (aka not clever) removal. 
    ASAP I need to store the SNPalign 
    """
    db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
    timestamp = ""
    metadata = None
    tree     = None
    tree_leaves = None    # dictionary with {name:node}
    tree_needs_update = False # by default missing tree leaves do not remove samples (since we can reestimate tree) 
    sequences = None      # dictionary with sename:SeqRecord(); peroba.utils and other lowlevel use a lot of lists of SeqRecord()

    def __init__ (self, treefile=None, raw_metadata_filelist=None, sequence_filelist = None, input_path = None, peroba_db_path = None): 
        if peroba_db_path: self.db_path = peroba_db_path
        else: self.db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
        self.reset_timestamp()
        if (treefile): self.add_treefile (treefile, use_db_path=input_path)
        if (raw_metadata_filelist): self.add_metadata_from_raw_csv_files (raw_metadata_filelist, use_db_path=input_path)
        if (sequence_filelist): self.add_sequences (sequence_filelist, use_db_path=input_path)

    def reset_timestamp(self, timestamp = None):
        if (timestamp):  self.timestamp = timestamp
        else:            self.timestamp = datetime.datetime.now().strftime("prb%Y%m%d_%H%M%S") 

# TODO: decide if reestimate tree;  checkpoint/resume
    def add_treefile (self, tree=None, use_db_path=False): 
        if tree is None or self.tree is not None: ## LOG: tree already exists; refuse to replace it (I cant merge them)
            return
        if isinstance (tree, str):
            if use_db_path is True:
                treefile = f"{self.db_path}/{tree}"
            elif use_db_path:
                treefile = f"{use_db_path}/{tree}"
            tree = ete3.Tree(treefile)
        self.tree = tree
        # check duplicated names (ete3 accepts without checking, but complains later)
        tree_length = len([leaf.name for leaf in tre.iter_leaves()])
        self.tree_leaves = {str(leaf.name):leaf for leaf in tre.iter_leaves()} # dup leaves will simply overwrite node information
        if (tree_length > len(self.tree_leaves)):
            tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates

        if self.data is not None:      self.merge_data_tree ()
        if self.sequences is not None: self.merge_sequence_tree ()

    def add_metadata_from_clean_dataframe (self, metadata = None, use_db_path = False):
        """ Assume all information is already formatted (from previous run, or curated by hand)
        use_db_path can be a custom directory, or "True" if default DB is to be used
        """
        if metadata is None or self.metadata is not None:
            return
        if isinstance (metadata, str):
            if use_db_path is True:
                metadatafile = f"{self.db_path}/{metadata}"
            elif use_db_path:
                metadatafile = f"{use_db_path}/{metadata}"
            metadata = df_read_genome_metadata (metadatafile, index_name = "peroba_seq_uid")
            # I guess we can assume df.index.name is now "peroba_seq_uid", without duplicates
        self.metadata = metadata
        if self.tree is not None:      self.merge_data_tree ()
        if self.sequences is not None: self.merge_data_sequence ()

    def add_metadata_from_raw_csv_files (self, metadata_filelist = None, use_db_path = False):
        """ Given list of CSV/TSV files, merge them as possible. 
        If all files are in same directory structure, you can set `use_db_path`; o.w. please give relative paths for each file
        """
        if isinstance (metadata_filelist, str):
            metadata_filelist = [metadata_filelist] ## must be a list even if a single file 
        if use_db_path is True:
            dbpath = f"{self.db_path}/"
        elif use_db_path:
            dbpath = f"{use_db_path}/"
        else:
            dbpath = ""
        df1 = None
        for f in metadata_filelist:
            filepath = dbpath + f
            if "tsv" in f: sep = "\t"
            else:          sep = ","
            df2 = df_read_genome_metadata (filepath, sep = sep, index_name = "peroba_seq_uid")
            if df1 is None: df1 = df2
            else: df1 = df_merge_metadata_by_index (df1, df2) 
        df1 = df_finalise_metadata (df1)
        self.add_metadata_from_clean_dataframe (df1)

    def add_sequences (self, sequence_filelist = None, use_db_path = False):
        if isinstance (sequence_filelist, str):
            sequence_filelist = [sequence_filelist]
        if use_db_path is True:
            dbpath = f"{self.db_path}/"
        elif use_db_path:
            dbpath = f"{use_db_path}/"
        else:
            dbpath = ""
        self.sequences = dict(); 
        for f in sequence_filelist:
            filepath = dbpath + f
            if   "bz2" in f: zipmode = "bz2"
            elif "gz"  in f: zipmode = "gz"
            else:            zipmode = None
            seqs = read_fasta (filepath, zip = zipmode, check_name = True) # list
            self.sequences.update({x.id:x for x in seqs})  # dictionary of SeqRecord() (so duplicates are simply overwritten)
        if self.metadata is not None: self.merge_data_sequence ()
        if self.tree is not None:     self.merge_sequence_tree ()

    def merge_data_sequence (self): # keep only intersection
        if self.sequences is None or self.metadata is None:
            return
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_seqs = [x for x in self.sequences.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove sequences absent from the metadata
        for seqname in id_seqs:
            if seqname not in id_meta:
                del self.sequences[seqname] # delete this dict element
        # remove table rows without sequence information
        samples_found_in_sequences = self.metadata.index.isin([x for x in self.sequences.keys()])
        self.metadata = self.metadata([samples_found_in_tree])

    def merge_sequence_tree(self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.sequences is None:
            return
        id_seqs = [x for x in self.sequences.keys()]   # copy needed since we may delete elements 
        id_tree = [x for x in self.tree_leaves.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove leaves without a sequence
        purged_bool = False
        for leafname in id_tree:
            if leafname not in id_seqs:
                purged_bool = True
                del self.tree_leaves[leafname] # delete this dict element
        if purged_bool: # some leaves were removed
            tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates

        ## sequences not in the tree
        id_tree = [x for x in self.tree_leaves.keys()] 
        sequences_without_leaf = [x for x in self.sequences.keys() if x not in id_tree]
        if len(sequences_without_leaf): # or we remove sequences or we infer new tree (later)
            if use_tree_to_remove_samples: 
                for seq in sequences_without_leaf:
                    del self.sequences[seq]
            else:
                self.tree_needs_update = True 

    def merge_data_tree(self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.metadata is None:
            return
        ## remove leaves not in metadata
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_tree = [x for x in self.tree_leaves.keys()]   # copy needed since I'll delete elements (modifying dict)
        purged_bool = False
        for leafname in id_tree:
            if leafname not in id_meta:
                purged_bool = True
                del self.tree_leaves[leafname] # delete this dict element
        if purged_bool: # some leaves were removed
            tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates

        ## remove table rows absent from tree or just mark tree as incomplete
        samples_found_in_tree = self.metadata.index.isin([x for x in self.tree_leaves.keys()])
        if not all(samples_found_in_tree): # or we purge metadata rows or we infer new tree (later)
            if use_tree_to_remove_samples:  self.metadata = self.metadata([samples_found_in_tree])
            else:   self.tree_needs_update = True 
    
    def save (prefix=None, directory=None):
        if prefix is None: 
            self.reset_timestamp()
            prefix = str (self.timestamp)
        if directory is None:
            directory = f"{self.db_path}/"
        path = directory + prefix
        self.metadata.to_csv (f"{path}.metadata.csv.bz2")
        self.tree.write(format=1, outfile=f"{path}.tree.nhx")
        with bz2.open(f"{path}.unaligned.fasta.bz2","w") as fw: 
            for name, rec in self.sequences.items():
                fw.write(f">{name}\n{rec.seq}\n")

class SNPalign:
    """ This holds an aligment of SNP-only sites, to speed up calculations
        These can become characters for ancestral reconstruction.
    """
    reference = "China/Wuhan-Hu-1/2019"
    lineage = []
    snp = None

    def __init__(self, lineage = None, reference = None):
        self.lineage = [x for x in lineage] # sequence name, lineage, uk_lineage, etc

    def do_something (self, something=None):
        self.snp = None
