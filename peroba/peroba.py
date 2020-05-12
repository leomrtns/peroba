#!/usr/bin/env python

import utils 

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
    tree_leaves = None # dictionary with {name:node}
    sequence = None

    def __init__ (self, db_path = None):
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        if db_path: self.db_path = db_path

    def reset_timestamp(self, timestamp = None):
        if (timestamp):  self.timestamp = timestamp
        else:            self.timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") 

    def add_treefile (self, tree=None, use_db_path=False): 
        if tree is None or self.tree: ## LOG: tree already exists; refuse to replace it (I cant merge them)
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

        if self.data:     self.merge_data_tree ()
        if self.sequence: self.merge_sequence_tree ()
        
    def merge_data_tree(self): # keep only intersection
        if not self.tree or not self.metadata:
            return
        ## remove leaves not in metadata
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_tree = self.tree_leaves.keys()   # copy needed since I'll delete elements (modifying dict)
        purged_bool = False
        for leafname in id_tree:
            if leafname not in id_meta:
                purged_bool = True
                del self.tree_leaves[leafname] # delete this dict element
        if purged_bool:
            tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates


    def add_metadata_from_clean_dataframe (self, metadata = None, use_db_path = False):
        """ Assume all information is already formatted (from previous run, or curated by hand)
        use_db_path can be a custom directory, or "True" if default DB is to be used
        """
        if metadata is None or self.metadata:
            return
        if isinstance (metadata, str):
            if use_db_path is True:
                metadatafile = f"{self.db_path}/{metadata}"
            elif use_db_path:
                metadatafile = f"{use_db_path}/{metadata}"
            metadata = df_read_genome_metadata (metadatafile, index_name = "peroba_seq_uid")
            # I guess we can assume df.index.name is now "peroba_seq_uid", without duplicates
        self.metadata = metadata
        if self.tree:     self.merge_data_tree ()
        if self.sequence: self.merge_data_sequence ()

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
            if df1: df1 = peroba.utils.df_merge_metadata_by_index (df1, df2)
            else:   df1 = df2
        df1 = df_finalise_metadata (df1)
        self.add_metadata_from_clean_dataframe (df1)


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
