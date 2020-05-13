#!/usr/bin/env python

#import utils 
from utils import *

class DataSeqTree:   ## python3 doesn't use base class object anymore ;)
    """ Basic structure with all available data: sequence, tree, and epi-metadata, with file location information
    All the merge/append between data of the same type should be done elsewhere: this class will drop data that does not
    have correspondence between data, seq, and tree (tree is an exception since we can generate it). 
    This class does check for duplicates, with arbitrary (aka not clever) removal.
    We can have one global DST and several "focused" or local ones, by changing the `prefix` they can be stored in same
    database (aka directory :D )
    """
    db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
    timestamp = ""
    metadata = None
    tree     = None
    tree_leaves = None    # dictionary with {name:node}
    tree_needs_update = False # by default missing tree leaves do not remove samples (since we can reestimate tree) 
    sequences = None      # dictionary with sename:SeqRecord(); peroba.utils and other lowlevel use a lot of lists of SeqRecord()
    f_prefix = None
    f_suffix = {
            "metadata": ".metadata.csv.gz",
            "tree": ".tree.nhx",
            "sequences": ".sequences.fasta.bz2"
            }
                

    def __init__ (self, prefix=None, timestamp = None, treefile=None, raw_metadata_filelist=None, 
            sequence_filelist = None, input_path = None, peroba_db_path = None): 
        """ timestamp overwrites whole file ID (prefix+time-of-day), while prefix changes only the 
        beginning (ToD still updated) 
        """
        if peroba_db_path: self.db_path = peroba_db_path
        else: self.db_path = "/usr/users/QIB_fr005/deolivl/Academic/Quadram/021.ncov/03.perobaDB"
        if prefix is None: self.f_prefix = "peroba."
        else:              self.f_prefix = prefix
        self.reset_timestamp(timestamp) ## uses f_prefix
        if (treefile): self.add_treefile (treefile, use_db_path=input_path)
        if (raw_metadata_filelist): self.add_metadata_from_raw_csv_files (raw_metadata_filelist, use_db_path=input_path)
        if (sequence_filelist): self.add_sequences (sequence_filelist, use_db_path=input_path)

    def reset_timestamp(self, timestamp = None):
        if (timestamp):  self.timestamp = timestamp
        else:            self.timestamp = self.prefix + datetime.datetime.now().strftime("%m%d%H%M%S") 

    def copy_from_DST (self, original = None, prefix = None): # untested
        if original is None: return
        if prefix is None and self.file_prefix is None: 
            self.prefix = original.prefix + datetime.datetime.now().strftime("%H%M%S") + "." 
        self.db_path   = original.db_path
        self.timestamp = original.timestamp[:-10] ## time part only
        self.metadata  = original.metadata.copy(deep=True)
        self.tree = ete3.Tree (original.tree)
        self.tree_leaves = copy.deepcopy (original.tree_leaves)
        self.tree_needs_update = original.tree_needs_update
        self.sequences = copy.deepcopy (original.sequences) # not deep enough, we need to create new SeqRecord to be truly independent
 
    def read_from_database (self, prefix=None, directory=None, subfolder=None): # untested
        if directory is None: directory = self.db_path + "/"
        if subfolder is None: directory = directory
        else:                 directory = f"{directory}/{subfolder}/"
        if prefix is None: # must find the most recent files; assumes default naming
            names = [int(x.replace(self.f_suffix["database"],"")[:-10]) 
                    for x in glob.glob(f"{directory}*" + self.f_suffix["database"])]
            prefix = self.f_prefix + str(sorted(names)[-1]) # highest timestamp => most recent?

        path = directory + prefix # directory + file prefix
        self.add_treefile (treefile=path + self.f_suffix["tree"])
        self.add_metadata_from_clean_dataframe (metadata = path + self.f_suffix["metadata"], use_db_path = False)
        self.add_sequences (sequence_filelist = path + self.f_suffix["sequences"], use_db_path = False)

    def save_to_database (self, prefix=None, directory=None, subfolder=None):
        """ save your initial big table etc, since we may modify them
        """
        if prefix is None: 
            self.reset_timestamp()
            prefix = str (self.timestamp)
        if directory is None: directory = self.db_path + "/"
        if subfolder is None: directory = directory
        else:                 directory = f"{directory}/{subfolder}/"
        pathlib.Path(directory).mkdir(parents=True, exist_ok=True) # python 3.5+
        path = directory + prefix # directory + file prefix
        self.metadata.to_csv (path + self.f_suffix["metadata"])
        self.tree.write(format=1, outfile=path + self.f_suffix["tree"])
        with bz2.open(path + self.f_suffix["sequences"],"wb") as fw: 
            for name, rec in self.sequences.items():
                seq = str(rec.seq)
                fw.write(str(f">{name}\n{seq}\n").encode())

    # TODO: decide if reestimate tree;
    def add_treefile (self, treefile=None, use_db_path=False, merge_now = False): 
        if treefile is None or self.tree is not None: return
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
        print_redblack ("LOG:: leaves in treefile:", len(self.tree_leaves))

        if merge_now:
            self.merge_data_tree ()
            self.merge_sequence_tree ()

    def add_metadata_from_clean_dataframe (self, metadata = None, use_db_path = False, merge_now = False):
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
            else: 
                metadatafile = metadata
            metadata = df_read_genome_metadata (metadatafile, index_name = "peroba_seq_uid")
            # I guess we can assume df.index.name is now "peroba_seq_uid", without duplicates
        #self.metadata = metadata.iloc[::5, :]  ## DEBUG ONLY
        self.metadata = metadata
        print_redblack ("LOG:: final metadata shape (rows x cols)", self.metadata.shape)
        if merge_now:
            self.merge_data_tree ()
            self.merge_data_sequence ()

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

    def add_sequences (self, sequence_filelist = None, use_db_path = False, merge_now = False):
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

        print_redblack ("LOG:: number of valid sequences:", len(self.sequences))
        if merge_now:
            self.merge_data_sequence ()
            self.merge_sequence_tree ()

    def merge_DST (self):
            self.merge_data_sequence ()
            self.merge_sequence_tree ()
            self.merge_data_tree ()

    def merge_data_sequence (self): # keep only intersection
        if self.sequences is None or self.metadata is None:
            return
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_seqs = [x for x in self.sequences.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove sequences absent from the metadata
        idx_sequence_only = set (id_seqs) - set (id_meta) # names unique to sequence, not found in metadata
        idx_matched = set(id_seqs) & set(id_meta)
        df_matched = self.metadata.loc[idx_matched,] ## these rows have well-behaved sequences, with same name 

        unique_id = "central_sample_id" # this is an alternative, unique ID, that is always (?) found within seq names
        _this_is_a_test_ = False  ## problem happens with "raw" elan.consensus.fasta but not with (1 week older) phylogenetics/alignment
        if _this_is_a_test_ and unique_id in self.metadata.columns: # somes thing is wrong if this column is missing from COGUK metadata
            idx_metadata_only = set (id_meta) - set (id_seqs) # names unique to metadata, not found in sequence
            ## dataframe with unique_id of rows without sequences: (notice the list [] of one element, o.w. returns Series
            df_u = self.metadata.loc[idx_metadata_only,[unique_id]] 
            found_indices = []
            # map current seqname to a string likely to have central_sample_id
            description_leaves = {sn:self.sequences[sn].description.replace("EPI_ISL","").replace("COG-UK/","")[:50] for sn in idx_sequence_only} 

            print ("DEBUG:: start test, length:", df_u.shape) # this block can be safely removed (bugged), loop below is slow
            desc_dict = {self.sequences[sn].description.replace("EPI_ISL","").replace("COG-UK/","")[:50]:sn for sn in idx_sequence_only} 
            df1 = pd.DataFrame.from_dict(desc_dict, orient='index', columns=["description"])
            df_u = df_u[df_u[unique_id].isin(df1.index)] ## FIXME: length zero (none found)
            print ("DEBUG:: end test, length:", df_u.shape)

            print ("DEBUG:: entering loop")
            for s_id, s_desc in description_leaves.items(): ## try to find sequence through "central_sample_id", a short code in metadata
                for index, row in df_u.iterrows(): # loop over rows of unmatched metadata
                    if index not in found_indices and str(row[unique_id]) in s_desc: # unique_id was found within the sequence long name
                        self.sequences[index] = self.sequences[s_id]
                        self.sequences[index].id = index
                        #SeqRecord(Seq.Seq(str(xq.seq),Alphabet.IUPAC.ambiguous_dna),id=str(index),description=xq.description)
                        # update dataframes to speed up next iteration
                        found_indices.append(index)
                        if len(found_indices)%10 == 0: print (".", end="")
                        break # next seqname
            print ("DEBUG:: left loop")
            df_matched = df_matched.append(self.metadata.loc[found_indices,]) ## these rows have well-behaved sequences, with same name 

        # sequences found will be new SeqRecords, old ones can be removed
        for seqname in idx_sequence_only:
            del self.sequences[seqname] # delete this dict element
        self.metadata = df_matched
        print_redblack ("LOG:: new dataframe size, after removing missing sequences:", self.metadata.shape)
        # remove table rows without sequence information
        samples_found_in_sequences = self.metadata.index.isin([x for x in self.sequences.keys()])
        self.metadata = self.metadata[samples_found_in_sequences]
        print_redblack ("LOG:: new table size is (after excluding missing sequences):", sum(samples_found_in_sequences))

    def merge_sequence_tree (self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.sequences is None:
            return
        id_seqs = [x for x in self.sequences.keys()]   # copy needed since we may delete elements 
        id_tree = [x for x in self.tree_leaves.keys()] # copy needed since I'll delete elements (modifying dict)
        ## remove leaves without a sequence
        leaves_to_delete = set (id_tree) - set (id_seqs)
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        print_redblack ("LOG:: number of leaves to be deleted (no sequence for it):", len(leaves_to_delete))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## sequences not in the tree
        id_tree = [x for x in self.tree_leaves.keys()] 
        sequences_to_remove = set (id_seqs) - set (id_tree)
        if len(sequences_to_remove): # or we remove sequences or we infer new tree (later)
            print_redblack ("LOG:: number of sequences that would be removed (not in tree):", len(sequences_to_remove))
            if use_tree_to_remove_samples: 
                for seq in sequences_to_remove:
                    del self.sequences[seq]
            else:
                self.tree_needs_update = True  

    def merge_data_tree (self, use_tree_to_remove_samples = False): # keep only intersection
        if self.tree is None or self.metadata is None:
            return
        ## remove leaves not in metadata
        id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
        id_tree = [x for x in self.tree_leaves.keys()]   # copy needed since I'll delete elements (modifying dict)
        leaves_to_delete = set(id_tree) - set(id_meta) # warning: loops can be quite slow
        for l in leaves_to_delete:
            del self.tree_leaves[l] # delete this dict element
        print_redblack ("LOG:: number of leaves to be deleted (not in metadata):", len(leaves_to_delete))
        if len(leaves_to_delete): # some leaves were removed
            self.tree.prune([node for node in self.tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        ## remove table rows absent from tree or just mark tree as incomplete
        samples_found_in_tree = self.metadata.index.isin([x for x in self.tree_leaves.keys()])
        if not all(samples_found_in_tree): # or we purge metadata rows or we infer new tree (later)
            print_redblack ("LOG:: new table size would be (after excluding missing leaves):", sum(samples_found_in_tree))
            if use_tree_to_remove_samples:  self.metadata = self.metadata[samples_found_in_tree]
            else:   self.tree_needs_update = True 
    
    def add_local_data_sequence (self, metadata_file=None, sequence_file=None, use_db_path = False, 
        priority_label="peroba_255", initial_priority = 127): 
        """
        Local (NORW) sequences are named NORW-EXXXX, which is central_sample_id. We match those to COGUK when possible
        This adds all our local (stablished) sequences to the database, i.e. not only the ones to query
        The priority column will have a value from 0 (sequence can be safely excluded) to 255 (sequence must remain);
        therefore local sequences must start with at least 127
        TODO: create a df["week"] s.t. increments every time we add new seqs
        """
        if priority_label is None: priority_label = "peroba_255"
        if initial_priority < 127: initial_priority = 127    
        if initial_priority > 255: initial_priority = 255 # no crime to have larger numbers, some code might truncateit ...
        if use_db_path is True: dbpath = f"{self.db_path}/"
        elif use_db_path:       dbpath = f"{use_db_path}/"
        else:                   dbpath = ""
        filepath = dbpath + sequence_file
        if   "bz2" in sequence_file: zipmode = "bz2"
        elif "gz"  in sequence_file: zipmode = "gz"
        else:                        zipmode = None
        ## read local csv, which may contain oinfo on rejected samples (e.g. as 2020.05.12 there are 
        ## 452 rows but only 302 sequences, since 150 were not sequenced yet or QC rejected)
        norwmeta = df_read_genome_metadata (metadata_file, index_name = "central_sample_id")
        matched = self.metadata[ self.metadata["central_sample_id"].isin(norwmeta.index) ]
        
        ## number of sequences is smaller than metadata; named by central_sample_id
        seq_matrix = read_fasta (filepath, zip = zipmode, check_name = False) 
        norwseqnames = [x.id for x in seq_matrix]
        
        ## temporarily use central_sample_id as index, so that we can merge_by_index
        matched.reset_index(inplace=True) ## downgrades current index to a regular column
        matched.set_index ("central_sample_id", append=True, drop = False, inplace = True) # append creates a new col w/ index
        norwmeta.dropna (axis=1, how='all', inplace=True) # currently, useless columns
        matched.dropna  (axis=1, how='all', inplace=True) # with NA only 

        ## rows of local metadata with sequence, and merge with corresponding elements from global metadata 
        norwmeta = norwmeta[ norwmeta.index.isin(norwseqnames) ]
        norwmeta = df_merge_metadata_by_index (norwmeta, matched) 

        ## revert index to same as global metadata ("peroba_seq_uid" usually)
        norwmeta.reset_index(drop=True, inplace=True) ## drop index completely, not even becomes a column
        norwmeta.set_index (self.metadata.index.names, drop = False, inplace = True)
        
        ## add missing (global) sequence_name with local names; o.w. change local names to follow COGUK
        norwmeta['sequence_name'].fillna(norwmeta.index.to_series(), inplace=True)
        for x in seq_matrix:
            x.id = norwmeta[norwmeta["central_sample_id"] == x.id].index.values[0]
        norwsequences = {x.id:x for x in seq_matrix} 
        self.sequences.update(norwsequences) # duplicates are replaced by local version 

        if priority_label not in list(norwmeta.columns):
            norwmeta[priority_label] = initial_priority
        else: ## some will have values from previous iterations, here we respect their priority
            norwmeta[priority_label].fillna(int(initial_priority), inplace=True)

        self.metadata = df_merge_metadata_by_index (norwmeta, self.metadata) 

    def align_sequences (self, reference_genome=None, trim=None, n_threads = 8):
        if reference_genome is None: reference_genome = f"{self.db_path}/data/MN908947.3.fas"
        align_list = minimap2_align_seqs([x for x in self.sequences.values()], n_threads = n_threads, reference_path = reference_genome)
        if trim is True: trim = [265, 29675]
        if isinstance (trim, list): 
            for x in align_list:
                x.seq = x.seq[trim[0]:trim[1]]
        self.sequences = {x.id:x for x in align_list}

    def rough_tree (self, n_threads = 8):
        tmpfile = "/tmp/roughttree.aln"
        aln = snpsites_from_alignment ([x for x in self.sequences.values()], n_threads = n_threads, outfile=tmpfile)
        tree = rapidnj_from_alignment (infile = tmpfile, n_threads = n_threads)
        os.system("rm -f " + tmpfile)
        self.tree = tree
        self.tree_needs_update = False

    def subsample (self): ## FIXME: just copy-paste so far
        # in order of resolution (34 lineages, 54 special, 479 uk, 1600 phylotypes):
        clade_cols = ["lineage", "special_lineage", "uk_lineage", "phylotype"]
        df = df.sort_values(by=['lineage_support', 'days_since_Dec19'], ascending=[False, True])
        df1 = df.groupby(["lineage", "uk_lineage"]).head(1) # for each group combination, pick first (highest support as we defined above)
        df1 = df.groupby(clade_cols[0]).head(4) 
        df2 = df.groupby(clade_cols[1]).head(3) 
        df3 = df.groupby(clade_cols[2]).filter(lambda x: len(x) > 3).head(1) # only large uk_lineages
        df4 = df[(~df["submission_org_code"].str.contains("NORW")) & (df["adm2"].str.contains("NORFOLK"))]
        df5 = df[df["uk_lineage"] == "x"].head(2); ## more sampling for those not in known lineages 
        #df4 = df[df["submission_org_code"].str.contains("NORW")] ## all
        #df5 = df[df["adm2"].str.contains("NORFOLK")] ## all from NORFOLK and SUFFOLK
        df0 = pd.concat([df1, df2, df3, df4, df5])
        df0 = df0.groupby(df0.index).aggregate("first"); # index is "seq_id"
        #df0['pango_id'] = df0.apply(lambda a: create_pango_id_from_df(a), axis=1) ## recreate using new index
    
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
