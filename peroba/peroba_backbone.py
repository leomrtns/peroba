from peroba import common ## matplotlib trick before anything
from peroba.utils import *
import peroba.regression as ml 
import logging, ete3, argparse, treeswift
import pkg_resources, gc  # garbage collector
import numpy as np, pandas as pd

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_backbone %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

current_working_dir = os.getcwd()
# TODO: we can create a "priority" column based on geography or phylogeny, to sort/preference

# order of preference for samples (True => smaller values preferred)
prefsort = [
    ['peroba_freq_acgt', False],
    [ 'peroba_freq_n', True],
    [ "PCR Ct value", True],
    [ "Coverage (X)", False],
    [ "submission_org_code", True],  # just to add preference for any value instead of nan (i.e. COGUK is preferred)
    [ "adm1", True], # adms and submission_org_code will be catergories, with first ones preferred
    [ "adm2", True],
    [ "country", True],
    [ 'lineage_support', False],
    [ 'collection_date', False] 
    ]

# keys will be categorical columns, values are in order of preference
favourite_cats = { 
        "adm1": ["UK-ENG"],
        "adm2": ["Norfolk", "Suffolk", "Cambridgeshire"],
        "submission_org_code": ["NORW", "CAMB", "SANG", "PHEC", "LOND"],
        "country": ["UK"]
        }

# NaNs will be replaced by "" in these columns (so that we can groupby)
fill_missing_cols = ["lineage", "acc_lineage", "del_lineage", "submission_org_code", 
        "submission_org", "adm0", "adm1", "adm2","uk_lineage","phylotype", "country"]

# cast type as numeric, to allow comparison/sorting
numeric_cols = ['peroba_freq_acgt', 'peroba_freq_n', "PCR Ct value", "Coverage (X)",'lineage_support']

class PerobaBackbone:
    g_csv = None
    g_seq = None # dictionary
    g_snp = None # must match local, therefore calculated first 
    g_names = None # keep a list of current backbone (instead of pruning csv and seq)
    l_csv = None
    l_seq = None # dictionary 
    l_snp = None
    trees = None # these are treeswift trees that will be modified, pruned
    unaligned = None # used only for replacing local seqs by global if these have higher quality
    usrtrees = None ## these are user-defined trees that we'll return enriched (minimal pruning)
    extended_mode = 0
    fast_mode_seqs = False
    #usrseqs = dict() 
    # subset of columns from that may be useful (drop others)
    cols = ["sequence_name", "central_sample_id", "submission_org_code", "submission_org", "collection_datetime", 
            "adm0", "adm1", "adm2", "acc_lineage", "del_lineage",  
            "country" , "cov_id", "sequencing_org", "sequencing_org_code", "sequencing_submission_date",
            "lineage", "lineage_support", "special_lineage", "uk_lineage", "phylotype",
            "peroba_freq_acgt", "peroba_freq_n", "peroba_seq_uid", "source_age", "source_sex", ## until here  from global, below is local
            "adm2_private", "Repeat Sample ID", "icu_admission", "PCR Ct value", "No. Reads", "Mapped Reads", 
            "No. Bases (Mb)", "Coverage (X)", "Average read length", "Basic QC", "High Quality QC", "Missing bases (N)",
            "Consensus SNPs"] 
    cols = ["sequence_name", "central_sample_id", "submission_org_code", "collection_date", 
            "adm0", "adm1", "adm2", "acc_lineage", "del_lineage",  "country" , "lineage", 
            "lineage_support", "uk_lineage", "phylotype", "peroba_freq_acgt", "peroba_freq_n", "peroba_seq_uid"] 
    sort_cols = None 
    sort_ascend = None

    def __init__ (self, peroba_db, global_level, is_fast = False): # not split between global and local yet; peroba_db = [metadata, sequences, tree, unaligned]
        self.g_csv = peroba_db[0]  # formatted metadata
        self.g_seq = {x.id:x for x in peroba_db[1]} # dictionary
        self.trees = [peroba_db[2]]  ## trees is a list; peroba_db[2] is treeswift already (assuming no dup names)
        if len(peroba_db[3]) > 0: # only if user wants to compare with local to keep higher quality
            self.unaligned = {x.id:x for x in peroba_db[3]} # dictionary

        cols = [x for x in self.cols if x in peroba_db[0].columns]
        self.g_csv = self.g_csv[cols] # remove other columns
        
        logger.info("Imported %s rows from database", str(self.g_csv.shape[0]))
        self.fast_mode_seqs = is_fast;
        if (self.fast_mode_seqs is True):
            logger.info ("Fast mode: will update only sequences given as 'local'; NORW samples will be treated as regular COGUK ones")
        else: 
            self.fast_mode_seqs = False # anything other than True is False
            logger.info ("Slow mode: will update all NORW sequences, treating them separately as other COGUK ones")

        self.extended_mode = global_level
        if self.extended_mode == 0:
            logger.info("Assuming all sequences are from NORW/COGUK (subsampling will be UK-centric)")
        elif self.extended_mode == 1:
            logger.info("Extended (global) mode: will use more GISAID data, with settings similar to the classic (COGUK) search")
        else:
            self.extended_mode = 2
            logger.info("Extended (global) mode: will use more GISAID data and perform a broader search")

    
    def add_local_data_and_sequences (self, csv=None, sequence=None):
        if not sequence:
            logger.warning("Nothing to merge without sequences")
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
                else: ## will be added through a new_rows dataframe
                    s_new.append (seq.id)
                    #self.g_csv.loc[str(seq.id)] = pd.Series({'sequence_name':seq.id, 'central_sample_id':seq.id,'submission_org_code':"NORW", 'submission_org':"Norwich"})
            else: # sequence has long, official name
                if "NORW" in seq.id: ## we only consider replacing local seqs, otherwise database might have newer  
                    s_long.append(seq.id)
                    seqs_in_global.append(seq.id)
                else:
                    seq.id = None
        gc.collect()

        if len(s_new) > 0:
            logger.warning("Sequences not found in global database will be added by hand:\n%s\n", "   ".join(s_new))
            new_rows = pd.DataFrame({'peroba_seq_uid':s_new, 'sequence_name': s_new, 'central_sample_id':s_new,
                'submission_org_code':"NORW" ,'submission_org':"Norwich"})
            new_rows = new_rows.groupby("peroba_seq_uid").aggregate("first") # duplicate sequences --> this makes peroba_seq_uid the index 
            self.g_csv = self.g_csv.combine_first(new_rows)
        if (self.fast_mode_seqs is True):
            if len (s_new) > 0:
                self.fast_mode_seqs = s_new
            else:
                logger.warning ("Fast mode selected but no new local sequences found (all already found on database); Changing back to slow mode")
                self.fast_mode_seqs = False

        logger.info("%s long-named sequences and %s short-named sequences found on database. %s new sequences (not in database).",
                str(len(s_long)), str(len(s_short)), str(len(s_new))) # obvious duplicates were checked by initial_quality_control()
        
        # if duplicates with both long and short names, by now both have long name
        logger.info("Checking sequences where both short- and long-named versions were input (by now all these have long version)")
        sequence, l_qual = common.remove_duplicated_sequences (sequence) # return a dict of seqs and one of qualities 
        # self.unaligned has unaligned global sequences, and is used only if user doesn't trust global seqs (and wants
        # to check back). In which case we must neglect the current aligned and align it again 
        if self.unaligned is not None and len(seqs_in_global) > 0: 
            logger.info("Replacing local and global sequences by one with better quality (more non-N bases)")
            self.unaligned, sequence, to_replace = common.sequence_dict_pair_with_better_quality (self.unaligned, sequence, q2=l_qual, matched=seqs_in_global)
            for seqname in to_replace[1]: # list of local seq names with better quality
                del self.g_seq[seqname] # will align and use local version instead 
            seqs_in_global = [x for x in seqs_in_global if x not in to_replace[1]] # must realign even if present in global
            self.unaligned.clear()
            self.unaligned = None
        
        if len(seqs_in_global) > 0:
            seq_list = [x for x in sequence.values() if x.id not in seqs_in_global]
        else:
            seq_list = [x for x in sequence.values()]

        if len(seq_list) > 0: # we have new sequences, not yet in COGUK 
            sqn = [x.id for x in seq_list]
            logger.info("Will update frequency info and align %s sequences", str(len(sqn)))
            sqn = self.g_csv["sequence_name"].isin(sqn)
            self.g_csv.loc[sqn, "peroba_freq_acgt"] = np.nan # recalculate
            self.g_csv.loc[sqn, "peroba_freq_n"] = np.nan
            self.g_csv, sequence = common.add_sequence_counts_to_metadata (self.g_csv, sequence, from_scratch=False) # seqname must be in index

            # merge sequences (align local first, since global are already aligned)
            ref_seq = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas")  
            aln = common.align_sequences_in_blocks (seq_list, reference_file = ref_seq, seqs_per_block = 1000)
            logger.info("Finished aligning %s local sequences", str(len(aln)))
            self.g_seq.update({x.id:x for x in aln})
        else:
            logger.info("No new sequences to be updated or aligned")

        if csv is not None:
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
        return

    def trim_sequences (self, trim = True):  
        if trim is True: trim = [265, 29675]
        if trim[0] < 1: trim[0] = 1
        if trim[1] > 29902: trim[1] = 29902 # genome size is 29903
        logger.info("Trimming genomes from site %s to site %s",str(trim[0]), str(trim[1]))
        if self.g_seq is not None:  # this should always work
            for x in self.g_seq.values():
                x.seq = x.seq[trim[0]:trim[1]]
        if self.l_seq is not None:  # usually this is None since we trim _before_ splitting 
            for x in self.l_seq.values():
                x.seq = x.seq[trim[0]:trim[1]]
    
    def sort_categories (self):
        df = self.g_csv.copy()
        pref = [x for x in prefsort if x[0] in df.columns] # sort only existing columns
        if self.extended_mode > 0: # not UK sequences; thus do not favour/sort geographically
            pref = [x for x in pref if x not in favourite_cats.values()] # exclude uk centric preferences
        self.sort_cols = [x[0] for x in pref] # follow this order (thus can't use dict...)
        self.sort_ascend = [x[1] for x in pref] # same order as sort_cols
        logger.info("Ordering metadata after categorisation") 
        if "adm2" in df.columns:
            df["adm2"] = df["adm2"].replace(["Unknown Source","Unknown", np.nan],"")
            df["adm2"] = df["adm2"].replace({"Greater London":"Greater_London", "Hertfordshire":"Herefordshire"})

        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col])

        for col in fill_missing_cols:
            if col in df.columns:
                df[col].fillna ("", inplace = True) 

        if self.extended_mode == 0: # favourite categories are NORW-centric 
            for col,favs in favourite_cats.items():
                if col in df.columns:
                    values = [x for x in df[col].unique() if x not in favs + [""]] 
                    values = favs + values  + [""] # make sure favs are first in list, and empty are last
                    df[col] = pd.Categorical(df[col], categories = values, ordered=True)

        df = df.sort_values(by=self.sort_cols, ascending=self.sort_ascend)
        self.g_csv = df

    def finalise_and_split_data_sequences (self, trim = True, strict = False):
        self.trim_sequences(trim=trim)
        self.sort_categories()
        logger.info("Finding SNPs to speed up calculations")
        snp_aln = snpsites_from_alignment ([x for x in self.g_seq.values()], strict=strict)
        logger.info("In total %s SNPs were found (i,e, alignment size)", len(snp_aln[0].seq))
        
        if (self.fast_mode_seqs is False):
            logger.info ("Splitting data into global and local (NORW)")
            self.l_csv = self.g_csv[  self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
            self.g_csv = self.g_csv[ ~self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
        else:
            logger.info ("Splitting data into local (not found on COGUK) and global (COGUK)")
            self.l_csv = self.g_csv[  self.g_csv["sequence_name"].isin(self.fast_mode_seqs) ] # list of names
            self.g_csv = self.g_csv[ ~self.g_csv["sequence_name"].isin(self.fast_mode_seqs) ]

        self.l_seq = {x:y for x,y in self.g_seq.items() if x in self.l_csv["sequence_name"]}
        self.g_seq = {x:y for x,y in self.g_seq.items() if x in self.g_csv["sequence_name"]}
        
        #l_snp = [x for x in snp_aln if x.id in self.l_csv["sequence_name"]]
        #idx = sorted_uncertainty_from_alignment (l_snp, max_freq_n = 0.01)
        #logger.info("Compact representation uses %s SNPs (columns from local seqs with fewer Ns)", len(idx))
        #snp_aln = alignment_from_column_index (snp_aln, idx)
        self.l_snp = {x.id:x for x in snp_aln if x.id in self.l_csv["sequence_name"]}
        self.g_snp = {x.id:x for x in snp_aln if x.id in self.g_csv["sequence_name"]}

        logger.info ("split data into %s global and %s local sequences", str(self.g_csv.shape[0]), str(self.l_csv.shape[0]))
        gc.collect() ## collect garbage

    def remove_seq_tree_based_on_metadata (self, seqnames = None, local = False): 
        if seqnames is not None and local:
            logger.warning("Removing local sequences selected by hand (i.e. not from metadata update); structure may be compromised")
        # default is to remove only global sequences
        if seqnames is None: 
            if local:
                logger.info("removing local sequences after metadata update")
                seqnames = self.l_csv["sequence_name"].tolist()
            else:
                logger.info("removing global sequences after global metadata update")
                seqnames = self.g_csv["sequence_name"].tolist()
        if local:
            newseqs = {x:y for x,y in self.l_seq.items() if x in seqnames}
            newsnps = {x:y for x,y in self.l_snp.items() if x in seqnames}
            seqnames += self.g_csv["sequence_name"].tolist() # add all global to avoid being pruned 
        else:
            newseqs = {x:y for x,y in self.g_seq.items() if x in seqnames}
            newsnps = {x:y for x,y in self.g_snp.items() if x in seqnames}
            seqnames += self.l_csv["sequence_name"].tolist() # add all local to avoid being pruned 
        newtrees = []
        for i,t in enumerate(self.trees): # computing-intensive step 
            remove_leaf = []
            for l in t.traverse_leaves():
                lab = l.get_label()
                if lab not in seqnames:
                    remove_leaf.append(lab)
            if len(remove_leaf) > 0:
                logger.warning("In tree id %s, %s leaves were pruned (absent from metadata or list)", str(i), str(len(remove_leaf)))
                newt = t.extract_tree_without(remove_leaf)
                if newt.num_nodes(internal=False) > 4: ## number of leaves
                    newtrees.append(newt)
            else:
                newtrees.append(t)
        return  newseqs, newsnps, newtrees

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
                logger.info("In tree id %s the number of leaves to be pruned is %s", str(i), len(remove_leaf))
                logger.debug("and they are:\n%s\n", "\n".join(remove_leaf))
                self.trees.append(t.extract_tree_without(remove_leaf))
            else:
                self.trees.append(t)
        # copy of user-defined trees that don't get pruned
        if len(self.trees) > 1:
            self.usrtrees = []
            for t in self.trees[1:]:
                self.usrtrees.append( treeswift.read_tree_newick (str(t)) ) ## copy 

    def remove_duplicates (self, blocks = 4, leaf_size = 500, radius=0.00001):
        logger.info("Removing duplicates (identical sequences) from global data")
        clusters = ml.list_duplicates (self.g_snp, blocks, leaf_size, radius)

        df = pd.DataFrame([[y,i] for i,x in enumerate(clusters) for y in x], columns=["sequence_name","peroba_tmp"])
        df = df.sort_values(by=["peroba_tmp"], ascending=True) # same sequence belongs to several clusters
        df = df.groupby("sequence_name").aggregate("first") # now each sequence belongs only to lowest cluster number
        logger.info("Number of unique sequences: %s", str(len(df["peroba_tmp"].unique())) )

        self.g_csv.reset_index (drop=False, inplace=True) ## merge will destroy index...
        df = self.g_csv.merge(df, on="sequence_name")
        df = df.sort_values(by=self.sort_cols, ascending=self.sort_ascend)

        ## create list of accepted seqnames instead of deleting rows; replace peroba_tmp by peroba_group
        #for c in [x for x in replace_duplicate_cols if x in df.cols]: # most common value amongst identical sequences
        #    df[c] = df.groupby(["peroba_tmp"])[c].transform(pd.Series.mode) # transform() keeps index

        df = df.groupby("peroba_tmp").aggregate("first") # only one sequence from each cluster, following self.order_col preference
        df.set_index ("peroba_seq_uid", drop = True, inplace = True)
        self.g_csv = df
        self.g_seq, self.g_snp, self.trees = self.remove_seq_tree_based_on_metadata()
    
    def reduce_redundancy (self, clade_rule = None):
        if clade_rule is None: # for each lineage level with at least x[0] samples, keep at most x[1]
            clade_rule = [ # rules can be repeated, for different thresholds; some samples fall into several 
                    ["lineage",   10, 50], # only those with >1 samples; then take up to 50
                    ["lineage",   40, 1000], # only those with >20 samples; then take up to 1000 
                    ["adm1",      10, 200],
                    ["adm2",       5, 100],
                    ["uk_lineage", 4, 20], 
                    ["uk_lineage",20, 500], 
                    ["phylotype",  4, 20],
                    ["acc_lineage", 1, 10],
                    ["del_lineage", 1, 10]
                    ]
        missing_rule = [ ## for samples without this information, we keep the top ones 
                ["lineage",    50], 
                ["uk_lineage", 500], 
                ["phylotype",  500], 
                ["acc_lineage", 1000], 
                ["del_lineage", 1000]
                ]
        if self.extended_mode == 1: # more permissive since many non-UK sequences don't have this info
            clade_rule = [[x[0], int(x[1]/2), 4 * x[2]] for x in clade_rule] # zero is fine 
            missing_rule = [[x[0], 16 * x[1]] for x in missing_rule]
        if self.extended_mode == 2: # more permissive since many non-UK sequences don't have this info
            clade_rule = [[x[0], int(x[1]/2), 8 * x[2]] for x in clade_rule] # zero is fine 
            missing_rule = [[x[0], 32 * x[1]] for x in missing_rule]

        df = self.g_csv
        logger.info("Subsampling redundant global samples (i.e. those fom same lineage etc.)") 
        df = df.sort_values(by=self.sort_cols, ascending=self.sort_ascend)
        
        dfcat = None
        for column, rule1, rule2 in clade_rule:
            if column in df.columns:
                # groupby().filter() creates another DF (excluding rows with less than rule1)
                # therefore a second groupby is needed: groupby().head()
                df1 = df.groupby(column).filter(lambda x: len(x.index) > rule1).groupby(column).head(rule2)
                if dfcat is None: dfcat = df1
                else: dfcat = pd.concat([dfcat, df1])
        for column, n_keep in missing_rule: 
            if column in df.columns:
                df1 = df[ df[column] == "" ].head(n_keep)  # undefined (missing) lineages 
                if dfcat is None: dfcat = df1
                else: dfcat = pd.concat([dfcat, df1])

        self.g_csv = dfcat.groupby(dfcat.index).aggregate("first") # many rows will be duplicated
        logger.info("After subsampling, global metadata has %s samples", self.g_csv.shape[0])
        self.g_seq, self.g_snp, self.trees = self.remove_seq_tree_based_on_metadata()

    def remove_low_quality (self, g_acgt = 0.75, l_acgt = 0.5, g_n = 0.1, l_n = 0.3):
        logger.info(f"Remove global sequences with proportion of ACGT less than  {g_acgt} or proportion of N higher than {g_n}")
        self.g_csv = self.g_csv.loc[ (self.g_csv["peroba_freq_acgt"] > g_acgt) & (self.g_csv["peroba_freq_n"] < g_n) ]
        self.g_seq, self.g_snp, self.trees = self.remove_seq_tree_based_on_metadata()

        logger.info(f"Also removing local sequences with proportion of ACGT less than {l_acgt} or proportion of N higher than {l_n}")
        self.l_csv = self.l_csv.loc[ (self.l_csv["peroba_freq_acgt"] > l_acgt) & (self.l_csv["peroba_freq_n"] < l_n) ]
        self.l_seq, self.l_snp, self.trees = self.remove_seq_tree_based_on_metadata(local=True)
        logger.info("After removal, global data has %s samples and local data has %s.", self.g_csv.shape[0], self.l_csv.shape[0])

    def find_neighbours_ball (self, blocks = 1500, leaf_size = 400, dist_blocks = 3, nn = 5, exclude = None):
        logger.info(f"Finding neighbours to local sequences, using a distance of {dist_blocks} to {blocks} segments")
        # search only amongst sequences with lineage (notice that sort_categories() replaces null by "")
        target = self.g_csv.loc[(self.g_csv["lineage"] != ""), "sequence_name"].tolist()
        tmp1 = sum(self.g_csv["lineage"] != "")
        if exclude is not None:
            target = [x for x in target if x not in exclude]
        g_seq = {x:self.g_snp[x] for x in target} 
        neighbours1 = ml.list_r_neighbours (g_seq, self.l_snp, blocks, leaf_size, dist_blocks)
       
        blocks = 2 * blocks; leaf_size = leaf_size/2
        logger.info("Found %s neighbours; now will find %s closest neighbours to %s on %s segments", 
                len(neighbours1), str(nn), "them" if len(neighbours1)>0 else "the queries", str(blocks))
        if len(neighbours1) > 0:
            aln_d = {x:self.g_snp[x] for x in neighbours1} 
        else:
            logger.warning ("No neighbours found within appropriate distance: using query sequences to find their nearest neighbours") 
            aln_d = self.l_snp
        # search amongst all sequences, not only those with lineage info
        neighbours = ml.list_n_neighbours (self.g_snp, aln_d, blocks, leaf_size, nn)
        neighbours = list(set(neighbours + neighbours1))
        gc.collect() ## collect garbage
        logger.info("Found %s neighbours through hashed nearest-neighbours", len(neighbours))
        return neighbours

    def find_neighbours_paf (self, blocks = 2000, leaf_size = 500, n_segments = 1, nn = 10, exclude = None):
        logger.info(f"Finding neighbours to unclassified local sequences, by alignment mapping of {n_segments} segment(s)")
        query = self.l_csv.loc[(self.l_csv["uk_lineage"] == ""), "sequence_name"].tolist()
        if len(query) < 1: 
            logger.info("Skip alignment mapping since all local sequences already have uk_lineage")
            return []
        if self.extended_mode > 0: # any sample with a global lineage
            target = self.g_csv.loc[(self.g_csv["lineage"] != ""),"sequence_name"].tolist()
        else: # only samples with a uk_lineage
            target = self.g_csv.loc[(self.g_csv["uk_lineage"] != ""),"sequence_name"].tolist()
        if exclude is not None:
            target = [x for x in target if x not in exclude]
        l_seq = {x:self.l_snp[x] for x in query}  # only local sequences still without uk_lineage
        g_seq = {x:self.g_snp[x] for x in target} # search only amongst those with uk_lineage
        neighbours1 = ml.list_paf_neighbours (g_seq, l_seq, n_segments = n_segments, n_threads = 2)
       
        logger.info("Found %s neighbours; now will find their %s closest neighbours on %s segments", 
                len(neighbours1), str(nn), str(blocks))
        aln_d = {x:self.g_snp[x] for x in neighbours1}
        neighbours = ml.list_n_neighbours (g_seq, aln_d, blocks, leaf_size, nn)
        neighbours = list(set(neighbours + neighbours1))
        logger.info("Found %s neighbours", len(neighbours))
        return neighbours

    def find_neighbours (self):
        if self.extended_mode == 2:
            n1 = self.find_neighbours_ball(blocks = 2000, leaf_size = 400, dist_blocks = 5, nn = 20) 
            n2 = self.find_neighbours_paf (blocks = 4000, leaf_size = 500, n_segments = 2, nn = 30, exclude = n1) 
        elif self.extended_mode == 1:
            n1 = self.find_neighbours_ball(blocks = 2000, leaf_size = 400, dist_blocks = 4, nn = 20) 
            n2 = self.find_neighbours_paf (blocks = 3000, leaf_size = 500, n_segments = 1, nn = 20, exclude = n1) 
        else: ## local (COGUK) mode
            n1 = self.find_neighbours_ball(blocks = 1500, leaf_size = 400, dist_blocks = 3, nn = 5) 
            n2 = self.find_neighbours_paf (blocks = 2000, leaf_size = 500, n_segments = 1, nn = 10, exclude = n1) 

        if self.fast_mode_seqs is not False: # Add all remaining NORW sequences as _global_ 
            n1 = list(set(n1 + n2))
            n2 = self.g_csv.loc[(self.g_csv["submission_org_code"].str.contains("NORW", na=False)), "sequence_name"].tolist()
        return list(set(n1 + n2))

def save_global_from_seqnames (bb, seqnames, prefix):
    seqs, snps, trees = bb.remove_seq_tree_based_on_metadata (seqnames)
    desc = common.save_sequence_dict_to_file (seqs, prefix + ".aln.xz")
    fname = prefix + ".trees.nhx" 
    with open(fname,"w") as fw:
        for t in trees:
            fw.write(str(t) + "\n")
    logger.info(f"Finished saving global tree(s) to file {fname}")
    
    # format(colour_string(desc).ljust(32, ' ')) # not using ljust anymore
    desc  = "{}\n Aligned sequences selected from global data\n".format(colour_string(desc))

    fname = os.path.basename(fname)
    desc += "{}\n Pruned trees, with leaves only from global data (please refer to 'user' trees below for phylo inference)\n".format(colour_string(fname))
    desc += " First tree is full tree from COGUK, remaining trees are provided from user, if any\n"
    return desc + "\n"

def save_user_trees (bb, seqnames, prefix, add_nj_tree = False):
    if bb.usrtrees is None:
        logger.info ("No user-defined trees, skipping creation of custom data set")
        return "" 
    logger.info ("Saving user-defined trees with added sequences")
    leafnames = set([l.get_label() for t in bb.usrtrees for l in t.traverse_leaves()])
    aln = {x:y for x,y in bb.g_seq.items() if x in leafnames}
    aln.update({x:y for x,y in bb.l_seq.items() if x in leafnames})
    trees = bb.usrtrees

    if len(aln) < len(leafnames):
        remove_leaf = [l.get_label() for t in bb.usrtrees for l in t.traverse_leaves() if l.get_label() not in aln.keys()]
        logger.info("In total %s user-defined leaves do not have a sequence, and will be pruned",str(len(remove_leaf)))
        logger.debug("And they are\n%s\n", "\t".join(remove_leaf))
        trees = [t.extract_tree_without(remove_leaf) for t in bb.usrtrees]

    ## add sequences from current analysis (all local plus global neighbours)
    aln.update({x:y for x,y in bb.g_seq.items() if x in seqnames})
    aln.update({x:y for x,y in bb.l_seq.items()})
    desc =  common.save_sequence_dict_to_file (aln, prefix + ".aln.xz")
    desc  = "{}\n alignment with sequences from all user-defined trees (that I have access to), as well as all \n".format(colour_string(desc))
    desc += " global and local sequences.\n"

    logger.info ("In total, %s sequences are part of the user-defined data set (user trees plus sequences found here)",
            str(len(aln)))

    fname = prefix + ".trees.nhx" 
    with open(fname,"w") as fw:
        treestring = None
        if add_nj_tree:
            logger.info ("Estimating NJ tree with all user-defined sequences (will be the first tree in file)")
            snp_aln = snpsites_from_alignment ([x for x in aln.values()])
            treestring = rapidnj_from_alignment (snp_aln, n_threads =12)
            fw.write(treestring + "\n")
        for t in trees:
            fw.write(str(t) + "\n")
    logger.info(f"Finished saving user-defined trees to file {fname}")

    fname = os.path.basename(fname)
    desc += "{}\n User-defined trees with subset of leaves for which I've found sequences\n".format(colour_string(fname))
    if (add_nj_tree):
        desc += " where the first tree is a NJ estimate from all 'user' sequences above\n"
    return desc + "\n"

def save_all_sequences (bb, seqnames, prefix, add_nj_tree = True):
    logger.info ("Saving data set with all local and selected global sequences")
    ## add sequences from current analysis (all local plus global neighbours)
    aln = {x:y for x,y in bb.g_seq.items() if x in seqnames}
    aln.update({x:y for x,y in bb.l_seq.items()})
    desc =  common.save_sequence_dict_to_file (aln, prefix + ".aln.xz")
    desc = "{}\n alignment with (all) local and (selected) global sequences. This is main output of program\n".format(colour_string(desc))
    logger.info ("Total of %s sequences will form the local+global data set", str(len(aln)))

    if add_nj_tree:
        logger.info ("Estimating NJ tree with local+global sequences")
        snp_aln  = [x for x in bb.g_snp.values() if x.id in seqnames] # list, not dict
        snp_aln += [x for x in bb.l_snp.values()]
        treestring = rapidnj_from_alignment (snp_aln, n_threads =12) # list input, string output
        fname = prefix + ".trees.nhx" 
        with open(fname,"w") as fw:
            fw.write(treestring + "\n")
        logger.info(f"Finished saving local+global tree to {fname}")
        fname = os.path.basename(fname)
        desc += "{}\n NJ estimate of all global and local sequences (useful for phylo inference)\n".format(colour_string(fname))
    else:
        desc += " (NJ estimation was not requested)\n"
    return desc + "\n"

def read_peroba_database (f_prefix, trust_global_sequences = False): 
    if f_prefix[-1] == ".": f_prefix = f_prefix[:-1] ## both `perobaDB.0621` and `perobaDB.0621.` are valid
    fname = f_prefix + common.suffix["metadata"]
    logger.info(f"Reading database metadata from \'{fname}\'")
    metadata = pd.read_csv (fname, compression="infer", index_col="peroba_seq_uid", dtype="unicode") 
    metadata = common.df_finalise_metadata (metadata) 

    fname = f_prefix + common.suffix["tree"]
    logger.info(f"Reading database tree from \'{fname}\'")
    treestring = open(fname).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = treeswift.read_tree_newick (treestring) 

    fname = f_prefix + common.suffix["alignment"]
    logger.info(f"Reading database alignment from \'{fname}\'")
    sequences = common.read_fasta (fname, check_name = False)

    unaligned = []
    if trust_global_sequences:
        logger.info(f"Will assume global sequences are 'better' than local when duplicates exist")
    else: 
        fname = f_prefix + common.suffix["sequences"]
        logger.info(f"Reading database unaligned sequences from \'{fname}\'")
        unaligned = common.read_fasta (fname, check_name = False)
    

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequences, tree, unaligned]

def sequences_initial_quality_control (sequences, min_length = 20000, max_freq_N = 0.5):
    pass_seqs = []
    pass_qual = []
    fail_names = [] 
    ## remove very low quality, and store quality values (some may be duplicates)
    for x in sequences:
        seq = str(x.seq)
        n1 = common.calc_freq_N_from_string (seq)
        l = len(seq)
        if n1 < max_freq_N and l > min_length:
            pass_seqs.append(x)
            pass_qual.append((1. - n1) * l) # number of non-N bases
        else:
            fail_names.append([x.id, "{:.2f}".format(n1*100), l])
    if len(fail_names) > 0:
        logger.warning ("From local sequences, %s passed and %s failed quality control", len(pass_seqs), len(fail_names))
        fail_detail = [f"Sequence {x[0]} has {x[1]}% of Ns and length {x[2]}" for x in fail_names]
        logger.debug("Sequences that failed basic quality control:\n%s\n","\n".join(fail_detail))
    else:
        logger.info ("Quality control of local sequences: all %s passed", len(pass_seqs))

    uniq_seqs = {}
    uniq_qual = {}
    duplicates = []
    for s,q in zip (pass_seqs, pass_qual):
        if s.id in uniq_seqs.keys(): # sequence is a duplicate
            duplicates.append(s.id)
            if uniq_qual[s.id] < q:
                uniq_qual[s.id] = q
                uniq_seqs[s.id] = s
        else:
            uniq_qual[s.id] = q
            uniq_seqs[s.id] = s
    duplicates = list(set(duplicates))
    if len(duplicates) > 0:
        logger.warning ("From local sequences, %s appear with same name more than once (the longer and with less Ns was chosen).", len(duplicates))
        logger.debug("List of duplicated names:\n%s\n","\t".join(duplicates))
    else:
        logger.info ("After checking, all local sequences have a unique name")

    return [x for x in uniq_seqs.values()]

def main_generate_backbone_dataset (database, csv, sequences, trees, prefix, global_level, fast_mode):
    # initialisation
    bb = PerobaBackbone (database, global_level, fast_mode)
    sequences = sequences_initial_quality_control (sequences) # very bad sequences break MAFFT
    bb.add_local_data_and_sequences (csv, sequences)
    bb.finalise_and_split_data_sequences()
    bb.add_trees (trees)
    bb.remove_low_quality()
    bb.remove_duplicates()
    bb.reduce_redundancy() ## add step to remove too distant 
    # all methods come here
    neighbours = bb.find_neighbours()
    # for each method we can use chosen neighbours to save global only or both
    if bb.fast_mode_seqs is False:
        all_seqs_fname = "norw-coguk"
        loc_seqs_fname = "norw.aln.xz"
    else:
        all_seqs_fname = "local-coguk"
        loc_seqs_fname = "local.aln.xz"
    description  = save_global_from_seqnames (bb, neighbours, prefix + "coguk")
    description += save_user_trees (bb, neighbours, prefix + "user", add_nj_tree = True)
    description += save_all_sequences (bb, neighbours, prefix + all_seqs_fname, add_nj_tree = True)
    # finally, save all NORW sequences (or local, if fast_mode)
    desc =  common.save_sequence_dict_to_file (bb.l_seq, prefix + loc_seqs_fname)
    description += "{}\n alignment with all local sequences only\n".format(colour_string(desc))

    print ("Finished. The output files are described below, where 'global' means COGUK and GISAID data which were ")
    print ("{} generated in NORW. Those, together with the extra sequences are being called 'local'.".format(colour_string("not", "red")))
    print (f"Files produced:\n{description}")
    
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
    It's recommended that you also include local sequences, even without CSV metadata. You can furthermore add a newick file with extra 
    trees (the tree from previous run, for instance).
    """, 
    usage='''peroba_backbone <perobaDB> [options]''')

    parser.add_argument('perobaDB')
    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
        help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-i', '--input', action="store", help="Directory where perobaDB files are. Default: working directory")
    parser.add_argument('-c', '--csv', metavar='csv', help="csv table with metadata from NORW")
    parser.add_argument('-s', '--sequences', metavar='fasta', nargs='+', help="extra files with local sequences (from NORW)")
    parser.add_argument('-t', '--trees', metavar='', help="file with (user-defined) trees in newick format to help produce backbone")
    parser.add_argument('-o', '--output', action="store", help="Output database directory. Default: working directory")
    parser.add_argument('-g', '--global_level', metavar='[0,1,2]', type=int, default=0, 
        help="how broad the search should be (default=0 wich means local (COGUK) new samples only)")
    parser.add_argument('-f', '--fast', default=False, action='store_true', 
        help="Fast mode (known NORW samples are added to backbone and not to query)")
    parser.add_argument('-r', '--trust', default=False, action='store_true', 
        help="Trust global sequences, skipping quality comparison for matches")

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
    database = read_peroba_database (os.path.join(input_d, args.perobaDB), args.trust) # something like "my_folder/perobaDB.0515"

    csv = None
    if (args.csv):
        fname = os.path.join(current_working_dir, args.csv)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.csv)
        if not os.path.exists(fname):
            logger.warning (f"Could not find local CSV file {args.csv}; Will proceed without it")
        else:
            logger.info("Reading CSV file with metadata from NORW")
            csv = common.df_read_genome_metadata (fname, index_name = "central_sample_id")
            csv = common.df_finalise_metadata (csv)

    sequences = None
    if (args.sequences):
        logger.info("Reading fasta files with sequences from NORW")
        if isinstance (args.sequences, str):
            seqfiles = [args.sequences]
        else:
            seqfiles = args.sequences
        sequences = []
        for f in seqfiles:
            fname = os.path.join(input_d, f)
            if not os.path.exists(fname):
                fname = os.path.join(current_working_dir, f)
            if not os.path.exists(fname):
                logger.warning (f"Could not find sequence file {f}; Will proceed without it")
            else:
                s = common.read_fasta (fname, check_name = False)
                sequences += s

    trees = None
    if (args.trees):
        fname = os.path.join(current_working_dir, args.trees)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.trees)
        if not os.path.exists(fname):
            logger.warning (f"Could not find tree file {args.trees}; Will proceed without it")
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

    main_generate_backbone_dataset (database, csv, sequences, trees, prefix, args.global_level, args.fast)


if __name__ == '__main__':
    main()
