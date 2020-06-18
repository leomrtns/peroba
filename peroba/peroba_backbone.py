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
    l_csv = None
    l_seq = None # dictionary 
    l_snp = None
    trees = None ## these are treeswift trees that will be modified, pruned
    usrtrees = None ## these are user-defined trees that we'll return enriched (minimal pruning)
    usrseqs = dict() ## all discarded sequences (which may be on user tree)
    # subset of columns from that may be useful (drop others)
    cols = ["sequence_name", "central_sample_id", "submission_org_code", "submission_org", "collection_datetime", 
            "adm0", "adm1", "adm2", "acc_lineage", "del_lineage",  
            "country" , "cov_id", "sequencing_org", "sequencing_org_code", "sequencing_submission_date",
            "lineage", "lineage_support", "special_lineage", "uk_lineage", "phylotype",
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
        
        logger.info("Imported %s rows from database", str(self.g_csv.shape[0]))
    
    def add_local_data_and_sequences (self, csv=None, sequence=None, replace = False):
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
                if "NORW" in seq.id: ## we only consider replacing local seqs, otherwise database migth have newer  
                    s_long.append(seq.id)
                    seqs_in_global.append(seq.id)
                else:
                    seq.id = None

        if len(s_new) > 0:
            logger.warning("Sequences not found in global database will be added by hand:\n%s\n", "   ".join(s_new))
            new_rows = pd.DataFrame({'peroba_seq_uid':s_new, 'sequence_name': s_new, 'central_sample_id':s_new,
                'submission_org_code':"NORW" ,'submission_org':"Norwich"})
            new_rows = new_rows.groupby("peroba_seq_uid").aggregate("first") # duplicate sequences --> this makes peroba_seq_uid the index 
            self.g_csv = self.g_csv.combine_first(new_rows)

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
        sqn = self.g_csv["sequence_name"].isin(sqn)
        self.g_csv.loc[sqn, "peroba_freq_acgt"] = np.nan # recalculate
        self.g_csv.loc[sqn, "peroba_freq_n"] = np.nan
        self.g_csv, sequence = common.add_sequence_counts_to_metadata (self.g_csv, sequence, from_scratch=False) # seqname must be in index
        ## FIXME :: very low quality sequences break mafft (>90% Ns for example)

        # merge sequences (align local first, since global are already aligned)
        ref_seq = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas")  
        aln = common.align_sequences_in_blocks (seq_list, reference_file = ref_seq, seqs_per_block = 1000)
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

        logger.info ("Splitting data into global and local (NORW)")
        self.l_csv = self.g_csv[  self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
        self.g_csv = self.g_csv[ ~self.g_csv["submission_org_code"].str.contains("NORW", na=False) ]
        self.l_seq = {x:y for x,y in self.g_seq.items() if x in self.l_csv["sequence_name"]}
        self.g_seq = {x:y for x,y in self.g_seq.items() if x in self.g_csv["sequence_name"]}
        self.l_snp = {x.id:x for x in snp_aln if x.id in self.l_csv["sequence_name"]}
        self.g_snp = {x.id:x for x in snp_aln if x.id in self.g_csv["sequence_name"]}
        logger.info ("split data into %s global and %s local sequences", str(self.g_csv.shape[0]), str(self.l_csv.shape[0]))

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
            self.usrseqs.update({x:y for x,y in self.l_seq.items() if x not in seqnames})
            seqnames += self.g_csv["sequence_name"].tolist() # add all global to avoid being pruned 
        else:
            newseqs = {x:y for x,y in self.g_seq.items() if x in seqnames}
            newsnps = {x:y for x,y in self.g_snp.items() if x in seqnames}
            self.usrseqs.update({x:y for x,y in self.g_seq.items() if x not in seqnames})
            seqnames += self.l_csv["sequence_name"].tolist() # add all local to avoid being pruned 
        newtrees = []
        for i,t in enumerate(self.trees):
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
                logger.debug("and they are:\n%s", "\n".join(remove_leaf))
                self.trees.append(t.extract_tree_without(remove_leaf))
            else:
                self.trees.append(t)
        # copy of user-defined trees that don't get pruned
        if len(self.trees) > 1:
            self.usrtrees = []
            for t in self.trees[1:]:
                self.usrtrees.append( treeswift.read_tree_newick (str(t)) ) ## copy 

    def remove_duplicates (self, blocks = 4, leaf_size = 500, radius=0.00001):
        logger.info("Removing duplicates (identical sequences)")
        clusters = ml.list_duplicates (self.g_snp, blocks, leaf_size, radius)

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
        self.g_seq, self.g_snp, self.trees = self.remove_seq_tree_based_on_metadata()
    
    def reduce_redundancy (self, clade_rule = None):
        if clade_rule is None: # for each lineage level with at least x[0] samples, keep at most x[1]
            clade_rule = [ # rules can be repeated, for different thresholds; some samples fall into several 
                    ["lineage",    1, 20], # only those with >1 samples; then take up to 10
                    ["lineage",   20, 500], # only those with >20 samples; then take up to 500 
                    ["acc_lineage",1, 50],
                    ["del_lineage",1, 50],
                    ["adm1",       5, 500],
                    ["uk_lineage", 2, 20], 
                    ["uk_lineage",20, 200], 
                    ["phylotype",  2, 20],
                    ["acc_lineage", 2, 10],
                    ["del_lineage", 2, 10]
                    ] 
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
        for column in ["lineage", "uk_lineage", "acc_lineage", "del_lineage"]:
            if column in df.columns:
                df1 = df[ df[column] == "" ].head(100)  # undefined (missing) lineages 
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

    def find_neighbours (self, blocks = 2000, leaf_size = 500, dist_blocks = 4, nn = 25):
        logger.info(f"Finding neighbours to local sequences, using a distance of {dist_blocks} to {blocks} segments")
        neighbours1 = ml.list_r_neighbours (self.g_snp, self.l_snp, blocks, leaf_size, dist_blocks)
       
        blocks = 2 * blocks; leaf_size = leaf_size/2
        logger.info("Found %s neighbours; now will find their %s closest neighbours on %s segments", 
                len(neighbours1), str(nn), str(blocks))
        aln_d = {x:self.g_snp[x] for x in neighbours1}
        neighbours = ml.list_n_neighbours (self.g_snp, aln_d, blocks, leaf_size, nn)
        neighbours = list(set(neighbours + neighbours1))
        logger.info("Found %s neighbours", len(neighbours))
        return neighbours

def save_global_from_seqnames (bb, seqnames, prefix):
    seqs, snps, trees = bb.remove_seq_tree_based_on_metadata (seqnames)
    desc = save_sequences (seqs, prefix)
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
    aln.update({x:y for x,y in bb.usrseqs.items() if x in leafnames})
    trees = bb.usrtrees

    if len(aln) < len(leafnames):
        remove_leaf = [l.get_label() for t in bb.usrtrees for l in t.traverse_leaves() if l.get_label() not in aln.keys()]
        logger.info("In total %s user-defined leaves do not have a sequence, and will be pruned",str(len(remove_leaf)))
        logger.debug("And they are\n%s\n", "\t".join(remove_leaf))
        trees = [t.extract_tree_without(remove_leaf) for t in bb.usrtrees]

    ## add sequences from current analysis (all local plus global neighbours)
    aln.update({x:y for x,y in bb.g_seq.items() if x in seqnames})
    aln.update({x:y for x,y in bb.l_seq.items()})
    desc = save_sequences (aln, prefix)
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
    desc = save_sequences (aln, prefix)
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
    return os.path.basename(fname)

def read_peroba_database (f_prefix): 
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
    logger.info(f"Reading database sequences from \'{fname}\'")
    sequences = common.read_fasta (fname, check_name = False)

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequences, tree]

def main_generate_backbone_dataset (database, csv, sequences, trees, replace, prefix):
    # initialisation
    bb = PerobaBackbone (database)
    bb.add_local_data_and_sequences (csv, sequences, replace)
    bb.finalise_and_split_data_sequences()
    bb.add_trees (trees)
    bb.remove_duplicates()
    ## these two reduce data set, but ideally should return a copy (TODO)
    bb.remove_low_quality()
    bb.reduce_redundancy() ## add step to remove too distant 
    # more methods come here
    neighbours = bb.find_neighbours()
    # for each method we can use chosen neighbours to save global only or both
    description  = save_global_from_seqnames (bb, neighbours, prefix + "coguk")
    description += save_user_trees (bb, neighbours, prefix + "user", add_nj_tree = True)
    description += save_all_sequences (bb, neighbours, prefix + "norw-coguk", add_nj_tree = True)
    # finally, save all NORW sequences
    desc = save_sequences (bb.l_seq, prefix + "norw")
    description += "{}\n alignment with all local sequences only\n".format(colour_string(desc))

    print ("Finished. The output files are described below, where 'global' means COGUK and GISAID data which were ")
    print ("{} generated in NORW. Those, together with the extra sequences are being called'local'.".format(colour_string("not", "red")))
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
    parser.add_argument('-t', '--trees', metavar='', help="file with (user-defined) trees in newick format to help produce backbone")
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
            logger.warning (f"Could not find local CSV file {args.csv}; Will proceed without it")
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
            logger.warning (f"Could not find sequence file {args.sequences}; Will proceed without it")
        else:
            logger.info("Reading fasta file with sequences from NORW")
            sequences = common.read_fasta (fname, check_name = False)

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

    main_generate_backbone_dataset (database, csv, sequences, trees, args.replace, prefix)


if __name__ == '__main__':
    main()
