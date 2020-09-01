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

class PerobaSubsample:
    csv = None
    csv_new = None # new columns only 
    seq = None # dictionary
    tree = None 
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

    def __init__ (self, peroba_db): # peroba_db = [metadata, aligned seqs, tree]  
        self.csv = peroba_db[0]  # formatted metadata
        self.csv0 = self.csv.copy()
        self.seq = {x.id:x for x in peroba_db[1]} # dictionary
        self.tree = peroba_db[2] 
        cols = [x for x in self.cols if x in peroba_db[0].columns]
        self.csv = self.csv[cols] # remove other columns
        logger.info("Imported %s rows from database", str(self.csv.shape[0]))

    def trim_sequences (self, trim = True):  
        if trim is True: trim = [265, 29675]
        if trim[0] < 1: trim[0] = 1
        if trim[1] > 29902: trim[1] = 29902 # genome size is 29903
        logger.info("Trimming genomes from site %s to site %s",str(trim[0]), str(trim[1]))
        for x in self.seq.values():
            x.seq = x.seq[trim[0]:trim[1]]
    
    def sort_categories (self):
        df = self.csv
        if "adm2" in df.columns:
            df["adm2"] = df["adm2"].replace(["Unknown Source","Unknown", np.nan],"")
            df["adm2"] = df["adm2"].replace({"Greater London":"Greater_London", "Hertfordshire":"Herefordshire"})
        
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col])
        
        for col in fill_missing_cols:
            if col in df.columns:
                df[col].fillna ("", inplace = True) 

        pref = [x for x in prefsort if x[0] in df.columns] # sort only existing columns

        ## local -- assume sequences from NORW or UK at least
        logger.info("Ordering metadata after categorisation -- UK-centric order") 
        sort_cols = [x[0] for x in pref] # follow this order (thus can't use dict...)
        sort_ascend = [x[1] for x in pref] # same order as sort_cols
        for col,favs in favourite_cats.items(): # favourite_cats happens to have local columns as keys
            if col in df.columns:
                values = [x for x in df[col].unique() if x not in favs + [""]] 
                values = favs + values  + [""] # make sure favs are first in list, and empty are last
                df[col] = pd.Categorical(df[col], categories = values, ordered=True)
        df = df.sort_values(by=sort_cols, ascending=sort_ascend)
        df['peroba_local_order'] = np.arange(len(df)) # new column
        self.csv = self.csv.combine_first(df[['peroba_local_order']])

        ## global -- assme seqs from outside UK, do not focus on uk_lineage or NORW regions
        logger.info("Ordering metadata after categorisation -- global order") 
        df = self.csv 
        pref = [x for x in pref if x not in favourite_cats.values()] # exclude uk centric preferences
        sort_cols = [x[0] for x in pref] # follow this order (thus can't use dict...)
        sort_ascend = [x[1] for x in pref] # same order as sort_cols
        df = df.sort_values(by=sort_cols, ascending=sort_ascend)
        df['peroba_global_order'] = np.arange(len(df)) # new column
        self.csv = self.csv.combine_first(df[['peroba_global_order']])

    def finalise_data_sequences (self, trim = True, strict = False):
        self.trim_sequences(trim=trim)
        self.sort_categories()
        #logger.info("Finding SNPs to speed up calculations")
        #snp_aln = snpsites_from_alignment ([x for x in self.seq.values()], strict=strict)
        #logger.info("In total %s SNPs were found (i,e, alignment size)", len(snp_aln[0].seq))
        #self.snp = {x.id:x for x in snp_aln}
        
        #l_snp = [x for x in snp_aln if x.id in self.l_csv["sequence_name"]]
        #idx = sorted_uncertainty_from_alignment (l_snp, max_freq_n = 0.01)
        #logger.info("Compact representation uses %s SNPs (columns from local seqs with fewer Ns)", len(idx))
        #snp_aln = alignment_from_column_index (snp_aln, idx)

    def remove_seq_tree_based_on_metadata (self, seqnames = None): 
        if seqnames is None: 
            logger.info("removing sequences after metadata update")
            seqnames = self.csv["sequence_name"].tolist()
        newseqs = {x:y for x,y in self.seq.items() if x in seqnames}
        logger.info("Finding tree leaves to be pruned (this step is slow)")
        remove_leaf = []
        for lab in self.tree.labels(internal=False): 
            if lab not in seqnames:
                remove_leaf.append(lab)
        if len(remove_leaf) > 0:
            logger.warning(" %s leaves from tree were pruned (absent from metadata or list)", str(len(remove_leaf)))
            newtree = self.tree.extract_tree_without(remove_leaf)
        else:
            newtree = self.tree
        return  newseqs, newtree

    def remove_low_quality (self, f_acgt = 0.8, f_n = 0.1):
        logger.info(f"Removing sequences with proportion of ACGT less than {f_acgt} or proportion of N higher than {f_n}")
        self.csv = self.csv.loc[ (self.csv["peroba_freq_acgt"] > f_acgt) & (self.csv["peroba_freq_n"] < f_n) ]
        if "host" in self.csv.columns:
            self.csv = self.csv[ self.csv["host"].str.contains("Human",case=False,na=True) ]
        
        logger.info("Database has now %s rows (after removing low-quality and non-human samples)", str(self.csv.shape[0]))
        self.seq, self.tree = self.remove_seq_tree_based_on_metadata()

    def remove_duplicates (self, blocks = 2, leaf_size = 500, radius=0.00001):
        logger.info("Removing duplicates (identical sequences) from data")
        clusters = ml.list_duplicates (self.seq, blocks, leaf_size, radius)

        df = pd.DataFrame([[y,i] for i,x in enumerate(clusters) for y in x], columns=["sequence_name","peroba_tmp"])
        df = df.sort_values(by=["peroba_tmp"], ascending=True) # same sequence belongs to several clusters
        df = df.groupby("sequence_name").aggregate("first") # now each sequence belongs only to lowest cluster number
        logger.info("Number of unique sequences: %s", str(len(df["peroba_tmp"].unique())) )

        self.csv.reset_index (drop=False, inplace=True) ## merge will destroy index...
        df = self.csv.merge(df, on="sequence_name")
        df = df.sort_values(by="peroba_global_order", ascending=True) # break ties by global order (i.e. not giving preference to UK)

        ## create list of accepted seqnames instead of deleting rows; replace peroba_tmp by peroba_group
        #for c in [x for x in replace_duplicate_cols if x in df.cols]: # most common value amongst identical sequences
        #    df[c] = df.groupby(["peroba_tmp"])[c].transform(pd.Series.mode) # transform() keeps index

        df = df.groupby("peroba_tmp").aggregate("first") # only one sequence from each cluster, following self.order_col preference
        df.set_index ("peroba_seq_uid", drop = True, inplace = True)
        self.csv = df
        self.seq, self.tree = self.remove_seq_tree_based_on_metadata()
    
    def reduce_redundancy (self, extended_mode, new_col_name):
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
        if extended_mode == 1: # more permissive since many non-UK sequences don't have this info
            clade_rule = [[x[0], int(x[1]/2), 4 * x[2]] for x in clade_rule] # zero is fine 
            missing_rule = [[x[0], 16 * x[1]] for x in missing_rule]

        if extended_mode == 2: # more permissive since many non-UK sequences don't have this info
            clade_rule = [[x[0], int(x[1]/2), 8 * x[2]] for x in clade_rule] # zero is fine 
            missing_rule = [[x[0], 32 * x[1]] for x in missing_rule]

        df = self.csv
        logger.info(f"Subsampling redundant samples (same lineage etc.) for column {new_col_name}") 
        if extended_mode == 0:
            df = df.sort_values(by="peroba_local_order", ascending=True)
        else:
            df = df.sort_values(by="peroba_global_order", ascending=True)
        
        dfcat = None
        for column, rule1, rule2 in clade_rule:
            if column in df.columns:
                # groupby().filter() creates another DF (excluding rows with less than rule1) therefore a second groupby is needed: groupby().head()
                df1 = df.groupby(column).filter(lambda x: len(x.index) > rule1).groupby(column).head(rule2)
                if dfcat is None: dfcat = df1
                else: dfcat = pd.concat([dfcat, df1])
        for column, n_keep in missing_rule: 
            if column in df.columns:
                df1 = df[ df[column] == "" ].head(n_keep)  # undefined (missing) lineages 
                if dfcat is None: dfcat = df1
                else: dfcat = pd.concat([dfcat, df1])

        dfcat = dfcat.groupby(dfcat.index).aggregate("first") # many rows will be duplicated before this groupby()
        logger.info("After subsampling, global metadata has %s samples", dfcat.shape[0])
        dfcat[new_col_name] = 1
        dfcat = dfcat[[new_col_name]]
        if (self.csv_new is None):
            self.csv_new = dfcat
        else:
            self.csv_new = self.csv_new.combine_first(dfcat)
        self.csv_new[new_col_name] = self.csv_new[new_col_name].fillna(0)
#        self.seq, self.snp, self.tree = self.remove_seq_tree_based_on_metadata()

    def find_pda_samples (self, fraction_remain = 0.5):
        n_remain = int (fraction_remain * self.csv.shape[0]) + 1
        if n_remain < 16:
            return
        logger.info(f"Finding the {n_remain} most distant leaves in the tree")
        leafnames = pda_names_from_tree (self.tree, n_remain = n_remain)
        df = pd.DataFrame({"peroba_seq_uid":leafnames, "peroba_pda":1})
        df.set_index("peroba_seq_uid", inplace=True)
        if (self.csv_new is None):
            self.csv_new = df
        else:
            self.csv_new = self.csv_new.combine_first(df)
        self.csv_new["peroba_pda"] = self.csv_new["peroba_pda"].fillna(0)


    def save_subsample (self, f_prefix):
        self.csv_new["peroba_subsample"] = 1 # this all-one column will be used once merged with big metadata 
        fname = f_prefix + common.suffix["subsample"]
        self.csv_new.to_csv(fname)


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
    logger.info(f"Reading database alignment from \'{fname}\'")
    sequences = common.read_fasta (fname, check_name = False)

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequences, tree]

def main_generate_subsample_dataset (f_prefix):
    # initialisation
    bb = PerobaSubsample (read_peroba_database (f_prefix))
    bb.finalise_data_sequences()
    bb.remove_low_quality()
    bb.remove_duplicates()
    bb.reduce_redundancy(0, "peroba_level_0") 
    bb.reduce_redundancy(1, "peroba_level_1") 
    bb.reduce_redundancy(2, "peroba_level_2") 
    bb.save_subsample (f_prefix)

    
class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="""
    peroba_subsample is the script that generates a smaller global backbone data set (COGUK+GISAID) from the full one.
    It depends on the prefix for a perobaDB set of files (from `peroba_database`), like "perobaDB.0519".
    """, 
    usage='''peroba_subsample <perobaDB> [options]''')

    parser.add_argument('perobaDB')
    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
        help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-i', '--input', action="store", help="Directory where perobaDB files are. Default: working directory")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    if args.input: input_d = os.path.join(current_working_dir, args.input)
    else: input_d = current_working_dir

    logger.info("Reading metadata, sequences, and tree from peroba_database")
    main_generate_subsample_dataset (os.path.join(input_d, args.perobaDB))

if __name__ == '__main__':
    main()
