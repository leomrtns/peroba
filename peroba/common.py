import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm, colors, patches

## nabil's needed columns: peroba_uk_lineage  peroba_lineage peroba_phylotype  peroba_special_lineage central_sample_id
## TODO: "UNKNOWN SOURCE" and "UNKNOWN" are the same adm2 (in cog) (fixed in backbone)
## TODO: columns w/ values in csv but not metadata (e.g. adm_private or this week's) may not ASreconstructed
## (not merged?)


import logging, ete3
import numpy as np, pandas as pd
from Bio import Seq, SeqIO
import datetime, sys, lzma, gzip, bz2, re, glob, collections, subprocess, os, itertools, pathlib

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

prefix = {
        "database": "perobaDB." ## timestamp comes here (watch for double dots)
        }
suffix = {
        "metadata": ".metadata.csv.gz",
        "raw": ".raw.csv.gz",
        "tree": ".tree.nhx",
        "sequences": ".sequences.fasta.xz",
        "alignment": ".sequences.aln.xz"
        }
        
#asr_cols = ["adm2", "uk_lineage", "lineage", "phylotype", "submission_org_code", "date_sequenced", "source_age", "source_sex", "collecting_org", "ICU_admission"]

# For nabil, these must be on our output csv: peroba_uk_lineage    peroba_lineage    peroba_phylotype    peroba_special_lineage  central_sample_id
# therefore these go to master sheet and come back without "peroba_" on following week
asr_cols = ["uk_lineage", "lineage", "phylotype"]  
# therefore we exclude the following columns from local (a.k.a. master table), since these must be from global metadata if present
remove_from_master_cols = ["uk_lineage", "lineage", "phylotype", "special_lineage", "acc_lineage", "del_lineage"]


def read_ete_treefile (treefile, multi = None):
    if multi is None: 
        multi = False
    if multi:
        logger.info("Reading file %s with set of trees and checking for duplicate names", treefile)
        treestring = [x.rstrip().replace("\'","").replace("\"","").replace("[&R]","") for x in open(fname)]
    else:
        logger.info("Reading tree file %s (only first tree will be used) and checking for duplicate names", treefile)
        treestring = [open(treefile).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")]
    trees = []
    for i,trs in enumerate (treestring): # check duplicated names (ete3 accepts without checking, but complains later)
        tre = ete3.Tree(trs)
        tree_length = len([leaf.name for leaf in tre.iter_leaves()])
        tree_leaves = {str(leaf.name):leaf for leaf in tre.iter_leaves()} # dup leaves will simply overwrite node information
        if (tree_length > len(tree_leaves)):
            tre.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
            logger.warning(f"Found duplicated leaf names in input treefile id={i}, will keep one at random")
        logger.info("%s leaves in treefile id=%s", len(tree_leaves), str(i))
        trees.append(tre)
    if multi: return trees
    else: return trees[0]

def read_fasta (filename, fragment_size = 0, check_name = False):
    unaligned = []
    if   "bz2" in filename[-5:]: this_open = bz2.open
    elif "gz"  in filename[-5:]: this_open = gzip.open
    elif "xz"  in filename[-5:]: this_open = lzma.open
    else:  this_open = open
    with this_open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record.seq  = record.seq.upper()
            # record.name = name_without_spaces.hash() ## IDEA
             # # seq names are a mess, it's better to map using metadata as key (however we have keys as "D02" etc
            # COG-UK/SHEF-D34A0/SHEF:20200423_1347_X1_FAN43269_d941f767|SHEF-D34A0|... 
            # COG-UK/OXON-AF08B/OXON
            # EPI_ISL_422016|hCoV-19/Wales/PHWC-26796/2020|Wales|WALES|2020-03-28|PHWC|...
            # hCoV-19/England/201140442/2020|PENDING|England|GREATER_LONDON|2020-03-10|PHEC|2
            # hCoV-19/England/CAMB-74D3D/2020||England|CAMBRIDGESHIRE|2020-03-19|  also hCoV Hcov hcoV etc.
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                seqname = record.description.split("|")[0]
                record.id = seqname.replace("hCoV-19/","")
            if len(record.seq) > fragment_size:
                unaligned.append(record)
    logger.info("Read %s sequences from file %s", str(len(unaligned)), filename)
    return unaligned

def df_read_genome_metadata (filename, primary_key = "sequence_name", rename_dict = "default", sep=',', compression="infer", 
        index_name = "peroba_seq_uid", remove_edin = True, exclude_na_rows = None, exclude_columns = "default"):
    ''' Read CSV data using "sequence_name" as primary key 
    TODO: there are sequences with same name (e.g. distinct assemblies), currently only one is kept. 
    I could create a column combination (multi-index?) to keep all
    '''
    if rename_dict == "default":
        rename_dict = {
            'sample_date':'collection_date',  # from cog_gisaid
            'strain':'sequence_name',         # others below from GISAID
            'date':'collection_date',
            'sex':'source_sex',
            'age':'source_age'}
    if exclude_na_rows == "default": ## better NOT to set this up here (only at finalise_data())
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    df1 = pd.read_csv (str(filename), compression=compression, sep=sep, dtype='unicode') # sep='\t' for gisaid

    # fix column names  
    if rename_dict:
        for old, new in rename_dict.items():
            if old in list(df1.columns): ## I expected pandas to do this check...
                df1.rename(columns={old:new}, inplace=True)
    if remove_edin: # EDINBURGH custom labels have an "edin" prefix
        no_edin = {x:x.replace("edin_", "") for x in list(df1.columns) if x.startswith("edin_")}

    # reduce table by removing whole columns and rows with missing information
    if exclude_na_rows: # list of important columns, such that rows missing it should be removed
        important_cols = [x for x in exclude_na_rows if x in list(df1.columns)]
        if important_cols:
            df1.dropna (subset = important_cols, inplace = True);

    #EX: exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]
    if exclude_columns: # list of irrelevant columns
        irrelevant_cols = [x for x in exclude_columns if x in list(df1.columns)]
        if irrelevant_cols:
            df1.drop (labels = irrelevant_cols, axis=1, inplace = True) # no sample with this information

     # This must be relaxed if we want global statistics (no purging of incomplete rows)
    if index_name not in list(df1.columns): # regular CSV file from external source
        df1.set_index (str(primary_key), drop = False, inplace = True) # dont drop the column to be used as index
        # now we have a column and an index with same name
        df1.dropna (subset=[str(primary_key)], inplace = True); 
        df1.rename_axis(str(index_name), inplace = True) # equiv. to df.index.name="peroba_seq_uid"
    else:
        df1.set_index (str(index_name), drop = True, inplace = True) # drop the column to avoid having both with same name
    
    if "collection_date" in df1.columns: ## keep most recent first (in case one has wrong date)
        df1["collection_date"] = pd.to_datetime(df1["collection_date"], infer_datetime_format=True, yearfirst=True, errors='coerce')
        df1.sort_values(by=["collection_date"], inplace = True, ascending=False)
    
    if "adm2" in df1.columns:
        df1['adm2'] = df1['adm2'].str.title()
    df1 = df1.groupby(df1.index).aggregate("first"); # duplicated indices are allowed in pandas
    return df1

def df_merge_metadata_by_index (df1, df2): # replaces NA whenever possible, using INDEX, and adds new columns from df2
    df2 = df2.combine_first(df1) # adds df1 columns; df1 informs missing values 
    df1 = df1.combine_first(df2) # improved df2 thus fills missing values from df1
    if "collection_date" in df1.columns: ## GISAID sometimes have old (i.e. missing) datetime from coguk
        df1.sort_values(by=["collection_date"], inplace = True, ascending=False) 
    # if before combine_first() there were dups they're kept by pandas
    df1 = df1.groupby(df1.index).aggregate("first"); 
    #df1.sort_index(inplace=True)
    return df1

def df_finalise_metadata (df, exclude_na_rows = None, exclude_columns = "default", remove_duplicate_columns = True):
    ''' if exclude_na_rows is the string "default" then we use sensible choices
    '''
    if exclude_na_rows == "default":
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    # df['days_since_Dec19'] = df['collection_date'].map(lambda a: get_days_since_2019(a, impute = True)) ## not here
    # df['days_since_Dec19'] = df['collection_date'].map(lambda a: datetime.datetime(a) - datetime.datetime(2019,12, 1).days) 
    df["collection_date"] = pd.to_datetime(df["collection_date"], infer_datetime_format=True, yearfirst=True, errors='coerce')
    cols = [x for x in ['lineage_support', 'collection_date'] if x in df.columns]
    df = df.sort_values(by=cols, ascending=False)
    
    ## CURRENTLY NOT WORKING since sequence_name was used as primary key... 
    # COGUK phylo analysis removes low-quality seqs, which are deleted from metadata as well. (COGUK phylo assigns decent senames)
    #   I give another chance for these deleted sequences, if their names correspond to sample_central_id 

    #df.reset_index (drop=False, inplace=True) ## drop=False will make peroba_seq_uid become a column
    #df['sequence_name'].fillna("sample_central_id", inplace=True)
    #df['peroba_seq_uid'].fillna("sample_central_id", inplace=True)
    #df.set_index ("peroba_seq_uid", drop = True, inplace = True) # drop deletes the extra peroba_seq_uid column

    # default values for missing rows (in case we don't want to remove those rows)
    #if not exclude_na_rows or "uk_lineage" not in exclude_na_rows:
    #    df["uk_lineage"] = df["uk_lineage"].replace(np.nan, "x", regex=True) 
    if not exclude_na_rows or "adm2" not in exclude_na_rows:
        df['adm2'].fillna(df.country, inplace=True)
    df['adm2'] = df['adm2'].str.title()

    if exclude_na_rows: # list of important columns, such that rows missing it should be removed
        important_cols = [x for x in exclude_na_rows if x in list(df.columns)]
        if important_cols:
            df.dropna (subset = important_cols, inplace = True);
    if exclude_columns: # list of irrelevant columns
        irrelevant_cols = [x for x in exclude_columns if x in list(df.columns)]
        if irrelevant_cols:
            df.drop (labels = irrelevant_cols, axis=1, inplace = True) # no sample with this information

    if remove_duplicate_columns: # who knows which column it will choose (i.e. column names)
        df = df.T.drop_duplicates().T # transpose, remove duplicate rows, and transpose again
    return df

def add_sequence_counts_to_metadata (metadata, sequences, from_scratch = None):
    if from_scratch is None: from_scratch = True
    """ counts the proportions of indel/uncertain bases (i.e `N` and `-`) as well as number of ACGT.
    Notice that they do not sum up to one since we have partial uncertainty (`W`,`S`, etc.) that can be used
    This should be done on UNALIGNED sequences
    """
    def calc_freq_N (index):
        if index not in sequences: return 1.
        if sequences[index] is None: return 1. ## missing sequences
        genome = str(sequences[index].seq); l = len(genome)
        if (l):
            number_Ns = sum([genome.upper().count(nuc) for nuc in ["N", "-"]])
            return number_Ns / l
        else: return 1.
    def calc_freq_ACGT (index):
        if index not in sequences: return 0.
        if sequences[index] is None: return 0. ## missing sequences
        genome = str(sequences[index].seq); l = len(genome)
        if (l):
            number_ACGTs = sum([genome.upper().count(nuc) for nuc in ["A", "C", "G", "T"]])
            return number_ACGTs / len(genome)
        else: return 0.

    if from_scratch or "peroba_freq_n" not in metadata.columns: # map() sends the index to lambda function
        metadata["peroba_freq_n"] = metadata.index.map(lambda x: calc_freq_N (x))  
    else: # only update null values 
        nilvalues = metadata[metadata["peroba_freq_n"].isnull()].index.map(lambda x: calc_freq_N (x))  
        if len(nilvalues) > 0:
            metadata.loc[metadata["peroba_freq_n"].isnull(), "peroba_freq_n"] = nilvalues 
    
    if from_scratch or "peroba_freq_acgt" not in metadata.columns: # map() sends the index to lambda function
        metadata["peroba_freq_acgt"] = metadata.index.map(lambda x: calc_freq_ACGT (x))  
    else: # only update null values 
        nilvalues = metadata[metadata["peroba_freq_acgt"].isnull()].index.map(lambda x: calc_freq_ACGT (x))  
        if len(nilvalues) > 0:
            metadata.loc[metadata["peroba_freq_acgt"].isnull(), "peroba_freq_acgt"] = nilvalues 

    return metadata, sequences

def transparent_cmap (color=None, cmap=None, final_alpha=None):
  # http://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
  if color is None:
    color = "blue"
  if (final_alpha is None) or (final_alpha < 0.01):
    final_alpha = 1.
  if cmap:
    mycmap = plt.get_cmap(cmap)
  else:
    from matplotlib.colors import colorConverter
    mycmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap', [colorConverter.to_rgba(color),colorConverter.to_rgba(color)],256)

  mycmap._init() # create the _lut array, with rgba values
  mycmap._lut[:,-1] = np.linspace(0, final_alpha, mycmap.N+3) # here it is progressive alpha array
  return mycmap

def continuous_cmap (cmap = None, list_size = None):
    """ if list_size is not null, returns a list of this size with colours from cmap; otherwise return the cmap
    """
    if cmap is None: cmap = "plasma"
    if list_size is None: list_size = 0
    ncls = {'Pastel1':9, 'Pastel2':8, 'Paired':12, 'Accent':8,'Dark2':8, 'Set1':9, 
            'Set2':8, 'Set3':12,'tab10':10, 'tab20':20, 'tab20b':20, 'tab20c':20}
    custom_cmap = cm.get_cmap(cmap)
    if cmap in ncls.keys():
        custom_cmap = colors.LinearSegmentedColormap.from_list("custom", [(x,custom_cmap(x)) for x in np.linspace(0, 1, ncls[cmap])])
    if list_size == 0:
        return custom_cmap
    return custom_cmap(np.linspace(0, 1, list_size))

