import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm, colors, patches

## nabil's needed columns: peroba_uk_lineage  peroba_lineage peroba_phylotype  peroba_special_lineage central_sample_id
## TODO: "UNKNOWN SOURCE" and "UNKNOWN" are the same adm2 (in cog) (fixed in backbone)
## TODO: columns w/ values in csv but not metadata (e.g. adm_private or this week's) may not ASreconstructed
## (not merged?)

import logging, ete3
import numpy as np, pandas as pd
from Bio import Seq, SeqIO
import random, datetime, sys, lzma, gzip, bz2, re, glob, collections, subprocess, itertools, pathlib, base64
import pandas_profiling # ProfileReport
from peroba.utils import * 

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

qi_logo_file = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/report/QIlogo.png") 
 
# For nabil, these must be on our output csv: peroba_uk_lineage    peroba_lineage    peroba_phylotype    peroba_special_lineage  central_sample_id
# therefore these go to master sheet and come back without "peroba_" on following week
asr_cols = ["uk_lineage", "lineage", "phylotype", "primary_uk_lineage"]  
# values imputed by peroba but that make no sense unless the sample is from the UK (with the "peroba_") 
uk_specific_cols = ["uk_lineage", "phylotype", "primary_uk_lineage"]  
# therefore we exclude the following columns from local (a.k.a. master table), since these must be from global metadata if present
remove_from_master_cols = ["uk_lineage", "lineage", "phylotype", "special_lineage", "acc_lineage", "del_lineage", "primary_uk_lineage"] 
# primary_uk_lineage is new on 2020.06.08

dtype_numeric_cols = [
        'gaps', 'length', 'edin_epi_week', 'epi_week', 'layout_insert_length', 'layout_read_length', 'missing', 
        'source_age', 'virus', 'lineage_support',  # upstream  (coguk+gisaid)
        "peroba_freq_acgt", "peroba_freq_n", "PCR Ct value", "Coverage (X)",  # local (peroba, NORW)
        "ct_1_ct_value", "ct_2_ct_value", "No. Reads", "Mapped Reads", "No. Bases (Mb)",
        "Average read length", "Missing bases (N)", "Consensus SNPs"]

dtype_datetime_cols = ["date_submitted", "collection_date", "received_date", "sequencing_submission_date", "start_time",
      "date_sequenced"] # from NORW metadata
# "Sequencing date" is coded (20200516) but not equal to "date_sequenced" 

## true false Basic QC , High Quality QC,
    
def metadata_to_html (df0, filename, description):
    df = set_dtypes_metadata (df0)
    df.dropna  (axis=1, how='all', inplace=True) # delete empty columns
    covrge = "Coverage (X)" ## must remove ending X (e.g. 2000X)
    if covrge in df.columns:
        df[covrge] =  pd.to_numeric(df[covrge].str[:-1])

    profile = pandas_profiling.ProfileReport (df, title = description, explorative=True)
    profile.set_variable("html.style.primary_color", "#00a990") # #bed12b avocado #00a990 bluegreen #333f48 dark green
    logo = base64.b64encode(open(qi_logo_file, 'rb').read()) 
    logo = "data:image/png;base64," + logo.decode('utf-8')
    profile.set_variable("html.style.logo", logo) 
    profile.set_variable("html.style.theme", "flatly") 
    profile.set_variable("correlations.cramers.calculate", False) # get a warning otherwise  
    profile.to_file(filename)

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

def align_sequences_in_blocks (sequences, reference_file, seqs_per_block = 2000): 
    if seqs_per_block < 100:   seqs_per_block = 100
    if seqs_per_block > 10000: seqs_per_block = 10000 # mafft chokes on large matrices (but 10k is fine btw)
    nseqs = len(sequences)
    aligned = []
    for i in range (0, nseqs, seqs_per_block):
        last = i + seqs_per_block
        if last > nseqs: last = nseqs
        aln = mafft_align_seqs (sequences[i:last], reference_file = reference_file)
        aligned += aln
        logger.info (f"First {last} sequences aligned")
    return aligned

def set_dtypes_metadata (df0):
    df = df0.copy()
    columns = [x for x in dtype_numeric_cols if x in df.columns]
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    columns = [x for x in dtype_datetime_cols if x in df.columns]
    for col in columns:
        df[col] = pd.to_datetime(df[col], infer_datetime_format=True, yearfirst=True, errors='coerce')
    return df

def replace_values_metadata (df0):
    df = df0.copy()
    if "source_sex" in df.columns:
        df["source_sex"] = df["source_sex"].replace(["Woman", "Female", "FEmale"],"F")
        df["source_sex"] = df["source_sex"].replace(["Male"],"M")
        df["source_sex"] = df["source_sex"].replace(["Unknown", "unknwon", "U"],"?")
    if "adm2" in df.columns:
        df['adm2'] = df['adm2'].str.title()
        df["adm2"] = df["adm2"].replace(["Unknown Source","Unknown"],"")
        df["adm2"] = df["adm2"].replace({"Greater_London":"Greater London"}) # "Hertfordshire" != "Herefordshire"
        #df['adm2'].fillna(df.country, inplace=True)
    if "is_icu_patient" in df.columns:
        df["is_icu_patient"] = df["is_icu_patient"].str.replace("Unknown","?")

    #df["uk_lineage"] = df["uk_lineage"].replace(np.nan, "x", regex=True) 
    return df
 
def rename_columns_metadata (df0):
    df = df0.copy()
    rename_dict = {
        'sample_date':'collection_date',  # from cog_gisaid
        'strain':'sequence_name',         # others below from GISAID
        'date':'collection_date',
        'sex':'source_sex',
        'age':'source_age'}
    for old, new in rename_dict.items():
        if old in df.columns and new not in df.columns: # cogul has a new column "sample_date
            df.rename(columns={old:new}, inplace=True)
    no_edin = {x:x.replace("edin_", "") for x in list(df.columns) if x.startswith("edin_")}
    df.rename(columns=no_edin, inplace=True)
    return df

def exclude_irrelevant_columns_metadata (df0):  # no samples with this information yet (all missing values) 
    df = df0.copy()
    exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "metadata" ] #  "travel_history" has a few
    exclude_columns = [x for x in exclude_columns if x in list(df.columns)]
    if exclude_columns:
        df.drop (labels = exclude_columns, axis=1, inplace = True)
    return df

def remove_rows_with_missing_metadata (df0):  # currently not used, since too restrictive 
    df = df0.copy()
    exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2']
    # reduce table by removing whole columns and rows with missing information
    exclude_na_rows = [x for x in exclude_na_rows if x in list(df.columns)]
    if exclude_na_rows:
        df.dropna (subset = exclude_na_rows, inplace = True);
    return df

def df_read_genome_metadata (filename, primary_key = "sequence_name", sep=',', index_name = "peroba_seq_uid", 
        exclude_na_rows = False, exclude_columns = True):
    ''' Read CSV data using "sequence_name" as primary key 
    TODO: there are sequences with same name (e.g. distinct assemblies), currently only one is kept. 
    I could create a column combination (multi-index?) to keep all
    '''

    df1 = pd.read_csv (str(filename), compression="infer", sep=sep, dtype='unicode') # sep='\t' for gisaid

    df1 = rename_columns_metadata (df1) # fix column names before working with them, below
    df1 = set_dtypes_metadata (df1)
    df1 = replace_values_metadata (df1)
    if exclude_columns:
        df1 = exclude_irrelevant_columns_metadata (df1)
    if exclude_na_rows: # list of important columns, such that rows missing it should be removed
        df1 = remove_rows_with_missing_data (df1)

     # This must be relaxed if we want global statistics (no purging of incomplete rows)
    if index_name not in list(df1.columns): # regular CSV file from external source
        df1.set_index (str(primary_key), drop = False, inplace = True) # dont drop the column to be used as index
        # now we have a column and an index with same name
        df1.dropna (subset=[str(primary_key)], inplace = True); 
        df1.rename_axis(str(index_name), inplace = True) # equiv. to df.index.name="peroba_seq_uid"
    else:
        df1.set_index (str(index_name), drop = True, inplace = True) # drop the column to avoid having both with same name
    
    if "collection_date" in df1.columns: ## keep most recent first (in case one has wrong date)
        df1.sort_values(by=["collection_date"], inplace = True, ascending=False)
    
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

def df_finalise_metadata (df, exclude_na_rows = False, exclude_columns = False, remove_duplicate_columns = True):
    ''' if exclude_na_rows is the string "default" then we use sensible choices
    '''
    if exclude_na_rows == "default":
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    cols = [x for x in ['lineage_support', 'collection_date'] if x in df.columns]
    df = df.sort_values(by=cols, ascending=False)
    
    ## CURRENTLY NOT WORKING since sequence_name was used as primary key... 
    # COGUK phylo analysis removes low-quality seqs, which are deleted from metadata as well. (COGUK phylo assigns decent senames)
    #   I give another chance for these deleted sequences, if their names correspond to sample_central_id 
    #df.reset_index (drop=False, inplace=True) ## drop=False will make peroba_seq_uid become a column
    #df['sequence_name'].fillna("sample_central_id", inplace=True)
    #df['peroba_seq_uid'].fillna("sample_central_id", inplace=True)
    #df.set_index ("peroba_seq_uid", drop = True, inplace = True) # drop deletes the extra peroba_seq_uid column

    if exclude_columns:
        df = exclude_irrelevant_columns_metadata (df)
    if exclude_na_rows: 
        df1 = remove_rows_with_missing_data (df1)

    if remove_duplicate_columns: # who knows which column it will choose (i.e. column names)
        df = df.T.drop_duplicates().T # transpose, remove duplicate rows, and transpose again
    return df

def calc_freq_N_from_string (genome):
    l = len(genome)
    if (l):
        number_Ns = sum([genome.upper().count(nuc) for nuc in ["N", "-"]])
        return number_Ns / l
    else: return 1.

def calc_freq_ACGT_from_string (genome):
    l = len(genome)
    if (l):
        number_ACGTs = sum([genome.upper().count(nuc) for nuc in ["A", "C", "G", "T"]])
        return number_ACGTs / l
    else: return 0.

def add_sequence_counts_to_metadata (metadata, sequences, from_scratch = None):
    if from_scratch is None: from_scratch = True
    """ counts the proportions of indel/uncertain bases (i.e `N` and `-`) as well as number of ACGT.
    Notice that they do not sum up to one since we have partial uncertainty (`W`,`S`, etc.) that can be used
    This should be done on UNALIGNED sequences
    """
    def calc_freq_N (index):
        if index not in sequences: return 1.
        if sequences[index] is None: return 1. ## missing sequences
        return calc_freq_N_from_string (str(sequences[index].seq)) 
    def calc_freq_ACGT (index):
        if index not in sequences: return 0.
        if sequences[index] is None: return 0. ## missing sequences
        return calc_freq_ACGT_from_string (str(sequences[index].seq)) 

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

def merge_global_and_local_metadata (metadata0, csv0):
    """
    Local (NORW) sequences are named NORW-EXXXX, which is central_sample_id. We match those to COGUK whenever possible.
    This adds all our local info to the database, i.e. not only the ones with corresponding sequence 
    (local csv includes info on 'rejected' or not sequenced samples)
    """
    metadata = metadata0.copy()
    csv = csv0.copy()

    logger.info("Merging global metadata (COGUK+GISAID usually) with local one (NORW)")
    remove_cols = [c for c in remove_from_master_cols if c in csv.columns]  ## imputed values from previous iteration
    if len(remove_cols):
        csv.drop (labels = remove_cols, axis=1, inplace = True) 

    matched = metadata[ metadata["central_sample_id"].isin(csv.index) ] # index of csv is "central_sample_id"
    csv["submission_org_code"] = "NORW"
    csv["submission_org"] = "Norwich"
    
    ## temporarily use central_sample_id as index, so that we can merge_by_index
    matched.reset_index(inplace=True) ## downgrades current index to a regular column
    matched.set_index ("central_sample_id", append=False, drop = True, inplace = True) # append creates a new col w/ index 
    ## merge csv with corresponding elements from global metadata (note that these are just intersection with csv)
    csv = df_merge_metadata_by_index (csv, matched) 
    
    # replace receive leaf names in case it's NORW-E996C 
    csv["peroba_seq_uid"] = csv["peroba_seq_uid"].fillna(csv.index.to_series())
    csv["sequence_name"] = csv["sequence_name"].fillna(csv.index.to_series())
    #csv[metadata.index.names[0]] = csv[metadata.index.names[0]].fillna(csv.index.to_series()) # same as above

    ## revert index to same as global metadata ("peroba_seq_uid" usually)
    csv.reset_index (drop=False, inplace=True) ## drop=True means drop index completely, not even becomes a column
    csv.set_index (metadata.index.names, drop = True, inplace = True) # drop to avoid an extra 'peroba_seq_uid' column
    metadata = df_merge_metadata_by_index (metadata, csv) 
    return metadata, csv

def transparent_cmap (color = None, cmap = None, final_alpha = None):
    # http://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
    if color is None: color = "blue"
    if (final_alpha is None) or (final_alpha < 0.01): final_alpha = 1.
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

def list_from_custom_colorset (start=0, n_base_colours=4, list_size = 8):
    base_set = ["F26440", "D4A0A7", "8FD5A6", "599B6C", "77A4BB", "EFC69B", "B1B7D1", "F6AA1C", # pallete of handcrafted colours
               "D8C99B", "FE938C", "99B2DD", "FFF9A5", "BAD8B6", "888098", "CACCC1", "5385AA", # such that sets of consecutive 
               "C48983", "666822", "E9CEBA", "56AACC", "FFBD44", "729FA4", "8C693C", "B04D2E", # fours should be fine
               "AE8EAF", "86927D", "E6B89C", "fed611", "73bfb8", "C4C8A6", "fef5ae", "BB9457", 
               "9B8870", "744D38", "991E4B", "74929F", "ACE0DC", "ADA982", "7A6157", "A2AABB"]
    n_base_set = len(base_set)
    if list_size < 1: list_size = 1
    if (n_base_colours < 1): n_base_colours = 1 # if you want a single colour?
    if start < 0: start = 0;
    start = start % n_base_set  # cycle 
    end = (start + n_base_colours) % n_base_set
    if start > end:
        pivot = end; end = start; start = pivot

    rgbcolours = [hex2rgb(x) for x in base_set[start:end]] # functions use RGB, not hex
    # create a linear (continuous) colormap and return list with discretised values
    custom_cmap = colors.LinearSegmentedColormap.from_list("custom", rgbcolours)
    if (list_size > n_base_colours):
        return custom_cmap(np.linspace(0, 1, list_size))
    elif list_size == 1: 
        return custom_cmap(0.5) # average colour in basic set
    else: # no deep thought here, just avoiding returning same set 
        cstm = custom_cmap(np.linspace(0, 1, list_size+1))
        return cstm[1:]

def hex2rgb (x, norm = True):
    if norm:
        return tuple(float(int(x[i:i+2], 16))/255. for i in (0, 2, 4))

def remove_duplicated_sequences (sequences): # input is list, returns a dict
    uniq_seqs = {}
    uniq_qual = {}
    duplicates = []
    for x in sequences:
        seq = str(x.seq)
        quality = len(seq) - sum([seq.upper().count(nuc) for nuc in ["N", "-"]])
        if x.id in uniq_seqs.keys(): # sequence name has been seen before 
            if uniq_qual[x.id] < quality: # replaces if better quality
                uniq_qual[x.id] = quality
                uniq_seqs[x.id] = x
            duplicates.append(x.id)
        else: # first time we see this sequence
            uniq_qual[x.id] = quality
            uniq_seqs[x.id] = x
    if len(duplicates)>0:
        logger.warning ("%s duplicate (i.e. with same name) sequences were resolved by choosing the one with highest quality", len(duplicates))
        duplicates = list(set(duplicates))
        logger.debug ("And the sequence names are:\n%s\n", "\n".join(duplicates))
    else:
        logger.info ("All sequences have distinct names")

    return uniq_seqs, uniq_qual

def sequence_dict_pair_with_better_quality (s1, s2, q1=None, q2=None, matched=None, description=["global","local"]):
    if matched is None:
        matched = list(set(s1.keys()) & set(s2.keys()))
    if q1 is None:
        q1 = {x:(len(str(s1[x].seq)) - sum([str(s1[x].seq).upper().count(nuc) for nuc in ["N", "-"]])) for x in matched}
    if q2 is None:
        q2 = {x:(len(str(s2[x].seq)) - sum([str(s2[x].seq).upper().count(nuc) for nuc in ["N", "-"]])) for x in matched}
    better = [[],[]]
    debug_seq = {} ## DEBUG
    for x in matched:
        if q1[x] > q2[x]: ## global is better
            better[0].append(x)
            debug_seq[f"{x}_G_pass_{q1[x]}"] = SeqRecord (Seq.Seq(str(s1[x].seq)), id=x)
            debug_seq[f"{x}_L_fail_{q2[x]}"] = SeqRecord (Seq.Seq(str(s2[x].seq)), id=x)
            s2[x].seq = s1[x].seq
        if q1[x] < q2[x]: ## local is better
            better[1].append(x)
            debug_seq[f"{x}_G_fail_{q1[x]}"] = SeqRecord (Seq.Seq(str(s1[x].seq)), id=x) 
            debug_seq[f"{x}_L_pass_{q2[x]}"] = SeqRecord (Seq.Seq(str(s2[x].seq)), id=x) 
            s1[x].seq = s2[x].seq
    logger.info("{} matched pairs, with {} better on {} and {} better on {} sequences".format(len(matched), 
        len(better[0]), description[0], len(better[1]), description[1]))
    fname = save_sequence_dict_to_file (debug_seq)
    logger.info(f"DEBUG::Quality-based comparison saved to file {fname}")
    return s1, s2, better

def save_sequence_dict_to_file (seqs, fname=None, use_seq_id = False):
    if fname is None: fname = "tmp." + '%012x' % random.randrange(16**12) + ".aln.xz"
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
            if use_seq_id is True: # default is to use dict key
                name = rec.id
            if rec:  ## missing/query sequences
                seq = str(rec.seq)
                fw.write(str(f">{name}\n{seq}\n").encode())
                rec.id = name ## make sure alignment will have same names
    logger.info(f"Finished saving sequences")
    return os.path.basename(fname)

def local_name_run (seq_name):
    x = seq_name.split(".")
    if len(x) > 1: return x[0], int(x[1])
    else:          return x[0], 0
