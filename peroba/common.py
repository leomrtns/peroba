#!/usr/bin/env python

## TODO: "UNKNOWN SOURCE" and "UNKNOWN" are the same adm2 (in cog)

import logging, ete3, argparse
import numpy as np, pandas as pd
from Bio import Seq, SeqIO, Align, AlignIO, Phylo, Alphabet, pairwise2
#from Bio.SeqRecord import SeqRecord
import datetime, time, codecs, sys, gzip, bz2, re, glob, pickle, collections, subprocess, os, errno, random, itertools, pathlib

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
        "sequences": ".sequences.fasta.bz2",
        "alignment": ".sequences.aln.bz2"
        }
        
#asr_cols = ["adm2", "uk_lineage", "lineage", "phylotype", "submission_org_code", "date_sequenced", "source_age", "source_sex", "collecting_org", "ICU_admission"]
asr_cols = ["adm2", "uk_lineage", "lineage", "phylotype", "special_lineage", "adm2_private"]

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

    df1['adm2'] = df1['adm2'].str.title()
    df1 = df1.groupby(df1.index).aggregate("first"); # duplicated indices are allowed in pandas
    return df1

def df_merge_metadata_by_index (df1, df2): # merge by replacing NA whenever possible, using INDEX
    df1 = df1.combine_first(df2)
    df1 = df2.combine_first(df1)
    # if before combine_first() there were dups they're kept by pandas
    df1 = df1.groupby(df1.index).aggregate("first"); 
    df1.sort_index(inplace=True)
    return df1

def df_finalise_metadata (df, exclude_na_rows = None, exclude_columns = "default", remove_duplicate_columns = True):
    ''' if exclude_na_rows is the string "default" then we use sensible choices
    '''
    if exclude_na_rows == "default":
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    # df['days_since_Dec19'] = df['collection_date'].map(lambda a: get_days_since_2019(a, impute = True)) ## not here
    df["collection_datetime"] = pd.to_datetime(df["collection_date"], infer_datetime_format=False, errors='coerce')
    df = df.sort_values(by=['lineage_support', 'collection_datetime'], ascending=[False, True])

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
