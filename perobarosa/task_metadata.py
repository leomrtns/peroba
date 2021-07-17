#!/usr/bin/env python
from perobarosa.utils_seq import *
import pandas as pd

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_TASK %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

peroba_xtracol = ["freq_ACGT", "freq_N", "seq_hash"]
peroba_columns = ["strain", "gisaid_id", "date", "location", "location_gisaid", "location_coguk", "age", "gender", "gisaid_clade", "pango_lineage", "timestamp"]
gisaid_columns = ["Virus name", "Accession ID", "Collection date", "Location", "Additional location information", "Patient age", "Gender", "Clade", "Pango lineage"]
epidem_columns = ["strain", "gisaid_epi_isl", "date", "age", "sex", "GISAID_clade", "pango_lineage", "region","country","division", "region_exposure","country_exposure","division_exposure"]
coguk_columns = ['sequence_name', 'gisaid_id', 'sample_date', 'country', 'adm1', 'NUTS1', 'source_age', 'source_sex','travel_history', 'lineage']

def metadata (metadata_file, defaults, alignment = None, csvfile = None, output = None, entry_timestamp = None, calc_seqs = False):
    if output is None: output = defaults["current_dir"] + "peroba_meta." +  defaults["timestamp"] + ".tsv.xz"
    if entry_timestamp is None: timestamp = datetime.datetime.now().strftime("%y%m%d")
    else:                       timestamp = str(entry_timestamp)
    has_extra_columns = False

    # read existing metadata file (we assume names are clean, i.e. valid fasta headers and same as in fasta file) 
    csv = None
    if (csvfile is not None):
        logger.info(f"Reading {csvfile} from previous round of peroba metadata")
        csv = pd.read_csv (csvfile, compression="infer", sep="\t", dtype='unicode')
        keep_columns = [x for x in peroba_columns  if x in csv.columns]
        if (len(keep_columns) < len(peroba_columns)): 
            logger.error ("Existing table %s does not look like a peroba metadata file, missing at least one column from\n%s",csvfile, "\n".join(peroba_columns))
            sys.exit(1)
        if ("strain" not in keep_columns):
            logger.error (f"Column 'strain' missing from existing table file ({csvfile} must be from previous call to peroba)")
            sys.exit(1)
        xtra_columns = [x for x in peroba_xtracol if x in csv.columns]
        if (len(xtra_columns) > 0):
            logger.info(f"Good, found extra columns (freq_ACGT etc.) in peroba metadata {csvfile}")
            keep_columns += xtra_columns
            has_extra_columns = True
        csv = csv[keep_columns]  ## reorder s.t. "strain" comes first etc
        csv = csv.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN s.t. new file can overwrite it

    # read new file from GISAID package / GISAID epidem / COGUK
    df = None
    logger.info(f"Reading {metadata_file} with new entries")
    df0 = pd.read_csv (metadata_file, compression="infer", sep="\t", dtype='unicode')

    keep_columns = [x for x in peroba_columns if x in df0.columns]
    if (len(peroba_columns) <= len(keep_columns)): 
        logger.info(f"Assuming {metadata_file} is already in peroba format")
        xtra_columns = [x for x in peroba_xtracol if x in df0.columns]
        if (len(xtra_columns) > 0):
            logger.info(f"Good, found extra columns (freq_ACGT etc.) in new {metadata_file}")
            keep_columns += xtra_columns
            has_extra_columns = True
        df = df0[keep_columns]
        if ("timestamp" in keep_columns): timestamp = None ## prevents overwritting

    keep_columns = [x for x in gisaid_columns if x in df0.columns]
    if (df is None and len(gisaid_columns) - len(keep_columns) < 2): 
        logger.info(f"{metadata_file} looks like GISAID metadata file from default 'package'")
        df = convert_from_gisaid_metadata (df0)
        if (entry_timestamp is None): timestamp = str(timestamp) + ".0"

    keep_columns = [x for x in epidem_columns if x in df0.columns]
    if (df is None and len(epidem_columns) - len(keep_columns) < 2): 
        logger.info(f"{metadata_file} looks like GISAID epidemiology metadata file")
        df = convert_from_epidem_metadata (df0)
        if (entry_timestamp is None): timestamp = str(timestamp) + ".1"

    if (df is None):
        df0 = pd.read_csv (metadata_file, compression="infer", sep=None, engine="python", dtype='unicode') ## COGUK is csv
    keep_columns = [x for x in coguk_columns if x in df0.columns]
    if (df is None and len(coguk_columns) - len(keep_columns) < 2): 
        logger.info(f"{metadata_file} looks like COGUK metadata file")
        df = convert_from_coguk_metadata (df0)
        if (entry_timestamp is None): timestamp = str(timestamp) + ".2"

    if (df is None or df.shape[0] == 0): 
        logger.warning("Could not convert metadata; column with the sequence names not found or unrecognised")
        sys.exit(1)

    df["strain"] = df["strain"].apply(clean_gisaid_name) ##  "value is trying to be set on a copy of a slice" complaint
    if (timestamp is not None):
        df["timestamp"] = timestamp
    logger.info("Read %d rows from table to be added", df.shape[0])

    # merge current and new metadata iff both seq name and gisaid ID are the same 
    if (csv is not None and csv.shape[0] > 0):
        df.set_index(["strain", "gisaid_id"], inplace=True, append=False, drop=False) ## append removes current index (lame counter) 
        csv.set_index(["strain", "gisaid_id"], inplace=True, append=False, drop=False) ## drop removes column, thus drop=F keeps both
        l1 = df.shape[0]
        l2 = csv.shape[0]
        df = csv.combine_first (df) ## df will only be added if absent from csv ## python suggests "sort=False" here
        df.reset_index(drop=True, inplace=True)
        logger.info("Current table has %d rows, merged table will have %d rows, with %d common ones. ", l1, df.shape[0], l1 + l2 - df.shape[0]);

    # check if several rows may have same sequence name
    n_unique = df.shape[0] - len(df["strain"].unique())
    if (n_unique > 0): 
        logger.warning((f"{n_unique} rows have same name; this is bad practice by GISAID but it's not a bug;" 
                "'peroba align' dedups sequences by default, but other software may fail if several sequences have same name"))
        dup_seqs = collections.Counter(df["strain"]).most_common(10)
        logger.info("Examples of duplicate names:\n%s   etc.\n", "\n".join([x[0] for x in dup_seqs]))

    # if metadata has extra columns with sequence stats, then we just update from missing seqs

    if alignment is not None:
        if has_extra_columns:
            missing_freq_set = set(df.loc[df["freq_ACGT"].isnull(), "strain"].unique())
            logger.info("Extra columns found; will calculate seq stats for %d samples", len(missing_freq_set))
        elif calc_seqs is True:
            logger.info("Extra columns not found but will calculate from scratch seq stats")
            missing_freq_set = set(df["strain"].unique())
        else:
            logger.info("Extra columns not found and no seq stats calculation will take place (add '-f' next time to force it)")
            missing_freq_set = None

    # read existing alignments
    aln_seqnames = set()
    if alignment is None: logger.warning (f"No alignment given; will output all entries")
    else:
        xtra_col_df = None
        logger.info (f"Reading alignments may take a while")
        for aln in alignment:
            logger.debug(f"Reading alignment {aln}") 
            seqnames, xtra_df = read_fasta_headers (aln, update_set = missing_freq_set)
            aln_seqnames.update(seqnames)
            if xtra_col_df is None:
                xtra_col_df = xtra_df
            elif xtra_df is not None:
                xtra_col_df = xtra_col_df.combine_first(xtra_df)  ## "strain" column is index of dataframe
        aln_seqnames = set (aln_seqnames) ## much faster lookup than list
        logger.info("Total of %d aligned sequences", len(aln_seqnames))
        if (len(aln_seqnames)):
            logger.warning("Keeping only rows from samples present in alignment")
            df = df[ df["strain"].isin(aln_seqnames) ]
            missing_seqs = set(aln_seqnames - set(df["strain"].unique()))
            if (len(missing_seqs) > 0):
                errfile = defaults["current_dir"] + "peroba_meta-missing_rows." +  defaults["timestamp"] + ".txt"
                logger.warning ("There are %d aligned sequences without metadata. List of sequences being saved to file %s", len(missing_seqs), errfile)
                with open_anyformat(errfile, "w") as fw:
                    for seqname in missing_seqs: 
                        fw.write(str(f"{seqname}\n").encode())
            if xtra_col_df is not None: ## at least one seq has new extra information
                logger.info("Will now merge info about %d new sequences", xtra_col_df.shape[0])
                xtra_col_df.index.name = "strain"
                df.set_index("strain", inplace=True, append=False, drop=False) ## drop removes column, thus drop=F keeps both
                df = df.combine_first (xtra_col_df)  ## obs: indices don't have to be unique (so several rows with same strain are updated)
                df.reset_index(drop=True, inplace=True) ## can drop the index column since it's duplicated with "drop=F" above
        else:
            logger.warning ("No sequence names found in alignments; metadata will have all rows")

    keep_columns = [x for x in peroba_columns + peroba_xtracol if x in df.columns]
    df = df[keep_columns] # reorder columns
    # some freq_ACGT are missing, this astype(int) does not work (NaN is a float)
    if "freq_ACGT" in keep_columns:
        df["freq_ACGT"] = df["freq_ACGT"].astype('float').astype(pd.Int32Dtype()) ## object -> float -> nullable integer (Int, not int)
    if "freq_N" in keep_columns:
        df["freq_N"] = df["freq_N"].astype('float').astype(pd.Int32Dtype())
    logger.info (f"Saving peroba metadata file into {output}")
    df.to_csv (output, sep="\t", index=False)
    return

def convert_from_gisaid_metadata (df0):
    df = df0
    rename_dict = {
            "Virus name":"strain", 
            "Accession ID":"gisaid_id", 
            "Collection date":"date", 
            "Location":"location", 
            "Additional location information":"location_gisaid", 
            "Patient age":"age", 
            "Gender":"gender", 
            "Clade":"gisaid_clade", 
            "Pango lineage":"pango_lineage"
            }
    for old, new in rename_dict.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old:new}, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN 

    ## add missing columns with empty data, to be able to merge
    missing_columns = [x for x in peroba_columns if x not in df.columns]
    for col in missing_columns:
        df[col] = np.nan

    if ("strain" not in df.columns):
        logger.error ("Column 'strain' is missing from metadata file")
        return None # allows for another format conversion to be tried 
    return df[peroba_columns] ## remove other columns

def convert_from_epidem_metadata (df0):
    df = df0
    rename_dict = {
            "gisaid_epi_isl":"gisaid_id", 
            "Location":"location", 
            "sex":"gender", 
            "GISAID_clade":"gisaid_clade" 
            }
    for old, new in rename_dict.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old:new}, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN 
    # create additional location info (if any value is np.nan, then the new cell is np.nan)
    valid_cols = [x for x in ["region","country","division"] if x in df.columns]
    if (len(valid_cols) == 3):
        df["location"] = df["region"] + " / " + df["country"] + " / " + df["division"]
    valid_cols = [x for x in ["region_exposure","country_exposure","division_exposure"] if x in df.columns]
    if (len(valid_cols) == 3):
        df["location_gisaid"] = df["region_exposure"].fillna('') + "/" + df["country_exposure"].fillna('') + "/" + df["division_exposure"].fillna('') 
        df["location_gisaid"] = df["location_gisaid"].replace(r'^\s*/\s*/\s*$', np.nan, regex=True) ## of both travel and nuts are empty 

    ## add missing columns with empty data, to be able to merge
    missing_columns = [x for x in peroba_columns if x not in df.columns]
    for col in missing_columns:
        df[col] = np.nan
    if ("strain" not in df.columns):
        logger.error ("Column 'strain' is missing from metadata file")
        return None # allows for another format conversion to be tried 
    return df[peroba_columns] ## remove other columns

def convert_from_coguk_metadata (df0):
    df = df0
    rename_dict = {
            "sequence_name":"strain", 
            "sample_date":"date", 
            "source_age":"age", 
            "source_sex":"gender", 
            "lineage":"pango_lineage"
            }
    for old, new in rename_dict.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old:new}, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN 
    valid_cols = [x for x in ["country","adm1"] if x in df.columns]
    if (len(valid_cols) == 2):
        df["location"] = df["country"].fillna('') + " / " + df["adm1"].fillna('')
        df["location"] = df["location"].replace(r'^\s*/\s*$', np.nan, regex=True)  
    valid_cols = [x for x in ["travel_history", "NUTS1"] if x in df.columns]
    if (len(valid_cols) == 2):
        df["location_coguk"] = df["travel_history"].fillna('') + "|" + df["NUTS1"].fillna('')
        df["location_coguk"] = df["location_coguk"].replace(r'^\s*\|\s*$', np.nan, regex=True) ## of both travel and nuts are empty 

    ## add missing columns with empty data, to be able to merge
    missing_columns = [x for x in peroba_columns if x not in df.columns]
    for col in missing_columns:
        df[col] = np.nan
    if ("strain" not in df.columns):
        logger.error ("Column 'strain' is missing from metadata file")
        return None # allows for another format conversion to be tried 
    return df[peroba_columns] ## remove other columns

def merge (metadata, alignment, defaults):
    csv_ofile = defaults["current_dir"] + "perobaDB." +  defaults["timestamp"] + ".tsv.xz"
    aln_ofile = defaults["current_dir"] + "perobaDB." +  defaults["timestamp"] + ".aln.xz"
    aln_efile = defaults["current_dir"] + "perobaDB-excluded." +  defaults["timestamp"] + ".aln.xz"
    logger.info(f"Reading metadata file {metadata}")
    csv = pd.read_csv (metadata, compression="infer", sep="\t", dtype='unicode')
    keep_columns = [x for x in peroba_columns  if x in csv.columns]
    if (len(keep_columns) < len(peroba_columns)): 
        logger.error ("Table %s is not a peroba metadata file, missing at least one column from\n %s",metadata, "\n".join(peroba_columns))
        sys.exit(1)
    xtra_columns = [x for x in peroba_xtracol if x in csv.columns]
    if (len(xtra_columns) == 3):
        logger.info(f"Found extra columns (freq_ACGT etc.) in peroba metadata; will check how complete it is")
        missing_freq_set  = set(df.loc[df["freq_ACGT"].isnull(), "strain"].unique())
        missing_freq_set |= set(df.loc[df["freq_N"].isnull(), "strain"].unique())
        missing_freq_set |= set(df.loc[df["seq_hash"].isnull(), "strain"].unique())
    else:
        missing_freq_set = None

    # split table into rows with and without gisaid_id, and merge same ones
    logger.info(f"Removing duplicate/redundant metadata rows")
    old_size = csv.shape[0]
    csv_1 = csv[csv["gisaid_id"].isnull()] 
    csv_2 = csv[~csv["gisaid_id"].isnull()]
    csv_2 = csv_2.groupby("gisaid_id").aggregate("first")
    csv_2.reset_index(drop=False, inplace=True) # gisaid_id becomes a column again (drop=True would delete it)
    csv = pd.concat([csv_1, csv_2])
    csv = csv.groupby("strain").aggregate("first") # new key is now 'strain', with column _lost_ (thus must reset_index(drop=False)
    include_set = set(csv.index.tolist()) ## set lookup much faster than list lookup
    logger.info(f"Metadata size decrease from {old_size} to {csv.shape[0]} rows")


    ofl = open_anyformat (aln_ofile, "w")
    efl = open_anyformat (aln_efile, "w") 
    logger.info(f"Selected sequences will be saved to {aln_ofile} and excluded will be saved to {aln_efile}")
    invalid = None
    for aln in alignment:
        logger.info(f"Reading alignment {aln}") 
        df, inval_df = partition_fasta_by_set (aln, ofl, efl, include_set, update_set = missing_freq_set)
        csv = csv.combine_first (df) ## df will only be added if absent from csv ## python suggests "sort=False" here
        if invalid is None:
            invalid = inval_df
        elif inval_df is not None and len (inval_df):
            invalid = pd.concat([invalid, inval_df])
    ofl.close()
    efl.close()
    
    if (invalid and len(invalid)):
        outcsv = defaults["current_dir"] + "peroba_align-excluded." +  defaults["timestamp"] + ".csv.gz"
        logger.info(f"Saving list of excluded sequences to {outcsv} table") 
        invalid.to_csv (outcsv, index=False)

    csv = csv[~csv["freq_ACGT"].isnull()] ## only rows with sequence information
    csv.reset_index(drop=False, inplace=True)
    csv = csv[peroba_columns + peroba_xtracol] # reorders columns
    logger.info(f"Saving metadata file {csv_ofile}")
    csv.to_csv (csv_ofile, sep="\t", index=False)

def stats (alignment, defaults):
    csv_ofile = defaults["current_dir"] + "peroba_stats." +  defaults["timestamp"] + ".tsv.xz"
    csv = None
    for aln in alignment:
        logger.info(f"Reading alignment {aln}") 
        stat_df = read_fasta_calc_stats_only (aln)
        if csv is None:
            csv = stat_df 
        elif stat_df is not None and len(stat_df):
            csv = pd.concat([csv, stat_df])
    if csv is not None and len(csv):
        logger.info(f"Saving metadata file {csv_ofile}")
        csv.index.name = "strain"
        csv.to_csv (csv_ofile, sep="\t", index=True)
    else:
        logger.warning("No stats were calculated, no sequences found")


