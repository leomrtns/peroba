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

peroba_columns = ["strain", "gisaid_id", "date", "location", "location_gisaid", "location_coguk", "age", "gender", "gisaid_clade", "pango_lineage", "timestamp"]
gisaid_columns = ["Virus name", "Accession ID", "Collection date", "Location", "Additional location information", "Patient age", "Gender", "Clade", "Pango lineage"]
epidem_columns = ["strain", "gisaid_epi_isl", "date", "age", "sex", "GISAID_clade", "pango_lineage", "region","country","division", "region_exposure","country_exposure","division_exposure"]
coguk_columns = ['sequence_name', 'gisaid_id', 'sample_date', 'country', 'adm1', 'NUTS1', 'source_age', 'source_sex','travel_history', 'lineage']

def update_metadata (metadata_file, defaults, alignment = None, csvfile = None, output = None, entry_timestamp = None):
    if output is None: output = defaults["current_dir"] + "peroba_meta." +  defaults["timestamp"] + ".tsv.xz"
    if entry_timestamp is None: timestamp = datetime.datetime.now().strftime("%y%m%d")
    else:                       timestamp = str(entry_timestamp)
  
    # read existing alignments
    aln_seqnames = set()
    if alignment is None: logger.info (f"No alignment given; will output all entries")
    else:                 
        logger.info (f"Will output only entries present in alignments")
        for aln in alignment:
            logger.debug(f"Reading alignment {aln}") 
            aln_seqnames.update(read_fasta_headers (aln))
        aln_seqnames = set (aln_seqnames) ## much faster lookup than list
        logger.info("Total of %d aligned sequences", len(aln_seqnames))

    # read existing metadata file (we assume names are clean, i.e. valid fasta headers and same as in fasta file) 
    csv = None
    if (csvfile is not None):
        csv = pd.read_csv (csvfile, compression="infer", sep="\t", dtype='unicode')
        keep_columns = [x for x in peroba_columns  if x in csv.columns]
        if (len(keep_columns) < len(peroba_columns)): 
            logger.error ("Existing table %s does not look like a peroba metadata file, missing at least one column from\n %s",csvfile, "\n".join(gisaid_keep_columns))
            sys.exit(1)
        if ("strain" not in keep_columns):
            logger.error (f"Column 'strain' missing from existing table file ({csvfile} must be from previous call to peroba)")
            sys.exit(1)
        csv = csv[keep_columns]  ## reorder s.t. "strain" comes first etc
        csv = csv.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN s.t. new file can overwrite it

    # read new file from GISAID package / GISAID epidem / COGUK
    df = None
    df0 = pd.read_csv (metadata_file, compression="infer", sep="\t", dtype='unicode')

    keep_columns = [x for x in peroba_columns if x in df0.columns]
    if (len(peroba_columns) <= len(keep_columns)): 
        logger.info(f"Assuming {metadata_file} is already in peroba format")
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
    logger.info("Read %d rows from new table", df.shape[0])

    # merge current and new metadata iff both seq name and gisaid ID are the same 
    if (csv is not None and csv.shape[0] > 0):
        df.set_index(["strain", "gisaid_id"], inplace=True, append=False, drop=False) ## append removes current index (lame counter) 
        csv.set_index(["strain", "gisaid_id"], inplace=True, append=False, drop=False) ## drop removes column, thus drop=F keeps both
        l1 = df.shape[0]
        l2 = csv.shape[0]
        df = csv.combine_first (df) ## df will only be added if absent from csv
        df.reset_index(drop=True, inplace=True)
        logger.info("Current table has %d rows, new table will have %d rows, with %d common ones. ", l1, df.shape[0], l1 + l2 - df.shape[0]);

    # check if several rows may have same sequence name
    n_unique = df.shape[0] - len(df["strain"].unique())
    if (n_unique > 0): 
        logger.warning((f"{n_unique} rows have same name; this is bad practice by GISAID but it's not a bug;" 
                "'peroba align' dedups sequences by default, but other software may fail if several sequences have same name"))
        dup_seqs = collections.Counter(df["strain"]).most_common(10)
        logger.info("Examples of duplicate names:\n%s   etc.\n", "\n".join([x[0] for x in dup_seqs]))

    if (len(aln_seqnames)):
        logger.info("keeping only rows from samples present in alignment")
        df = df[ df["strain"].isin(aln_seqnames) ]
        if (df.shape[0] < len(aln_seqnames)):
            errfile = defaults["current_dir"] + "peroba_meta-missing_rows." +  defaults["timestamp"] + ".txt"
            with open_anyformat(errfile, "w") as fw:
                for seqname in set(aln_seqnames - set(df["strain"].unique())):
                    fw.write(str(f"{seqname}\n").encode())
            logger.warning ("Some sequences don't have metadata: %d sequences but only %d table rows. List of sequences with missing metadata on file %s", 
                    len(aln_seqnames), df.shape[0], errfile)


    keep_columns = [x for x in peroba_columns if x in df.columns]
    df = df[keep_columns] # reorder columns
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
        
        
def update_metadata_BKP (metadata_file, defaults, alignment = None, csvfile = None, output = None, entry_timestamp = None):
    if output is None: output = defaults["current_dir"] + "gisaid_meta." +  defaults["timestamp"] + ".tsv.gz"
    if entry_timestamp is None: entry_timestamp = datetime.datetime.now().strftime("%y%m%d")
  
    # read existing alignments
    aln_seqnames = set()
    if alignment is None: logger.info (f"No alignment given; will output all entries")
    else:                 
        logger.info (f"Will output only entries present in alignments")
        for aln in alignment:
            logger.debug(f"Reading alignment {aln}") 
            aln_seqnames.update(read_fasta_headers (aln))
        aln_seqnames = set (aln_seqnames) ## much faster lookup than list
        logger.info("Total of %d aligned sequences", len(aln_seqnames))

    # read existing metadata file (we assume names are clean, i.e. valid fasta headers and same as in fasta file) 
    if (csvfile is not None):
        csv = pd.read_csv (csvfile, compression="infer", sep="\t", dtype='unicode')
        keep_columns = [x for x in gisaid_keep_columns if x in csv.columns]
        if (len(keep_columns) < len(gisaid_keep_columns)): 
            logger.error ("Existing table %s does not look like an updated metadata file, missing at least one column from\n %s",csvfile, "\n".join(gisaid_keep_columns))
            sys.exit(1)
        if ("Virus name" not in keep_columns):
            logger.error (f"Column 'Virus name' missing from existing table file ({csvfile} must be from previous call of peroba)")
            sys.exit(1)
        csv = csv[keep_columns]  ## reorder s.t. Virus name comes first etc
        csv = csv.replace(r'^\s*$', np.nan, regex=True) ## replace empty strings for NaN s.t. new file can overwrite it

    # read new file from GISAID
    df = pd.read_csv (metadata_file, compression="infer", sep="\t", dtype='unicode')
    keep_columns = [x for x in gisaid_keep_columns if x in df.columns]
    if (len(keep_columns) < 2): 
        logger.warning(f"{metadata_file} does not look like GISAID metadata file from default 'package'")
        df, keep_columns = convert_from_gisaid_epidemiology_metadata (df)
        if (len(keep_columns) < 2): 
            logger.error ("Conversion did not work, I do not recognise this metadata file")
            sys.exit(1)
    if ("Virus name" not in keep_columns):
        logger.error ("Column 'Virus name' missing from metadata file, probably not the default GISAID tsv")
        sys.exit(1)
    df = df[keep_columns]
    df["Virus name"] = df["Virus name"].apply(clean_gisaid_name)
    df["peroba_timestamp"] = entry_timestamp

    # merge current and new metadata iff both seq name and gisaid ID are the same 
    if (csv.shape[0] > 0):
        df.set_index(["Virus name", "Accession ID"], inplace=True, append=False, drop=False) ## append removes current index (lame counter) 
        csv.set_index(["Virus name", "Accession ID"], inplace=True, append=False, drop=False) ## drop removes column, thus drop=F keeps both
        l1 = df.shape[0]
        l2 = csv.shape[0]
        df = csv.combine_first (df) ## df will only be added if absent from csv
        df.reset_index(drop=True, inplace=True)
        logger.info("Current table has %d rows, new GISAID table has %d rows, with %d common ones", l1, l2, l1 + l2 - df.shape[0]);

    # check if several rows may have same sequence name
    n_unique = df.shape[0] - len(df["Virus name"].unique())
    if (n_unique > 0): 
        logger.warning((f"{n_unique} rows have same name, bad practice by GISAID but not a bug;" 
                "fasta file may cause other software to complain (although 'peroba align' dedups by default)"))
        dup_seqs = collections.Counter(df["Virus name"]).most_common(10)
        logger.info("Example of duplicate names:\n%s   etc.\n", "\n".join([x[0] for x in dup_seqs]))

    if (len(aln_seqnames)):
        logger.info("keeping only rows from samples present in alignment")
        df = df[ df["Virus name"].isin(aln_seqnames) ]
        if (df.shape[0] < len(aln_seqnames)):
            errfile = defaults["current_dir"] + "gisaid_meta-missing_rows." +  defaults["timestamp"] + ".txt"
            with open_anyformat(errfile, "w") as fw:
                for seqname in set(aln_seqnames - set(df["Virus name"].unique())):
                    fw.write(str(f"{seqname}\n").encode())
            logger.warning ("Some sequences don't have metadata: %d sequences but only %d table rows. Sequences with missing metadata on file %s", 
                    len(aln_seqnames), df.shape[0], errfile)


    keep_columns = [x for x in gisaid_keep_columns if x in df.columns]
    df = df[keep_columns] # reorder columns
    df.to_csv (output, sep="\t", index=False)
    return

def convert_from_gisaid_epidemiology_metadata_BKP (df):
    rename_dict = {
            "strain":"Virus name", 
            "gisaid_epi_isl":"Accession ID",
            "date":"Collection date", 
            "age":"Patient age", 
            "sex":"Gender", 
            "GISAID_clade":"Clade", 
            "pango_lineage":"Pango lineage"
            }
    logger.warning ("Will convert columns, assuming metadata comes from GISAID's epidemiology file");
    for old, new in rename_dict.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old:new}, inplace=True)
    valid_cols = [x for x in ["region","country","division"] if x in df.columns]
    if (len(valid_cols) == 3):
        df["Location"] = df["region"] + "/" + df["country"] + "/" + df["division"]
    valid_cols = [x for x in ["region_exposure","country_exporsure","division_exposure"] if x in df.columns]
    if (len(valid_cols) == 3):
        df["Additional location information"] = df["region_exposure"] + "/" + df["country_exporsure"] + "/" + df["division_exposure"]

    # extra columns (like "country" etc.) will be removed later, by df=df[keep_columns]
    keep_columns = [x for x in gisaid_keep_columns if x in df.columns]
    return df, keep_columns 

