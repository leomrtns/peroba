import os, glob, gzip, pickle, subprocess
from Bio import Seq, AlignIO
from Bio.SeqRecord import SeqRecord
import pandas as pd, datetime, re, gspread
from scipy import stats
from oauth2client.service_account import ServiceAccountCredentials


upload_fasta_dirs = "/cephfs/covid/bham/climb-covid19-alikhann/upload/*"
civet_metadata = "/cephfs/covid/bham/results/phylogenetics/latest/civet/cog/cog_global*metadata.csv"
civet_aln = "/cephfs/covid/bham/results/phylogenetics/latest/civet/cog/cog_global*alignment.fasta"
local_dir = "/cephfs/covid/bham/climb-covid19-deoliveiramartinsl/local/peroba_uvaia/" 
ref_wuhan  = local_dir + "MN908947.3.fas" 
merged_fasta  = local_dir + "norw.fasta.gz" 
recent_secret = local_dir + "recent.p" 
uv_aln = local_dir + "uv.aln" 
uv_tbl = local_dir + "uv.tbl" 
merged_tbl_filename = local_dir + "merged.csv.gz"
google_credentials = local_dir + ".credentials.json"

cols_in_fulltab = ["phylotype","uk_lineage", "lineage", "adm1"] # pd.Series.mode can't handle only NaN (unless dropna=False)

def save_fasta_recent (recent_fas, outfile):
    with gzip.open(outfile,"wb") as fw:
        for sequence in recent_fas: 
            fw.write(str(f">{sequence.id}\n{sequence.seq}\n").encode())

def load_fasta_from_recent_dirs (recent_dirs, civet_names):
    if recent_dirs is None: return None
    aligned = []
    old_names = []
    for d in recent_dirs:
        seqname_long = glob.glob(d + "/*") # list, each value has whole path
        seqname = [os.path.basename(x) for x in seqname_long]
        filname = [glob.glob(x + "/*.fa*")[0] for x in seqname_long]
        # exclude recent sequences that are already on civet
        new_names = [[sname, fname] for sname, fname in zip (seqname, filname) if sname not in civet_names]
        for sname, fname in new_names:
            singleseq = AlignIO.read(fname, "fasta")[0] # 
            aligned.append (SeqRecord(singleseq.seq, id=sname, description=sname))
        # recent sequences that are already on civet
        old_names += [sname for sname in seqname if sname in civet_names]
    return aligned, old_names

def update_recent_fasta (recent_dirs, outfile, civet_names):
    recent_fas, old_names = load_fasta_from_recent_dirs (recent_dirs, civet_names)
    if (len(recent_fas) > 0):
        save_fasta_recent (recent_fas, outfile)
    #if filecmp.cmp (newfile, outfile):   ## TRUE if equal
    return old_names, len(recent_fas)

def get_recent_norw_civet_files (norw, civet, n_dirs):
    if n_dirs < 1: n_dirs = 1
    recent_dirs = sorted (glob.glob(norw), key=os.path.getmtime)[-n_dirs:]
    recent_civet_meta = glob.glob (civet)[0]
    return [recent_dirs, recent_civet_meta]

def check_if_nothing_changed (recent_nc_files, secret):
    try:
        saved_nc_files = pickle.load(open(secret, "rb"))
    except:
        print ("falied to open secret file; will start from scratch")
        return False
    nothing_changed = [str(s1) == str(s2) for s1,s2 in zip (recent_nc_files[0], saved_nc_files[0])]
    if (all(nothing_changed)): nothing_changed = (str(recent_nc_files[1]) == str(saved_nc_files[1]))
    if (nothing_changed): return True
    return False

def uvaia (reference, merged, uv_align, civet_align, uv_table):
    runstr = f"uvaialign -r {reference} {merged} > {uv_align}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    runstr = f"uvaia -t 16 -n 1 -r {civet_align} {uv_align} > {uv_table}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)

def update_peroba_sheet (credentials, merged_df):
    scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(credentials, scope)
    client = gspread.authorize(creds)
    # 1. first sheet with lineages only
    sheet = client.open("peroba").sheet1
    full_tab = merged_df.reset_index (drop=False) # groupby creates index, we don't need it here 
    full_tab = full_tab.fillna("") # google doesnt like NaN but we do (for combine_first, for instance)
    sheet.update([full_tab.columns.values.tolist()] + full_tab.values.tolist())

    # 2. get master sheet (and add index to combine_first() ) 
    master_tbl = pd.DataFrame(client.open("SARCOV2-Metadata").sheet1.get_all_records())
    master_tbl.set_index("central_sample_id", inplace=True)

    # 3. combine tables with modified uvaia column names
    cols_to_rename = {x:"uvaia_"+x for x in ["phylotype","uk_lineage", "lineage", "adm1"]}
    merged_df.rename(columns=cols_to_rename, inplace=True)
    #merged_df.rename(columns={"acgt_distance":"uvaia_closest_distance", "row_class":"uvaia_group"}, inplace=True)
    merged_df.rename(columns={"snp_distance":"uvaia_closest_distance", "row_class":"uvaia_group"}, inplace=True)
    merged_df = merged_df.combine_first (master_tbl) 
    merged_df = merged_df.fillna("") # google sheet doesnt like NaN and it cannot handle the index
    merged_df.reset_index (drop=False, inplace=True)
    sheet3 = client.open("peroba").worksheet("merged_with_master") # open second sheet and save from scratch
    sheet3.update([merged_df.columns.values.tolist()] + merged_df.values.tolist())

def get_list_of_civet_seqnames (cfile): ## also returns metadata 
    cols_to_remove = ["secondary_identifier","is_community","is_hcw","is_travel_history","travel_history", "epi_week","is_surveillance"]
    civ_tab = pd.read_csv (glob.glob(cfile)[0], compression="infer", dtype="unicode")
    civ_tab = civ_tab.drop(cols_to_remove, axis=1)
    return civ_tab["central_sample_id"].tolist(), civ_tab

def merge_civet_uvaia (civ_tab, ufile):
    uva_tab = pd.read_csv (ufile, compression="infer", dtype="unicode") # index_col="peroba_seq_uid"
    uva_tab.columns = uva_tab.columns.str.strip()  ## removes tabs/spaces uvaia adds for pretty stdout formatting
    uva_tab = uva_tab.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    for col in ["valid_sites","ACGT_matches","prop_char_matches","partial_matches"]:
        uva_tab[col] = pd.to_numeric(uva_tab[col])
    
    # some query sequences are known to civet (called "old_names" outside this function):
    known_seqs = list(set(uva_tab["query sequence"]) & set(civ_tab["central_sample_id"]))
    # now uva_tab has only unknown query seqs
    uva_tab = uva_tab[~uva_tab["query sequence"].isin(known_seqs)]
    uva_tab = uva_tab.merge(civ_tab, left_on="reference sequence", right_on="sequence_name", how="left")
    # update uva_tab with phylotype et al info
    uva_tab["acgt_distance"] = uva_tab["valid_sites"] - uva_tab["ACGT_matches"]
    uva_tab["snp_distance"] = round (uva_tab["valid_sites"] * (1.0 - uva_tab["ACGT_matches"]),1)
    for col in cols_in_fulltab:
        uva_tab[col].fillna("",inplace=True)
    # one row per new query sequence
    full_tab = uva_tab.groupby('query sequence').agg(acgt_distance=("acgt_distance", pd.Series.mean), snp_distance=("snp_distance", pd.Series.mean))
    for col in cols_in_fulltab:
        df = uva_tab.groupby('query sequence').agg(cname=(col, lambda x: stats.mode(x.tolist())[0]))
        df.rename(columns={"cname":col}, inplace=True)
        full_tab = full_tab.combine_first(df)
    full_tab.reset_index(drop=False, inplace=True)
    full_tab.rename(columns={"query sequence":"central_sample_id"}, inplace=True)
    full_tab["row_class"] = "uvaia"
    ## 2. recent sequences which are already on civet DB (should have been handled before running uvaia, btw) 
    df = civ_tab[civ_tab["central_sample_id"].isin(known_seqs)]
    df = df[cols_in_fulltab + ["central_sample_id"]]
    df["row_class"] = "civet"
    df["acgt_distance"] = 0
    df["snp_distance"] = 0
    full_tab = pd.concat([full_tab,df])
    return full_tab  ## does not have index (only after groupby)
    
def update_merged_civet (full_tab, civ_tab, ofile, old_names):
    ## 1. recent sequences, but which are already on civet DB 
    df = civ_tab[civ_tab["central_sample_id"].isin(old_names)]
    df = df[cols_in_fulltab + ["central_sample_id"]]
    df["row_class"] = "civet"
    df["acgt_distance"] = 0
    df["snp_distance"] = 0
    if full_tab is None: # in case uvaia wasn't needed 
        full_tab = df
    else:
        full_tab = pd.concat([full_tab,df])
    ## 2. sequences that are local to Norfolk region
    local_seqs = civ_tab.loc[(civ_tab["adm2"].str.contains("Norfolk",flags=re.IGNORECASE, na=False)),"central_sample_id"].tolist()
    local_seqs += civ_tab.loc[(civ_tab["central_sample_id"].str.contains("NORW",flags=re.IGNORECASE, na=False)),"central_sample_id"].tolist()
    df = civ_tab[civ_tab["central_sample_id"].isin(local_seqs)]
    df = df[cols_in_fulltab + ["central_sample_id"]]
    df["row_class"] = "local"
    df["acgt_distance"] = 0
    df["snp_distance"] = 0
    full_tab = pd.concat([full_tab,df])
    
    full_tab = full_tab.groupby("central_sample_id").aggregate("first") # also makes central_sample_id the index
    full_tab.to_csv(ofile) # using index from groupby to avoid default numeric one
    return full_tab

def run_all ():
    recent_nc_files = get_recent_norw_civet_files (norw = upload_fasta_dirs, civet = civet_metadata, n_dirs = 3) 
    if (check_if_nothing_changed (recent_nc_files, recent_secret)): return 
    list_civet_names, civ_table = get_list_of_civet_seqnames (civet_metadata)
    old_names, aln_length = update_recent_fasta (recent_dirs = recent_nc_files[0], outfile = merged_fasta, civet_names = list_civet_names)
    merged_df = None
    if (aln_length):
        uvaia (ref_wuhan, merged_fasta, uv_aln, civet_aln, uv_tbl)
        merged_df = merge_civet_uvaia (civ_table, uv_tbl)
    merged_df = update_merged_civet (merged_df, civ_table, merged_tbl_filename, old_names)
    #update_peroba_sheet (google_credentials, merged_df)
    # only save dump if everything went smooth
    pickle.dump (recent_nc_files, open(recent_secret, "wb"))

run_all()
