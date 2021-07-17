import os, logging, xxhash, numpy as np, pandas as pd
from Bio import Seq, SeqIO
import random, datetime, sys, re, glob, collections, subprocess, itertools, pathlib, base64
import lzma, gzip, bz2


logger = logging.getLogger(__name__) 
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_UTIL %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

# GISAID sequence names for non-humans start with following prefixes (like hCoV-19/env, hCoV-19/mink etc):
non_human_prefix = ["canine", "leopard", "monkey", "hamster", "gorilla", "mouse", "tiger", "bat", "pangolin", "lion", "dog", "cat", "mink", "env"]

def clean_gisaid_name (description):
    seqname = re.sub("\s+","", description.split("|")[0])
    seqname = seqname.replace("'","-") # Coted-Ivoire
    return seqname.replace("hCoV-19/","")

def get_extra_cols_from_record (seq, trim=500):
    seq = seq.upper()
    vals =  [
            int(sum([seq.count(nuc) for nuc in ["A", "C", "G", "T"]])),
            int(sum([seq.count(nuc) for nuc in ["N", "-"]])),
            str(xxhash.xxh32_hexdigest(str(seq[trim:-trim]))) # alias of xxh32().hexdigest()
            ]
    return vals

def read_fasta_calc_stats_only (filename, check_name = True, trim=500):
    n_valid = 0
    seq_info = dict()

    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                record.id = clean_gisaid_name (record.description) 
            seq_info[record.id] = get_extra_cols_from_record (record.seq, trim)
            n_valid += 1
            if (not n_valid%100000):  logger.info(f"Stats about {n_valid} sequences calculated")
    logger.info("Calculated %s sequence stats from file %s", str(len(seq_info)), filename)
    if not seq_info:
        return None
    else:
        return pd.DataFrame.from_dict (seq_info, orient='index', columns = ["freq_ACGT", "freq_N", "seq_hash"])

def read_fasta_headers (filename, check_name = True, trim=500, update_set = None):
    seqnames = []
    n_valid = 0
    if update_set is not None: seq_info = dict()  ## default is to calculate only for list (used also by align tasks, which don't care about table)

    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                record.id = clean_gisaid_name (record.description) 
            seqnames.append(record.id)
            if update_set is not None and record.id in update_set:
                seq_info[record.id] = get_extra_cols_from_record (record.seq, trim)
                n_valid += 1
                if (not n_valid%100000):  logger.info(f"Stats about {n_valid} sequences calculated")

    logger.info("Read %s sequence names from file %s", str(len(seqnames)), filename)
    if update_set is None or not seq_info:
        info_df = None
    else:
        logger.info("Calculated %s sequence stats from file %s", str(len(seq_info)), filename)
        info_df = pd.DataFrame.from_dict (seq_info, orient='index', columns = ["freq_ACGT", "freq_N", "seq_hash"])
    return seqnames, info_df 

def read_fasta_as_list (filename, fragment_size = 0, check_name = False):
    unaligned = []
    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record.seq  = Seq.Seq(str(record.seq.upper()).replace(".","N")) # one damn sequence has dots 
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                record.id = clean_gisaid_name (record.description) 
            if len(record.seq) > fragment_size:
                unaligned.append(record)
    logger.info("Read %s sequences from file %s", str(len(unaligned)), filename)
    return unaligned

def partition_fasta_by_set (infile, ofl, efl, include_set, trim=500, update_set = None):
    seq_info = dict() ## default is to calculate for all values
    invalid = {"taxon":[],"excluded":[]}
    n_valid = 0
    with open_anyformat (infile, "r") as ifl:
        for record in SeqIO.parse(ifl, "fasta"):
            if ("|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                record.id = clean_gisaid_name (record.description) 
            if record.id in include_set:
                n_valid += 1
                if (not n_valid%250000):  logger.info(f"{n_valid} sequences saved")
                if update_set is None or record.id in update_set:
                    seq_info[record.id] = get_extra_cols_from_record (record.seq, trim)
                ofl.write(str(f">{record.id}\n{record.seq}\n").encode())
            else:
                invalid["taxon"].append(record.id)
                invalid["excluded"].append(f"duplicate")
                efl.write(str(f">{record.id}\n{record.seq}\n").encode())
    logger.info(f"{n_valid} sequences saved from file")
    info_df = pd.DataFrame.from_dict (seq_info, orient='index', columns = ["freq_ACGT", "freq_N", "seq_hash"])
    return info_df, pd.DataFrame(invalid)

def read_fasta_new_only (fastafile, prev_seqnames = [], min_length = 10000, ambig = 0.8, check_name = False):
    sequences = []
    invalid = {"taxon":[],"excluded":[]}
    n_valid = 0
    with open_anyformat (fastafile, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                record.id = clean_gisaid_name (record.description) 
            if record.id not in prev_seqnames:
                this_len = len(record.seq)
                if this_len > min_length:
                    record.seq = record.seq.upper()
                    freq_ACGT = sum([record.seq.count(nuc) for nuc in ["A", "C", "G", "T"]]) / this_len 
                    if (1 - freq_ACGT) < ambig:
                        #record.seq  = Seq.Seq(this_seq)
                        sequences.append(record)
                        n_valid += 1
                        if (not n_valid%10000): 
                            logger.info(f"{n_valid} sequences read")
                    else:
                        invalid["taxon"].append(record.id)
                        invalid["excluded"].append(f"acgt={freq_ACGT}")
                        logger.debug (f"Sequence {record.id} has {freq_ACGT} of ambiguous (non-ACGT) sites")
                else:
                    invalid["taxon"].append(record.id)
                    invalid["excluded"].append(f"len={this_len}")
                    logger.debug (f"Sequence {record.id} too short, has only {this_len} sites")
    if len (invalid["taxon"]): logger.warning ("Number of sequences excluded due to short length or highly ambiguous: %s", len(invalid["taxon"]))
    return sequences, pd.DataFrame(invalid)

def align_mafft_in_blocks (sequences, reference_file, seqs_per_block = 2000):    # list not dict
    if seqs_per_block < 100:   seqs_per_block = 100
    if seqs_per_block > 50000: seqs_per_block = 50000 # mafft chokes on large matrices (but 10k is fine btw)
    nseqs = len(sequences)
    aligned = []
    for i in range (0, nseqs, seqs_per_block):
        last = i + seqs_per_block
        if last > nseqs: last = nseqs
        aln = mafft_align_seqs (sequences[i:last], reference_file = reference_file)
        aligned += aln
        logger.info (f"First {last} sequences aligned")
    return aligned

def mafft_align_seqs (sequences=None, infile = None, outfile = None, reference_file = None, prefix = "/tmp/", exclude_reference = True):    # list not dict
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me a fasta object or file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = f"{prefix}/mafft.fasta"
    else: ifl = infile # if both infile and sequences are present, it will save (overwrite) infile
    if outfile is None: ofl = f"{prefix}/mafft.aln"
    else: ofl = outfile # in this case it will not exclude_reference
    if sequences: SeqIO.write(sequences, ifl, "fasta") ## else it should be present in infile

    reference_file = os.path.expanduser(reference_file) 
    runstr = f"mafft --auto --keeplength --thread -1 --addfragments {ifl} {reference_file} > {ofl}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    aligned = AlignIO.read(ofl, "fasta")

    if infile is None:  os.system("rm -f " + ifl)
    if outfile is None: os.system("rm -f " + ofl)

    if exclude_reference:
        refseq = AlignIO.read(reference_file, "fasta")
        refseqname = refseq[0].id
        aligned = [x for x in aligned if x.id != refseqname]
        if outfile is not None:
            SeqIO.write(aligned, ofl, "fasta")
    return aligned

def uvaia_align_seqs (sequences=None, ambiguous=None, infile = None, outfile = None, reference_file = None, prefix = "/tmp/"):    # list not dict
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me a fasta object or file")
    uniq = '%012x' % random.randrange(16**12)
    if prefix is None: prefix = "./"
    if infile is None: ifl = f"{prefix}{uniq}.fasta.gz"  # save as gz since uvaia can handle it
    else: ifl = infile # if both infile and sequences are present, it will save (overwrite) infile
    if outfile is None or outfile == "FILE": ofl = f"{prefix}{uniq}uvaia.aln" # create file but not delete it, instead returning it
    else: ofl = outfile 
    if (ambiguous is None): ambiguous = 0.1

    with open_anyformat (ifl, "w") as fw: 
        for x in sequences:
            fw.write(str(f">{x.id}\n{x.seq}\n").encode())

    runstr = f"uvaialign -a {ambiguous} -r {reference} {ifl} > {ofl}" # ~/bin hardcoded since climb crontab has issues
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    aligned = AlignIO.read(ofl, "fasta")

    if infile is None:  os.system("rm -f " + ifl)
    if outfile is None: os.system("rm -f " + ofl)
    if outfile == "FILE": return ofl
    else:                 return aligned

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

def remove_duplicated_sequences_list (sequences): # input is list, returns a dict
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
        logger.info ("Checked for duplicates but all sequences have distinct names")

    return uniq_seqs, uniq_qual

def save_sequence_dict_to_file (seqs, fname=None, use_seq_id = False):
    if fname is None: fname = "tmp." + '%012x' % random.randrange(16**12) + ".aln.xz"
    logger.info(f"Saving sequences to file {fname}")
    with open_anyformat (fname, "w") as fw: 
        for name, rec in seqs.items():
            if use_seq_id is True: # default is to use dict key
                name = rec.id
            if rec:  ## missing/query sequences
                seq = str(rec.seq)
                fw.write(str(f">{name}\n{seq}\n").encode())
                rec.id = name ## make sure alignment will have same names
    logger.info(f"Finished saving sequences")
    return os.path.basename(fname)

def open_anyformat (fname, mode = "r"):
    if (mode == "r"): openmode = "rt"
    else:             openmode = "wb"
    if   fname.endswith(".bz2"): this_open = bz2.open #if "bz2" in filename[-5:]: this_open = bz2.open
    elif fname.endswith(".gz"):  this_open = gzip.open
    elif fname.endswith(".xz"):  this_open = lzma.open
    else:  
        this_open = open
#      if (mode == "w"): openmode = "w"  ## only raw file for writting doesn't need "wb"
    return this_open (fname, openmode) 

