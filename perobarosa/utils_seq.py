import os
import logging, xxhash

import numpy as np, pandas as pd
from Bio import Seq, SeqIO
import random, datetime, sys,  re, glob, collections, subprocess, itertools, pathlib, base64
import lzma, gzip, bz2


logger = logging.getLogger(__name__) 
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='fileutils %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

def read_fasta_as_list (filename, fragment_size = 0, check_name = False):
    unaligned = []
    if filename.endswith(".bz2"): this_open = bz2.open #if "bz2" in filename[-5:]: this_open = bz2.open
    if filename.endswith(".gz"):  this_open = gzip.open
    if filename.endswith(".xz"):  this_open = lzma.open
    else:  this_open = open
    with this_open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record.seq  = Seq.Seq(str(record.seq.upper()).replace(".","N")) # one damn sequence has dots 
            if (check_name and "|" in record.description): # consensus from COGUK and GISAID names: `hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020-03-25`
                seqname = record.description.split("|")[0]
                record.id = seqname.replace("hCoV-19/","")
            if len(record.seq) > fragment_size:
                unaligned.append(record)
    logger.info("Read %s sequences from file %s", str(len(unaligned)), filename)
    return unaligned

def align_mafft_in_blocks (sequences, reference_file, seqs_per_block = 2000):    # list not dict
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
    mode = "wb"
    if fname.endswith(".bz2"): this_open = bz2.open #if "bz2" in filename[-5:]: this_open = bz2.open
    if fname.endswith(".gz"):  this_open = gzip.open
    if fname.endswith(".xz"):  this_open = lzma.open
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
