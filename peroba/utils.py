#!/usr/bin/env python
from Bio import Seq, SeqIO, Align, AlignIO, Alphabet 
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
import numpy as np
import datetime, sys, gzip, lzma, bz2, re, subprocess, os, itertools, ete3

#from Bio.Phylo import draw, TreeConstruction  #TreeConstruction.DistanceCalculator, TreeConstruction.DistanceTreeConstructor
#from Bio import Phylo, pairwise2
#from Bio.Align import Applications
#from Bio.Blast import NCBIXML
#import seaborn as sns, pandas as pd
#from sklearn import manifold, metrics, cluster, neighbors, decomposition, preprocessing
#import skbio, parasail, dendropy, time, codecs, joypy, pathlib, random, errno, collections, pickle, glob

def print_redblack(textr="", textb=""):
    print ('\x1b[0;1;31;1m'+ str(textr) + '\x1b[0;1;30;1m'+ str(textb) + '\x1b[0m')

def viralMSA (sequences, reference="MT072688", filename=None): ## alt reference is "RAMBAULT" 
    ifile = "/tmp/in.fas"
    odir  = "/tmp/out"
    ofile = "/tmp/out/in.fas.aln" # output dir + file.aln
    SeqIO.write(sequences, ifile, "fasta")
    runstr = "ViralMSA.py -e leo -s " + ifile + " -r " + reference + " -o " + odir
    print (runstr)
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    
    aligned = AlignIO.read(ofile, "fasta")
    os.system("rm -f " + ifile)
    if (filename):
        os.system("mv " + ofile + " " + filename + "; bzip2 -f " + filename) ## always save in bzip2, overwrite unzipped ones 
    os.system("rm -rf " + odir)
    return aligned

def get_days_since_2019 (x, impute = False):
    if x == np.nan or not x or x == "None" or len(str(x)) < 4 or "XX" in str(x):
        return np.nan
    if impute:  # TODO: split("-")
        if len(str(x)) < 7:
            x = x[:4] + "-01-01"
        elif len(str(x)) < 9:
            x = x[:7] + "-01"
    return (datetime.datetime.strptime(str(x),"%Y-%m-%d") - datetime.datetime(2019, 12, 1)).days

def align_seqs (sequences=None, infile=None, outfile=None):
    print ("started aligning...", flush=True, end=" ")
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me an alignment object or file")
        return [] ## OTOH if both are present then infile is overwritten with contents of sequences[]
    if infile is None:
        ifl = "/tmp/unaligned.fas"
    else:
        ifl = infile
    if outfile is None:
        ofl = "/tmp/mafft.aln"
    else:
        ofl = outfile
    SeqIO.write(sequences, ifl, "fasta")
    proc_run = subprocess.check_output("mafft --ep 0.3 --op 3.0 --auto " + ifl + " > " + ofl, shell=True, universal_newlines=True)
    aligned = AlignIO.read(ofl, "fasta")
    print ("Finished",flush=True)
    if infile is None:
        os.system("rm -f " + ifl)
    if outfile is None:
        os.system("rm -f " + ofl)
    return aligned

def lowlevel_parse_cigar (s):
    out = list(); ind = len(s)-1
    while ind >= 0:
        let = s[ind]; ind -= 1; num = ''
        while s[ind] not in {'M','D','I','S','H','=','X'}:  num += s[ind]; ind -= 1
        out.append((let, int(num[::-1])))
    return out[::-1]

def minimap2_align_seqs (sequences=None, infile = None, reference_path = None, prefix = "/tmp/", n_threads = 4):
    """ todo: we dont need the ref genome .mmi, we can use the ref genome directly, and use seqio.write() """
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me a fasta object or file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = f"{prefix}/minimap2.fasta"
    if sequences:      SeqIO.write(sequences, ifl, "fasta") ## else it should be present in infile
    samfile = f"{prefix}/minimap2.sam"

    reference_path = os.path.expanduser(reference_path) ## FIXME: I should always do this to avoid tilde curse
    index_path = '%s.mmi' % reference_path
    if not os.path.isfile(index_path):
        runstr = f"minimap2 -t {n_threads} -d {index_path} {reference_path}"
        proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    

    runstr = f"minimap2 -ax asm5 -t {n_threads} -o {samfile} {index_path} {ifl}" # preset for div < 1%
    #runstr = f"minimap2 -a -t {n_threads} -o {prefix}/minimap2.sam -d {index_path} {ifl}" # no presets
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    
    ## SAM to FASTA (by https://github.com/niemasd)
    ref_seq = []
    for line in open(reference_path):
        l = line.strip()
        if len(l) == 0:   continue
        if l[0] != '>':   ref_seq.append(l)
    ref_seq = ''.join(ref_seq)
    align = []
    for line in open(samfile):
        l = line.rstrip('\n')
        if len(l) == 0 or l[0] == '@':  continue
        parts = l.split('\t')
        flags = int(parts[1])
        if flags != 0 and flags != 16:  continue
        ID = parts[0].strip()
        ref_ind = int(parts[3])-1
        cigar = parts[5].strip()
        seq = parts[9].strip()
        edits = lowlevel_parse_cigar (cigar)
        this_sequence_name = ID
        if ref_ind > 0:  this_sequence = '-' * ref_ind # write gaps before alignment
        else:            this_sequence = ""  ## no gaps
        ind = 0; seq_len = ref_ind
        for e, e_len in edits:
            if e == 'M' or e == '=' or e == 'X': # (mis)match)
                this_sequence += seq[ind:ind+e_len]; ind += e_len; seq_len += e_len
            elif e == 'D':                       # deletion (gap in query)
                this_sequence += '-' * e_len; seq_len += e_len
            elif e == 'I':                       # insertion (gap in reference; ignore)
                ind += e_len
            elif e == 'S' or e == 'H':           # starting/ending segment of query not in reference (i.e., span of insertions; ignore)
                ind += e_len
        if seq_len < len(ref_seq):
            this_sequence += '-' * (len(ref_seq)-seq_len) # write gaps after alignment
        align.append (SeqRecord(Seq.Seq(this_sequence,Alphabet.IUPAC.ambiguous_dna), id=this_sequence_name, description=this_sequence_name))
     
    if infile is None:
        os.system(f"rm -f {ifl}")
    os.system(f"rm -f {samfile}")
    return align

def snpsites_from_alignment (sequences = None, strict = False, infile = None, outfile = None, prefix = "/tmp/"):
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me an alignment object or file")
    if prefix is None: prefix = "./"
    if strict is False: strict = ""
    else:               strict = "-c"
    if infile is None: ifl = prefix + "/wholegenome.aln"
    else:              ifl = infile
    if outfile is None: ofl = prefix + "/snpsites.aln"
    else:               ofl = outfile
    if sequences:       SeqIO.write(sequences, ifl, "fasta") # otherwise we assume sequence file exists 

    runstr = f"snp-sites {strict} {ifl} > {ofl}" # TODO: allow snp-sites -v for VCF
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    
    snps = AlignIO.read(ofl, "fasta")

    if infile is None:  os.system("rm -f " + ifl)
    #else:               os.system("bzip2 -f " + ifl) ## all fasta files shall be bzipped
    if outfile is None: os.system("rm -f " + ofl)

    return snps

def rapidnj_from_alignment (sequences = None, infile = None, outfile = None, prefix = "/tmp/", n_threads = 4):
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me an alignment object or file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = prefix + "/rapidnj.aln"
    else:              ifl = infile
    if outfile is None: ofl = prefix + "/rapidnj.tree"
    else:               ofl = outfile
    if sequences:       SeqIO.write(sequences, ifl, "fasta") # otherwise we assume sequence file exists 

    runstr = "rapidnj " + ifl + " -i fa -t d -n -c " + str(n_threads) + " -x " + ofl 
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    
    treestring = open(ofl).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = ete3.Tree (treestring)
    if infile is None:  os.system("rm -f " + ifl)
    #else:               os.system("bzip2 -f " + ifl) ## all fasta files shall be bzipped 
    if outfile is None: os.system("rm -f " + ofl)
    return tree

def pda_from_tree (tree, infile = None, outfile = None, prefix = "/tmp/", n_remain = 500):
    if (tree is None) and (infile is None):
        print ("ERROR: You must give me a tree object or file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = prefix + "bigtree.tre"
    else:              ifl = infile
    if outfile is None: ofl = prefix + "pda.tre"
    else:               ofl = outfile
    tree.write(format=1, outfile=ifl)
    n_remain = str(n_remain)
    runstr = f"iqtree -te {ifl} -k {n_remain}; tail -n 6 {ifl}.pda | head -1 > {ofl}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    
    treestring = open(ofl).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree_out = ete3.Tree (treestring)
    if infile is None:  os.system("rm -f " + ifl)
    if outfile is None: os.system("rm -f " + ofl)
    return tree_out

def improve_tree_from_align (tree, align, if_tre = None, of_tre = None, a_file = None, prefix = "/tmp/", n_threads = 4, params=None):
    if (tree is None) and (if_tre is None):
        print ("ERROR: You must give me a tree object or file")
    if (align is None) and (a_file is None):
        print ("ERROR: You must give me an alignment or file")
    if params is None: params = "-m HKY+G -me 0.05 -blmin 0.000005 -blmax 4"
    if prefix is None: prefix = "./"
    if if_tre is None: ifl = prefix + "in_iqtre.tre"
    else:              ifl = if_tre
    if of_tre is None: ofl = prefix + "out_iqtre.tre"
    else:              ofl = of_tre
    if a_file is None: aln = prefix + "seq.aln"
    else:              aln = a_file

    if align: SeqIO.write(align, aln, "fasta") ## else it should be present in infile
    if tree:  tree.write(format=1, outfile=ifl) ## to recycle the file make sure tree and align are None

    n_threads = str (n_threads)
    runstr = f"iqtree -g {ifl} -s {aln} {params} -ninit 1 -nt {n_threads} -redo; mv {aln}.treefile {ofl}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    

    treestring = open(ofl).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree_out = ete3.Tree (treestring)

#    os.system(f"rm -f {aln}.* ifl")
#    if a_file is None: os.system(f"rm -f {aln}")
#    if of_tre is None: os.system(f"rm -f {ofl}")
    return tree_out

iupac_dna = {''.join(sorted(v)):k for k,v in Seq.IUPAC.IUPACData.ambiguous_dna_values.items()}

def consensus_from_alignment (align): ## IUPAC ambiguity codes
    xaln = [SeqRecord(Seq.Seq(str(rec.seq).replace("-","N") ,Alphabet.IUPAC.ambiguous_dna), id=rec.id, description=rec.description) for rec in align]
    summary_align = AlignInfo.SummaryInfo(Align.MultipleSeqAlignment(align)) # must be an MSA, not a list
    pssm = summary_align.pos_specific_score_matrix(chars_to_ignore=["-"])
    consensus = [];
    # pssm example: {'-':3, 'A':0, 'T':4.0, 'G':0, 'C':2.0, 'N':1} per column, means 3 seqs have "-", 4 have "T"...
    for score in pssm: # we don't care about frequency, only presence
        # base can be "R", then iupac.dna_values[R] = [A,G]
        acgt_list = [x for base, count in score.items() for x in Seq.IUPAC.IUPACData.ambiguous_dna_values[base] if count > 0]
        consensus.append(iupac_dna[ ''.join(sorted(set(acgt_list))) ])
    return Seq.Seq(''.join(consensus),Alphabet.IUPAC.ambiguous_dna)
# profile could be [[x]*count for base,count...] with collection; or dict{A:0, C:0, ...} with dict[base]+=count

def distance_from_consensus (query, consensus):
    x = str(query.seq).replace("-","N")
    d1 = d2 = counter = 0
    for s1, s2 in zip(x, str(consensus.seq)):
        if s1 != "N" and s2 != "N":
            counter += 1
            if s1 != s2:
                d1 += 1 ## any difference
                l = len(set(Seq.IUPAC.IUPACData.ambiguous_dna_values[s1]).intersection(Seq.IUPAC.IUPACData.ambiguous_dna_values[s2]))
                if (l > 0): ## incompatible (e.g. W and A are compatible b/c W = A or T)
                    d2 +=1
    return [d1/counter, d2/counter] ## or return 3 values

def find_similar_scores (seqlist, threshold=1):
    scoremat = create_banded_score_matrix (seqlist)
    siblings = {}
    valid = np.repeat(1,scoremat.shape[0])
    for i in range(scoremat.shape[0]-1):
        for j in range(i+1, scoremat.shape[0]):
            if (valid[i]>0) and (valid[j]>0) and \
                ((abs(scoremat[i,i] - scoremat[i,j]) <= threshold) or \
                 (abs(scoremat[j,j] - scoremat[i,j]) <= threshold)):
                valid[j] = 0 # avoids d(a,b)=1 d(b,c)=1 but d(a,c)=2 being on same bucket since b becomes invalid
                if i in siblings:
                    siblings[i].append(j)
                else:
                    siblings[i] = [j]
    return siblings, valid

def score_to_distance_matrix_fraction (scoremat, mafft = False):
    distmat = np.zeros(scoremat.shape)
    offset = scoremat.min() - 1.
    scoremat -= offset
    if mafft is False: ## distance = fraction of indels
        for i,j in itertools.combinations(range(distmat.shape[0]),2):
            distmat[i,j] = distmat[j,i] = (scoremat[i,i] + scoremat[j,j]) / scoremat[i,j] - 2. # avoids division by zero
    else: # distance as Satoh in MAFFT 
        for i,j in itertools.combinations(range(distmat.shape[0]),2):
            distmat[i,j] = distmat[j,i] = 1. - scoremat[i,j]/min(scoremat[i,i],scoremat[j,j])
    return distmat

## rapidnj adds single quotes to names
## "China/Wuhan-Hu-1/2019" (NC_045512.2) identical to MN908947 acc to NCBI and is the longest, without Ns
