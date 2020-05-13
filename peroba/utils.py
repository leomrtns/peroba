#!/usr/bin/env python
from matplotlib import cm, colors # colormap
from Bio import Seq, SeqIO, Align, AlignIO, Phylo, Alphabet, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo, Applications
from Bio.Blast import NCBIXML
from Bio.Phylo import draw, TreeConstruction  #TreeConstruction.DistanceCalculator, TreeConstruction.DistanceTreeConstructor
from pastml.acr import acr
import ete3 # also used by pastml but may not be accessible
# https://bioinformatics.stackexchange.com/questions/4337/biopython-phylogenetic-tree-replace-branch-tip-labels-by-sequence-logos

import numpy as np, pandas as pd, seaborn as sns
from sklearn import manifold, metrics, cluster, neighbors, decomposition, preprocessing
import skbio, parasail, dendropy, datetime, time, codecs, joypy
import sys, gzip, bz2, re, glob, pickle, collections, subprocess, os, errno, random, itertools, pathlib

def print_redblack(textr="", textb=""):
    print ('\x1b[0;1;31;1m'+ str(textr) + '\x1b[0;1;30;1m'+ str(textb) + '\x1b[0m')

def read_fasta (filename, zip = "bz2", check_name = False, debug = False):
    unaligned = []
    if (zip == "bz2"):
        this_open = bz2.open
    elif (zip == "gz"):
        this_open = gzip.open
    else:
        this_open = open
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
            unaligned.append(record)
    if (debug):
        print_redblack ("first record length: " + str(len(unaligned[0].seq)))
        print (unaligned[0])

    return unaligned

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

def snpsites_from_alignment (sequences = None, infile = None, outfile = None, prefix = "/tmp/"):
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me an alignment object or file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = prefix + "/wholegenome.aln"
    else:              ifl = infile
    if outfile is None: ofl = prefix + "/snpsites.aln"
    else:               ofl = outfile
    if sequences:       SeqIO.write(sequences, ifl, "fasta") # otherwise we assume sequence file exists 

    runstr = f"snp-sites {ifl} > {ofl}" # TODO: allow snp-sites -v for VCF
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)    
    snps = AlignIO.read(ofl, "fasta")

    if infile is None:  os.system("rm -f " + ifl)
    #else:               os.system("bzip2 -f " + ifl) ## all fasta files shall be bzipped
    if outfile is None: os.system("rm -f " + ofl)

    return snps

iupac_dna = {''.join(sorted(v)):k for k,v in Seq.IUPAC.IUPACData.ambiguous_dna_values.items()}

## TODO: using reference MN908947.3, trim to [265, 29675]
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
    tree = ete3.Tree (ofl)
    if infile is None:  os.system("rm -f " + ifl)
    #else:               os.system("bzip2 -f " + ifl) ## all fasta files shall be bzipped 
    if outfile is None: os.system("rm -f " + ofl)
    return tree

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
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2', 'days_since_Dec19']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    df1 = pd.read_csv (str(filename), compression=compression, sep=sep) # sep='\t' for gisaid

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

    if index_name not in list(df1.columns): # regular CSV file from external source
        df1.set_index (str(primary_key), drop = False, inplace = True) # dont drop the column to be used as index
        df1.dropna (subset=[str(primary_key)], inplace = True); # now we have a column and an index with same name
        df1.rename_axis(str(index_name), inplace = True) # equiv. to df.index.name="peroba_seq_uid"
    else:
        df1.set_index (str(index_name), drop = True, inplace = True) # drop the column to avoid having both with same name

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
        exclude_na_rows = ['collection_date','sequence_name','lineage', 'adm2', 'days_since_Dec19']
    if exclude_columns == "default":
        exclude_columns = ["is_travel_history","is_surveillance","is_hcw", "is_community", "outer_postcode", "travel_history"]

    df['days_since_Dec19'] = df['collection_date'].map(lambda a: get_days_since_2019(a, impute = True))
    df = df.sort_values(by=['lineage_support', 'days_since_Dec19'], ascending=[False, True])
    df["collection_datetime"] = pd.to_datetime(df["collection_date"], infer_datetime_format=False, errors='coerce')

    # default values for missing rows (in case we don't want to remove those rows)
    if not exclude_na_rows or "uk_lineage" not in exclude_na_rows:
        df["uk_lineage"] = df["uk_lineage"].replace(np.nan, "x", regex=True) 
    if not exclude_na_rows or "adm2" not in exclude_na_rows:
        df['adm2'].fillna(df.country, inplace=True)

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

def get_ancestral_trait_subtrees (tre, csv,  tiplabel_in_csv = "sequence_name", elements = 1, 
        trait_column ="adm2", trait_value = "NORFOLK", n_threads = 4, method = "DOWNPASS"):
    '''
    Returns ancestral nodes predicted to have a given trait_value (e.g. "NORFOLK") for a given trait column (e.g. "adm2").
    Also returns nodes scrictly monophyletic regarding value.
    If using parsimony then it's better to store only nodes with a single state (or create a binary state?) otherwise
    we may end up chosing too close to root node...    MAP works well to find almost monophyletic nodes, but is quite slow
    max likelihood methods: pastml.ml.MPPA, pastml.ml.MAP, pastml.ml.JOINT,
    max parsimony methods: pastml.parsimony.ACCTRAN, pastml.parsimony.DELTRAN, pastml.parsimony.DOWNPASS
    '''
    if elements < 1: elements = 1 ## inferred state cardinality (how many values are allowed on matching internal node?)
    if method not in ["MPPA", "MAP", "JOINT", "ACCTRAN", "DELTRAN", "DOWNPASS"]:
        method = "DOWNPASS"
    if tiplabel_in_csv: # o.w. we assume input csv already has it
        csv_column = csv[[tiplabel_in_csv, trait_column]] # two columns now, tiplabel is removed below
        csv_column.set_index("sequence_name", drop = True, inplace = True) # acr needs index mapping ete3 leaves
    else:
        csv_column = csv[[trait_column]]  ## nodes will have e.g. n.adm2 (access through getattr(n.adm2)
    if isinstance (trait_value, list):
        trait_value = trait_value[0] ## This function can handle only one value; for grouping try the get_binary version
    ## Ancestral state reconstruction of given trait
    result = acr (tre, csv_column, prediction_method = method, force_joint=False, threads=n_threads) ## annotates tree nodes with states (e.g. tre2.adm2)

    ## Find all internal nodes where trait_value is possible state (b/c is seen at tips below)
    matches = filter(lambda n: not n.is_leaf() and trait_value in getattr(n,trait_column) and
            len(getattr(n,trait_column)) <= elements, tre.traverse("preorder"))
    # print ([x.__dict__ for x in matches]) dictionary of attributes; ete3 also has n.features[] with dict keys

    stored_leaves = set () # set of leaf names (created with get_cached_content)
    subtrees = [] # list of non-overlapping nodes
    node2leaves = tre.get_cached_content(store_attr="name") # set() of leaves below every node; store leaf name only
    for xnode in matches:
        if not bool (stored_leaves & node2leaves[xnode]): # both are sets; bool is just to be verbose
            stored_leaves.update (node2leaves[xnode]) # update() is append() for sets ;)
            subtrees.append(xnode)
    mono = tre.get_monophyletic (values=trait_value, target_attr = trait_column) # from ete3 
    return subtrees, mono, result

def get_binary_trait_subtrees (tre, csv,  tiplabel_in_csv = "sequence_name", elements = 1, 
          trait_column ="adm2", trait_value = "NORFOLK", n_threads = 4, method = "DOWNPASS"):
    '''
    Instead of reconstructing all states, we ask if ancestral state is value or not. Still we allow for `elements` > 1
    (2 in practice), so that we accept "yes and no" ancestral nodes. 
    You can group trait values into a list 
    '''
    if elements < 1: elements = 1 ## inferred state cardinality (1 makes much more sense, but you can set "2" as well)
    if method not in ["MPPA", "MAP", "JOINT", "ACCTRAN", "DELTRAN", "DOWNPASS"]:
        method = "DOWNPASS"
    if tiplabel_in_csv: # o.w. we assume input csv already has it
        csv_column = csv[[tiplabel_in_csv, trait_column]] # two columns now, tiplabel is removed below
        csv_column.set_index("sequence_name", drop = True, inplace = True) # acr needs index mapping ete3 leaves
    else:
        csv_column = csv[[trait_column]]  ## nodes will have e.g. n.adm2 (access through getattr(n.adm2)
    if isinstance (trait_value, str): # assume it's a list, below
        trait_value = [trait_value]

    # transform variable into binary
    new_trait = str(trait_column) + "_is_" + "_or_".join(trait_value)
    csv_column[new_trait] = csv_column[trait_column].map(lambda a: "yes" if a in trait_value else "no")

    csv_column.drop(labels = [trait_column], axis=1, inplace = True)
    ## Ancestral state reconstruction of given trait
    result = acr (tre, csv_column, prediction_method = method, force_joint=False, threads=n_threads) ## annotates tree nodes with states (e.g. tre2.adm2)
    ## Find all internal nodes where trait_value is possible state (b/c is seen at tips below)
    matches = filter(lambda n: not n.is_leaf() and "yes" in getattr(n,new_trait) and # traits are sets (not lists)
            len(getattr(n,new_trait)) <= elements, tre.traverse("preorder"))

    stored_leaves = set () # set of leaf names (created with get_cached_content)
    subtrees = [] # list of non-overlapping nodes
    node2leaves = tre.get_cached_content(store_attr="name") # set() of leaves below every node; store leaf name only
    for xnode in matches:
      if not bool (stored_leaves & node2leaves[xnode]): # both are sets; bool is just to be verbose
          stored_leaves.update (node2leaves[xnode]) # update() is append() for sets ;)
          subtrees.append(xnode)
    mono = tre.get_monophyletic (values = "yes", target_attr = new_trait) # from ete3 
    return subtrees, mono, result

def colormap_from_dataframe (df, column_list, column_names, cmap_list = None):
    ''' returns a dictionary of lists, where key is the tip label (should be index of csv)
    column_list and column_names must have same length, mapping to names on CSV and on plot
        the default colormap are qualitative so I guess they fail if #elements>#colours...
    '''
    if cmap_list is None:
        cmap_list = ["Accent", "Dark2", "jet", "hsv", "viridis", "plasma", "cividis", "rainbow", "nipy_spectral"]
    if isinstance (cmap_list, str): # assume it's a list, below
        cmap_list = [cmap_list]
    if len(column_list) != len(column_names):
        print ("ops, list of columns (from table) must have an equivalent name (that describes it)")
        return
    ## iterate over columns (phenotypes), to generate colormaps
    d_col = {}
    for i, (c, cname) in enumerate(zip(column_list, column_names)):
        uniq = df[c].unique()
        cmap_elem = cmap_list[i%len(cmap_list)] # cycle over given colormaps
        if cmap_elem in ['Pastel1', 'Pastel2', 'Paired', 'Accent','Dark2', 'Set1', 
                'Set2', 'Set3','tab10', 'tab20', 'tab20b', 'tab20c']: # only 8~20 colours
            customCMap = colors.LinearSegmentedColormap.from_list("custom", [(x,cm.get_cmap(cmap_elem)(x)) for x in np.linspace(0, 1, 8)])
        else:
            customCMap = cm.get_cmap(cmap_elem)
        colorlist = customCMap(np.linspace(0, 1, len(uniq)))
        d_col[cname] = {name:colors.to_hex(col) for name,col in zip(uniq, colorlist)} ## each value is another dict from csv elements to colors
    
    col_missing = colors.to_hex([1,1,1,0]) ## transparent, default to white
    ## iterate over rows (sequences)
    d_seq = {}
    d_seq_lab = {}
    for seq in df.index: # frown upon by pandas wizards
        d_seq[seq] = []
        d_seq_lab[seq] = []
        for c, cname in zip (column_list, column_names):
            if df.loc[seq,c] in d_col[cname]: # valid value for column
                d_seq[seq].append(d_col[cname][ df.loc[seq,c] ])
                d_seq_lab[seq].append(str(df.loc[seq,c]))
            else:
                d_seq[seq].append(col_missing)
                d_seq_lab[seq].append(str(" "))
    return [d_seq, d_seq_lab, d_col, column_names] 

def return_treestyle_with_columns (cmapvector):
    '''
    Need column names again to print header in order 
    '''
    [d_seq, d_seq_lab, d_col, column_names] = cmapvector
    label_font_size = 6

    # default node
    ns1 = ete3.NodeStyle()
    ns1["size"] = 1 ## small dot 
    ns1["shape"] = "square" ## small dot 
    ns1["fgcolor"] = "101010" ## small dot 
    ns1["hz_line_type"]  = ns1["vt_line_type"]  = 0 # 0=solid, 1=dashed, 2=dotted
    ns1["hz_line_color"] = ns1["vt_line_color"] = "#0c0c0c"
    ns2 = ete3.NodeStyle()
    ns2["size"] = 0 ## no dot
    ns2["shape"] = "circle" ## small dot 
    ns2["hz_line_type"]  = ns1["vt_line_type"]  = 0 # 0=solid, 1=dashed, 2=dotted
    ns2["hz_line_color"] = ns1["vt_line_color"] = "#0c0c0c"

    ## prepare table and other node information
    def tree_profile_layout (node):
        if node.is_leaf(): # the aligned leaf is "column 0", thus traits go to column+1
            node.set_style(ns1) ## may be postponed to when we have ancestral states
            ete3.add_face_to_node(ete3.AttrFace("name", fsize=label_font_size, text_suffix="   "), node, 0, position="aligned")
            for column, (rgb_val, lab) in enumerate(zip(d_seq[node.name], d_seq_lab[node.name])): ## colour of csv.loc[node.name, "adm2"]
                label = {"text": lab[:10], "color":"Black", "fontsize": label_font_size-1}
                ete3.add_face_to_node (ete3.RectFace (50, 10, fgcolor=rgb_val, bgcolor=rgb_val, label = label), node, column+1, position="aligned")
        else:
            node.set_style(ns2) 
            #node.img_style['hz_line_color'] = node.img_style['vt_line_color'] =

    ts = ete3.TreeStyle()
    ts.draw_guiding_lines = True # dotted line between tip and name
    ts.guiding_lines_color = "#bdb76d"
    ts.guiding_lines_type = 2 #  0=solid, 1=dashed, 2=dotted
    ts.layout_fn = tree_profile_layout
    ts.branch_vertical_margin = 0
    ts.min_leaf_separation = 1 # Min separation, in pixels, between two adjacent branches
    ts.scale = 2000000 # 2e6 pixels per branch length unit (i.e. brlen=1 should be how many pixels?)
    ts.show_scale = False
    show_branch_length = True
    ts.show_leaf_name = False # we handle this in the layout function

    ## STILL dont know how to do it
    #ts.legend.add_face(CircleFace(10, "red"), column=0)
    #ts.legend.add_face(TextFace("0.5 support"), column=1)
    #ts.legend_position = 3 #  TopLeft corner if 1, TopRight if 2, BottomLeft if 3, BottomRight if 4
    
    for col, label in enumerate([""] + column_names): # the first are tip labels
        labelFace = ete3.TextFace(label, fsize=10, fgcolor="DimGray") # fsize controls interval betweel columns
        labelFace.rotation = 290
        labelFace.vt_align = 1  # 0 top, 1 center, 2 bottom
        labelFace.hz_align = 2  # 0 left, 1 center, 2 right 
        ts.aligned_header.add_face(labelFace, col)

    return ts

## rapidnj adds single quotes to names
## "China/Wuhan-Hu-1/2019" (NC_045512.2) identical to MN908947 acc to NCBI and is the longest, without Ns
