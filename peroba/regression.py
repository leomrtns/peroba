#!/usr/bin/env python

from peroba.utils import *
import logging, ete3, argparse
import numpy as np
from Bio import Seq, SeqIO, Align, AlignIO, Phylo, Alphabet, pairwise2
import time, codecs, random
from sklearn.neighbors import NearestNeighbors, BallTree ## KDTree does NOT accept hamming
from sklearn.cluster import OPTICS
import xxhash

logger = logging.getLogger(__name__) 
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

def list_duplicates (seq_dict, blocks = 4, leaf_size = 500, radius=0.00001):
    aln = [x for x in seq_dict.values()]
    genome_size = len(aln[0].seq) ## only works for aligned sequences
    block_size = int(genome_size / blocks)
    if (block_size < 1): block_size = 1
    hashes = [[xxhash.xxh32(str(aln[j].seq[i:i+block_size])).intdigest() for i in range(0,genome_size,block_size)] for j in range(len(aln))]
    btre = BallTree(np.array(hashes), leaf_size=leaf_size, metric='hamming') # create a neighbours tree (KDTree doesnt accept hamming)
    idx = btre.query_radius(hashes, r=radius, return_distance=False) ## return_distace=False makes it faster; r=0.0001 is essentially zero
    #clusters = [[aln[j].id for j in x] for x in idx if len(x)>1] # only those with duplicates (len>1) are returned
    clusters = [[aln[j].id for j in x] for x in idx] # return all 
    del aln, hashes, btre, idx
    return clusters

def list_r_neighbours (g_seq, l_seq, blocks = 1000, leaf_size = 500, dist_blocks = 1):
    g_aln = [x for x in g_seq.values()]
    l_aln = [x for x in l_seq.values()]
    genome_size = len(g_aln[0].seq) ## only works for aligned sequences
    block_size = int(genome_size / blocks)
    if (block_size < 1): block_size = 1
    if dist_blocks < 1 : dist_blocks = 1
    dist_blocks += 0.1 # to ensure those within dist are returned (i.e. "<=" and not strictly "<")
    radius = dist_blocks / float(blocks)
    logger.info("Creating a hashed genome with blocks of %s bases",str(block_size))
    g_hash = [[xxhash.xxh32(str(g_aln[j].seq[i:i+block_size])).intdigest() for i in range(0,genome_size,block_size)] for j in range(len(g_aln))]
    btre = BallTree(np.array(g_hash), leaf_size=leaf_size, metric='hamming') # create a neighbours tree of global sequences

    # hashed query sequences
    l_hash = [[xxhash.xxh32(str(l_aln[j].seq[i:i+block_size])).intdigest() for i in range(0,genome_size,block_size)] for j in range(len(l_aln))]
    clusters = []
    n_iter = 0
    while len(clusters) < 1 and n_iter < 10:
        logger.info("Iter {}: trying to find neighbours with distance smaller than {:.6f}".format(n_iter, radius))
        idx = btre.query_radius(l_hash, r=radius, return_distance=False) # gives global neighbours to each local sequence; return_distance is expensive
        clusters = list(set([g_aln[j].id for x in idx for j in x])) # one-dimentional list of all global neighbours
        radius = radius * 1.4; 
        n_iter += 1;
    if (len(clusters) < 1):
        logger.warning ("Could not find neighours close enough;")
    del g_aln, g_hash, l_aln, l_hash, btre, idx
    return clusters

def list_n_neighbours (g_seq, l_seq, blocks = 1000, leaf_size = 200, nn = 10):
    g_aln = [x for x in g_seq.values()]
    l_aln = [x for x in l_seq.values()]
    genome_size = len(g_aln[0].seq) ## only works for aligned sequences
    block_size = int(genome_size / blocks)
    if (block_size < 1): block_size = 1
    logger.info("Creating a hashed genome with blocks of %s bases",str(block_size))
    g_hash = [[xxhash.xxh32(str(g_aln[j].seq[i:i+block_size])).intdigest() for i in range(0,genome_size,block_size)] for j in range(len(g_aln))]
    btre = BallTree(np.array(g_hash), leaf_size=leaf_size, metric='hamming') # create a neighbours tree of global sequences

    logger.info("And finding %s closest neighbours",str(nn)) 
    if (nn<2): logger.warning("Closest neighbour will be itself, if already on BallTree; useful only on two independent sets")

    l_hash = [[xxhash.xxh32(str(l_aln[j].seq[i:i+block_size])).intdigest() for i in range(0,genome_size,block_size)] for j in range(len(l_aln))]
    dist, idx = btre.query (l_hash, k=nn, return_distance=True) # return_distance is free 
    clusters = list(set([g_aln[j].id for x in idx for j in x])) # one-dimentional list of all global neighbours
    del g_aln, g_hash, l_aln, l_hash, btre, idx
    return clusters

def list_paf_neighbours (g_seq, l_seq, n_segments = 1, n_best = 20, n_threads = 2): # nthreads increase peak memory 
    if n_segments > 10: n_segments = 10
    if n_segments < 1:  n_segments = 1
    g_aln = [x for x in g_seq.values()] # calling function is faster if working with dictionary than list (alignment)
    l_aln = [x for x in l_seq.values()]
    genome_size = len(g_aln[0].seq) ## only works for aligned sequences
    block_size = int(genome_size / n_segments)
    #(str(g_aln[j].seq[i:i+block_size])) 
    #for j in range(len(g_aln))]
    clusters = []
    for i in range(0, genome_size, block_size):
        loc  = [SeqRecord(Seq.Seq(str(rec.seq[i:i+block_size])),id=rec.id,description=rec.description) for rec in l_aln]
        glob = [SeqRecord(Seq.Seq(str(rec.seq[i:i+block_size])),id=rec.id,description=rec.description) for rec in g_aln]
        clusters += minimap2_find_neighbours (query_aln = loc, target_aln = glob, batch_size = 64, n_best=n_best, n_threads = n_threads) # returns list of  seqnames
    return list(set(clusters))
