"""Host and flu share few ELM sequences,
   so I'm making clusters of ELM sequences
   to be able to make comparisons between
   host and flu.

   I want to see how the closest distances
   are distributed.
"""
import Bio.Cluster
import Levenshtein, sys
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

def get_closest_distances(flu_strings, host_strings):
    """For each flu sequence, find the closet distance to a human string"""

    min_distances = []

    for flu_s in flu_strings:
        dis = Levenshtein.distance(flu_s, 
                                   host_strings[0])
        for host_s in host_strings[1:]:
            d = Levenshtein.distance(flu_s, 
                                     host_s)
            dis = min([dis, d])
        print dis

def print_closest_distances(elm, flu_strings, host_strings):
    """For each flu sequence, find the distance to a human string.
       If it is 0,1,2,3 print it."""

    for flu_s in flu_strings:
        for host_s in host_strings:
            dis = Levenshtein.distance(flu_s, 
                                       host_s)
            if dis < 4:
                print elm + '\t' + flu_s + '\t' + host_s + '\t' + str(dis)

def count_flu(protein2counts, elm2seq):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                elm, sequence = elmSeq.split(':')
                counts[elmSeq] += protein2counts[protein][seq][elmSeq]
                elm2seq[elm][sequence] = True
    return counts

def get_flu_counts(afile, proteins):
    """Make protein_name -> seq_name -> elm_seq_counts"""

    counts = {}
    with open(afile) as f:
        for line in f:
            (protein, st, stp,
             elm, seq, junk) = line.strip().split('\t')
            name = protein.split('.')[-1]
            elmSeq = elm + ':' + seq
            if name in proteins:
                if name not in counts:
                    counts[name] = {}
                if protein not in counts[name]:
                    counts[name][protein] = {}
                if elmSeq not in counts:
                    counts[name][protein][elmSeq] = 0
                counts[name][protein][elmSeq] += 1
    return counts

def mk_vec(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for this host's counts"""
    
    vec = []
    for elmseq in all_elmSeqs:
        if elmseq in counts:
            vec.append(counts[elmseq])
        else:
            vec.append(float(0))
    return vec

def mk_count_vecs(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for all hosts"""

    vecs = {}
    for host in counts:
        vecs[host] = mk_vec(counts[host],
                            all_elmSeqs)
    return vecs

def mk_count_dists(vecs):
    """change count vectors into distributions"""

    dists = {}
    for host in vecs:
        dists[host] = utils.getDistFromCount(vecs[host])
    return dists

hosts = global_settings.TEST_GENOMES
flus = ('human', 'chicken')

# count elm:seq occurence
flu_elm2seq = defaultdict(dict)
pre_flu_counts = {}
flu_counts = {}
for flu in flus:
    pre_flu_counts[flu] = get_flu_counts('results/' + flu + '.H5N1.elms', 
                                         global_settings.FLU_PROTEINS)

flu_counts['human'] = count_flu(pre_flu_counts['human'], flu_elm2seq)
flu_counts['chicken'] = count_flu(pre_flu_counts['chicken'], flu_elm2seq)

host_elm2seq = defaultdict(dict)
for host in hosts:
    with open('working/runs/Jun25/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            host_elm2seq[elm][seq] = True

for elm in flu_elm2seq:
    print_closest_distances(elm, flu_elm2seq[elm], host_elm2seq[elm])



