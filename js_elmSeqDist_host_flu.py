"""Use Jensen-Shannon divergence to make a dendrogram for eukaryotic hosts"""
import itertools, sys, os, utils, random, global_settings, numpy, utils_graph
from collections import defaultdict

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

hosts = ('Sus_scrofa', 'Gallus_gallus')
flus = ('swine', 'chicken', 'human', 'duck', 'equine')

# count elm:seq occurence
flu_counts = {}
host_counts = {}
all_elmSeqs = {}
for flu in flus:
    flu_counts[flu] = defaultdict(utils.init_zero)
    with open('results/' + flu + '.H9N2.elms') as f:
        for line in f:
            (protein, st, stp, elm, seq, junk) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm == 'TRG_ENDOCYTIC_2':
                all_elmSeqs[elmSeq] = True
                flu_counts[flu][elmSeq] += 1
print flu_counts
host_all_Seqs = {}
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('working/Jun29/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm == 'TRG_ENDOCYTIC_2':
                #all_elmSeqs[elmSeq] = True
                host_all_Seqs[elmSeq] = True
                host_counts[host][elmSeq] += int(count)
print len(utils_graph.intersectLists([host_all_Seqs,
                                      all_elmSeqs]))
flu_vecs = mk_count_vecs(flu_counts, all_elmSeqs)
flu_dists = mk_count_dists(flu_vecs)
host_vecs = mk_count_vecs(host_counts, all_elmSeqs)
host_dists = mk_count_dists(host_vecs)

js_distances = defaultdict(dict)
for host in hosts:
    for flu in flus:
        js_dis = utils.jensen_shannon_dists(host_dists[host],
                                            flu_dists[flu])
        js_distances[host][flu] = js_dis
        print host, flu, js_dis
