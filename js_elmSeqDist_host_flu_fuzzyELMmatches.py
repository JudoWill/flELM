"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

# this comes from my scratch experiments
distance_file = '../../scratch/distances'

# map each host elmSeq to a flu elmSeq
host_flu_elmSeq_mapping = {}
host_flu_elmSeq_best_dis = {}
mapped = {}
with open(distance_file) as f:
    for line in f:
        (elm, flu_seq, host_seq, distance) = line.strip().split('\t')
        dis = int(distance)
        host_elmSeq = elm + ':' + host_seq
        flu_elmSeq = elm + ':' + flu_seq
        
        if host_elmSeq not in host_flu_elmSeq_best_dis:
            host_flu_elmSeq_best_dis[host_elmSeq] = {}
            host_flu_elmSeq_mapping[host_elmSeq] = {}
        if flu_elmSeq not in host_flu_elmSeq_best_dis[host_elmSeq]:
            host_flu_elmSeq_best_dis[host_elmSeq][flu_elmSeq] = dis
            host_flu_elmSeq_mapping[host_elmSeq][flu_elmSeq] = True
        if dis < host_flu_elmSeq_best_dis[host_elmSeq][flu_elmSeq]:
            # replace
            host_flu_elmSeq_best_dis[host_elmSeq][flu_elmSeq] = dis
            host_flu_elmSeq_mapping[host_elmSeq] = {}
            host_flu_elmSeq_mapping[host_elmSeq][flu_elmSeq] = True
            mapped[flu_elmSeq] = True
        elif dis == host_flu_elmSeq_best_dis[host_elmSeq][flu_elmSeq]:
            # append
            host_flu_elmSeq_mapping[host_elmSeq][flu_elmSeq] = True
            mapped[flu_elmSeq] = True

def count_flu(protein2counts, all_elmSeqs):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                counts[elmSeq] += protein2counts[protein][seq][elmSeq]
                all_elmSeqs[elmSeq] = True
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

hosts = ('H_sapiens', 'Gallus_gallus')
flus = ('human', 'chicken')
proteins = ('hemagglutinin', 'neuraminidase', 'nucleocapsid protein',
            'matrix protein 1', 'nonstructural protein 1', 'matrix protein 2',
            'nonstructural protein 2', 'polymerase PA', 'polymerase PB2',
            'polymerase PB1', 'PB1-F2 protein')

# count elm:seq occurence
flu_counts = {}
pre_flu_counts = {}
host_counts = {}
all_elmSeqs = {}
for flu in flus:
    pre_flu_counts[flu] = get_flu_counts('results/' + flu + '.H5N1.elms', 
                                         proteins)

# sample protein sequenes from chicken
new_chicken_counts = {}
new_chicken_proteins = {}
for protein in proteins:
    c = len(pre_flu_counts['human'][protein].keys())
    new_chicken_proteins[protein] = random.sample(pre_flu_counts['chicken'][protein], c)

for protein in new_chicken_proteins:
    new_chicken_counts[protein] = {}
    for seq in new_chicken_proteins[protein]:
        new_chicken_counts[protein][seq] = pre_flu_counts['chicken'][protein][seq]

flu_counts['human'] = count_flu(pre_flu_counts['human'], all_elmSeqs)
flu_counts['chicken'] = count_flu(new_chicken_counts, all_elmSeqs)

for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elmSeq in host_flu_elmSeq_mapping:
                for a_flu_elmSeq in host_flu_elmSeq_mapping[elmSeq]:
                    host_counts[host][a_flu_elmSeq] += int(count)

flu_vecs = mk_count_vecs(flu_counts, mapped)
flu_dists = mk_count_dists(flu_vecs)
host_vecs = mk_count_vecs(host_counts, mapped)
host_dists = mk_count_dists(host_vecs)

js_distances = defaultdict(dict)
for host in hosts:
    for flu in flus:
        js_dis = utils.jensen_shannon_dists(host_dists[host],
                                            flu_dists[flu])
        js_distances[host][flu] = js_dis
        print host, flu, js_dis
