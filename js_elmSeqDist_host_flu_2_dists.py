"""Make distributions of chicken/human host/flu ELMseq counts"""
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

def count_0s(ls):
    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

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
    # flu_counts[flu] = defaultdict(utils.init_zero)
    # with open('results/' + flu + '.H5N1.elms') as f:
    #     for line in f:
    #         (protein, st, stp,
    #          elm, seq, junk) = line.strip().split('\t')
    #         elmSeq = elm + ':' + seq
    #         all_elmSeqs[elmSeq] = True
    #         flu_counts[flu][elmSeq] += 1

for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            all_elmSeqs[elmSeq] = True
            host_counts[host][elmSeq] += int(count)

flu_vecs = mk_count_vecs(flu_counts, all_elmSeqs)
flu_dists = mk_count_dists(flu_vecs)
host_vecs = mk_count_vecs(host_counts, all_elmSeqs)
host_dists = mk_count_dists(host_vecs)

def print_it(name, vec):
    print name, float(count_0s(vec))/float(len(vec))

print_it('chicken_flu', flu_vecs['chicken'])
print_it('human flu', flu_vecs['human'])
print_it('H sapiens', host_vecs['H_sapiens'])
print_it('Gallus gallus', host_vecs['Gallus_gallus'])

flu_dists['chicken'].sort()
print flu_dists['chicken'][-3], flu_dists['chicken'][-2], flu_dists['chicken'][-1]
