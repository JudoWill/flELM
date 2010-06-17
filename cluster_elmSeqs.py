"""Host and flu share few ELM sequences,
   so I'm making clusters of ELM sequences
   to be able to make comparisons between
   host and flu.

   I'll start with 5 clusters where applicable.
"""
import Bio.Cluster
import Levenshtein
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

def mk_dis_mat(strings):
    """Make distance matrix"""

    dis_mat = []
    for s1 in strings:
        for s2 in strings:
            if s1 == s2:
                break
            else:
                dis = Levenshtein.distance(s1, s2)
                dis_mat.append(dis)
    return numpy.array(dis_mat)

def check_clusters(cluster2seq, host_counts, flu_counts, elm):
    """make sure flu sequences are not in clusters alone"""

    host_seqs = {}
    for host in host_counts:
        for elmSeq in host_counts[host]:
            helm,hseq = elmSeq.split(':')
            if helm == elm:
                host_seqs[hseq] = True
    flu_seqs = {}
    for flu in flu_counts:
        for elmSeq in flu_counts[flu]:
            helm,hseq = elmSeq.split(':')
            if helm == elm:
                flu_seqs[hseq] = True

    for cluster in cluster2seq:
        found_host = False
        found_flu = False
        for seq in cluster2seq[cluster]:

            if seq in host_seqs:
                found_host = True
            if seq in flu_seqs:
                found_flu = True
        if found_flu and not found_host:
            return False
    return True

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
                elm, sequence = elmSeq.split(':')
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

flu_counts['human'] = count_flu(pre_flu_counts['human'], all_elmSeqs)
flu_counts['chicken'] = count_flu(pre_flu_counts['chicken'], all_elmSeqs)

for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            all_elmSeqs[elmSeq] = True
            host_counts[host][elmSeq] += int(count)

elm2seq = defaultdict(list)
for elmSeq in all_elmSeqs:
    elm, seq = elmSeq.split(':')
    elm2seq[elm].append(seq)

num_clusters = 5
elm2cluster2seq = {}
for elm in elm2seq:
    elm2cluster2seq[elm] = defaultdict(dict)
    strings = elm2seq[elm]
    if len(strings) > num_clusters:
        ans, error, nfound = Bio.Cluster.kmedoids(mk_dis_mat(strings),
                                                  nclusters=num_clusters,
                                                  npass=100)
        for seq, cluster in zip(strings, ans):
            elm2cluster2seq[elm][cluster][seq] = True
            print('%s\t%s\t%d' %
                  (elm, seq, cluster))
    else:
        ans, error, nfound = Bio.Cluster.kmedoids(mk_dis_mat(strings),
                                                  nclusters=2,
                                                  npass=100)
        for seq, cluster in zip(strings, ans):
            elm2cluster2seq[elm][cluster][seq] = True
            print('%s\t%s\t%d' %
                  (elm, seq, cluster))

for elm in elm2cluster2seq:
    print elm, check_clusters(elm2cluster2seq[elm], host_counts,
                              pre_flu_counts, elm)
# for s, cluster in zip(strings, ans):
#     print s, cluster



