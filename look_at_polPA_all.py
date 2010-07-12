"""Check host/flu hypothesis for ELM/seq/length trios."""
from collections import defaultdict
import os, utils, utils_graph, sys

good_elms = utils_graph.getNodes('working/Jul7/good_phylogeny_elms')
use_protein = sys.argv[1]
print 'RUN', use_protein

def print_it(title, uniq, chicken_host_freqs, human_host_freqs):
    """Print counts of ELM/seqs that fit hypothesis"""

    if len(uniq) == 0: return

    print title
    sums = defaultdict(utils.init_zero)
    totals = defaultdict(utils.init_zero)
    seen_seqs = {}
    
    for elmseq in uniq:
        protein, elm, seq, length = elmseq.split(':')
        key = elm + ':' + seq + ':' + length
        chicken_host_freq = float(0)
        human_host_freq = float(0)
        if key in chicken_host_freqs:
            chicken_host_freq = chicken_host_freqs[key]
        if key in human_host_freqs:
            human_host_freq = human_host_freqs[key]
        if not (chicken_host_freq == float(0) and human_host_freq == float(0)) and elm in good_elms and chicken_host_freq != float(1) and human_host_freq != float(1):
                diff = chicken_host_freq-human_host_freq
                if diff > float(0):
                    sums[elm] += 1
                totals[elm] += 1
            #seen_seqs[seq] = True

    full_sum = 0
    full_total = 0
    for elm in totals:
        full_sum += sums[elm]
        full_total += totals[elm]
    print full_sum, full_total-full_sum, float(100)*float(full_sum)/float(full_total)

def get_host_freqs(afile):
    freqs = {}
    len_counts = defaultdict(utils.init_zero)
    with open(afile) as f:
        for line in f:
            elm_len, seq, count, freq = line.strip().split('\t')
            elm, length = elm_len.split(':')
            freqs[elm + ':' + seq + ':' + length] = float(freq)
    return freqs

def get_freqs(afile, seq_percents):
    """Look at flu seq % coverage"""

    seen = {}
    with open(afile) as f:
        for line in f:
            protein, elmseq, freq = line.strip().split('\t')
            key = ':'.join([protein, elmseq.split(':')[1],
                            str(len( elmseq.split(':')[1]))])
            if key not in seen:
                seq_percents[key].append(float(freq))
                seen[key] = True

def get_annotations(afile):
    """Get ELM seqs that are conserved"""

    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            protein, elm = line.strip().split('\t')
            if protein == use_protein:
                seq = elm + ':' + str(len(elm.split(':')[1]))
                d[protein][seq] = True
    return d

def get_uniq(uniq, ls1, ls2):
    """Of the conserved ELM seqs, which ones are unique?"""

    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[protein + ':' + elm] = True

uniq_human = {}
with open('working/Jul7/mammal_uniq') as f:
    for line in f:
        protein, elm = line.strip().split('\t')
        if protein == use_protein:
            length = str(len(elm.split(':')[1]))
            uniq_human[':'.join((protein, elm, length))] = True

uniq_bird = {}
with open('working/Jul7/bird_uniq') as f:
    for line in f:
       protein, elm = line.strip().split('\t')
       if protein == use_protein:
           length = str(len(elm.split(':')[1]))
           uniq_bird[':'.join((protein, elm, length))] = True

control = {}
with open('working/Jul7/control') as f:
    for line in f:
       protein, elm = line.strip().split('\t')
       if protein == use_protein:
           length = str(len(elm.split(':')[1]))
           control[':'.join((protein, elm, length))] = True

# human_elmseqs_file = 'working/Jul7/mammal_elms'
# bird_elmseqs_file = 'working/Jul7/bird_elms'

# human_elmseqs = get_annotations(human_elmseqs_file)
# bird_elmseqs = get_annotations(bird_elmseqs_file)

# uniq_human = {}
# uniq_bird = {}
# get_uniq(uniq_human, human_elmseqs, bird_elmseqs)
# get_uniq(uniq_bird, bird_elmseqs, human_elmseqs)

# notUniq_human = {}
# for protein in human_elmseqs:
#     for elm in human_elmseqs[protein]:
#         key = protein + ':' + elm
#         if key not in uniq_human:
#             notUniq_human[key] = True

# notUniq_bird = {}
# for protein in bird_elmseqs:
#     for elm in bird_elmseqs[protein]:
#         key = protein + ':' + elm
#         if key not in uniq_bird:
#             notUniq_bird[key] = True

print  len(uniq_bird), len(uniq_human)

dir = 'working/Jul7'
years = range(2000,2011,1)

human_host_file = 'working/Jul7/elmdict_Sus_scrofa.RWlenInit'
chicken_host_file = 'working/Jul7/elmdict_Gallus_gallus.RWlenInit'
human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)

print_it('BIRD', uniq_bird, chicken_host_freqs,
         human_host_freqs)
print_it('NOT UNIQ BIRD', control, chicken_host_freqs,
         human_host_freqs)
print '\n'
print_it('HUMAN', uniq_human, human_host_freqs,
         chicken_host_freqs)
print_it('NOT UNIQ HUMAN', control, human_host_freqs,
         chicken_host_freqs)
print '\n'
