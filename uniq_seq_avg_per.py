"""Find the average # of proteins
   each uniq seq is on in bird
   and human strains. Look for
   a correlation between host and
   flu."""
from collections import defaultdict
import utils, os, sys

use_host = sys.argv[1]

def get_annotations(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            protein, elm = line.strip().split('\t')
            d[protein][elm] = True
    return d

def get_important(ls1, ls2):
    """Get annotations important on flu"""

    elms = {}
    for ls in (ls1, ls2):
        for protein in ls:
            for elm in ls[protein]:
                elms[protein + ':' + elm] = True
    return elms

def get_uniq(uniq, ls1, ls2):
    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[protein + ':' + elm] = True

def get_freqs(afile, seq_percents):

    with open(afile) as f:
        for line in f:
            protein, elmseq, freq = line.strip().split('\t')
            key = ':'.join([protein, elmseq])
            seq_percents[key].append(float(freq))

def get_host_freqs(afile):
    freqs = {}
    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            freqs[elm + ':' + seq] = float(freq)
    return freqs

dir = 'working/Jul1_year'
years = range(2000,2011,1)

# bird
birds = ('chicken',)
strains = ('H9N2', 'H5N1')
seq_percents_bird = defaultdict(list)
total_bird = 0
for bird in birds:
    for strain in strains:
        for year in years:
            file = os.path.join(dir, '.'.join((bird, strain, str(year),
                                               'elmseqs.conservation')))
            if os.path.exists(file):
                total_bird += 1
                get_freqs(file, seq_percents_bird)

# human
hosts = ('human',)
'H1N1', 'H3N2', 'H3N8'
strains = ('H5N1', )
seq_percents_human = defaultdict(list)
total_human = 0
for host in hosts:
    for strain in strains:
        for year in years:
            file = os.path.join(dir, '.'.join((host, strain, str(year),
                                               'elmseqs.conservation')))
            if os.path.exists(file):
                total_human += 1
                get_freqs(file, seq_percents_human)

human_elmseqs_file = 'working/Jul1_year/mammal_elmseqs'
bird_elmseqs_file = 'working/Jul1_year/bird_elmseqs'

human_elmseqs = get_annotations(human_elmseqs_file)
bird_elmseqs = get_annotations(bird_elmseqs_file)

uniq = {}
get_uniq(uniq, human_elmseqs, bird_elmseqs)
get_uniq(uniq, bird_elmseqs, human_elmseqs)
impt_elms = get_important(human_elmseqs, bird_elmseqs)
# uniq = {}
if use_host == 'mammal':
    #get_uniq(uniq, human_elmseqs, bird_elmseqs)
    sp = seq_percents_human
elif use_host == 'bird':
    #get_uniq(uniq, bird_elmseqs, human_elmseqs)
    sp = seq_percents_bird

# for elmseq in seq_percents_bird:
#     if elmseq in uniq:
#         print elmseq + '\t' + str(sum(seq_percents_bird[elmseq])/float(total_bird)) + '\t' + str(sum(seq_percents_human[elmseq])/float(total_human))

human_host_file = 'working/Jun29/elmdict_H_sapiens.simple'
chicken_host_file = 'working/Jun29/elmdict_Gallus_gallus.simple'
human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)

for elmseq in uniq:
    dels = elmseq.split(':')
    new_seq = utils.mk_sub(dels[2])
    key = dels[1] + ':' + new_seq
    human_host = 0
    chicken_host = 0
    if key in chicken_host_freqs:
        chicken_host = chicken_host_freqs[key]
    if key in human_host_freqs:
        human_host = human_host_freqs[key]
    bird_flu_frac = sum(seq_percents_bird[elmseq])/float(total_bird)
    human_flu_frac = sum(seq_percents_human[elmseq])/float(total_human)
    if human_host != 0:
        print elmseq + '\t' + str(bird_flu_frac/human_flu_frac) + '\t' + str(chicken_host/human_host)


# human_freqs = get_freqs(freq_file_human)
# chicken_freqs = get_freqs(freq_file_chicken)
# v1 = []
# v2 = []
# for key in human_freqs:
#     if key in uniq and key in chicken_freqs:
#         hf = human_freqs[key]
#         cf = chicken_freqs[key]
#         if abs(hf-cf) > 5:
#             sp = key.split(':')
#             elmseq = sp[1] + ':' + sp[2]
#             hhf = 0
#             chf = 0
#             if elmseq in human_host_freqs:
#                 hhf = human_host_freqs[elmseq]
#             if elmseq in chicken_host_freqs:
#                 chf = chicken_host_freqs[elmseq]
#             if chf != 0 and hhf != 0:
#                 v1.append(hf/cf)
#                 v2.append(float(hhf)/float(chf))
#                 print key
# print utils_stats.pcc(v1,v2)
# with open('working/vals', 'w') as f:
#     for a,b in zip(v1,v2):
#         f.write(str(a) + '\t' + str(b) + '\n')



