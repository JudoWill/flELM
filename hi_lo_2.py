"""Do sequences uniq to either flu
   show a preference for hi/low
   distributions"""
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
birds = ('chicken','duck')
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
strains = ('H5N1', 'H1N1', 'H3N2', 'H3N8')
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

uniq_human = {}
uniq_bird = {}
get_uniq(uniq_human, human_elmseqs, bird_elmseqs)
get_uniq(uniq_bird, bird_elmseqs, human_elmseqs)
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

human_host_file = 'working/Jun29/elmdict_H_sapiens.init'
chicken_host_file = 'working/Jun29/elmdict_Gallus_gallus.init'
human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)

if sys.argv[1] == 'bird':
    with open('working/hi', 'w') as hi:
        for elmseq in uniq_bird:
            #'X' not in elmseq and 'J' not in elmseq and 'B' not in elmseq and
            if sum(seq_percents_bird[elmseq])/float(total_bird) - sum(seq_percents_human[elmseq])/float(total_human) > float(20):
                dels = elmseq.split(':')
                new_seq = utils.mk_sub(dels[2])
                key = dels[1] + ':' + dels[2]
                #frac = 0
                if key in chicken_host_freqs:
                    frac = chicken_host_freqs[key]
                    hi.write(str(frac) + '\n')

    with open('working/low', 'w') as low:
        for elmseq in uniq_human:
            # 'X' not in elmseq and 'J' not in elmseq and 'B' not in elmseq and
            if sum(seq_percents_human[elmseq])/float(total_human) - sum(seq_percents_bird[elmseq])/float(total_bird) > float(30):
                 dels = elmseq.split(':')
                 new_seq = utils.mk_sub(dels[2])
                 key = dels[1] + ':' + dels[2]
                 if key in chicken_host_freqs:
                     frac = chicken_host_freqs[key]
                     low.write(str(frac) + '\n')
elif sys.argv[1] == 'human':
    with open('working/hi', 'w') as hi:
        for elmseq in uniq_human:
            if sum(seq_percents_human[elmseq])/float(total_human) - sum(seq_percents_bird[elmseq])/float(total_bird) > float(30):
                dels = elmseq.split(':')
                new_seq = utils.mk_sub(dels[2])
                key = dels[1] + ':' + dels[2]
                #frac = 0
                if key in human_host_freqs:
                    frac = human_host_freqs[key]
                    cfrac = human_host_freqs[key]
                    hi.write(str(frac) + '\t' + str(cfrac) + '\n')
    with open('working/low', 'w') as low:
        for elmseq in uniq_bird:
            if sum(seq_percents_bird[elmseq])/float(total_bird) - sum(seq_percents_human[elmseq])/float(total_human) > float(40):
                 dels = elmseq.split(':')
                 new_seq = utils.mk_sub(dels[2])
                 key = dels[1] + ':' + dels[2]
                 if key in human_host_freqs:
                     frac = human_host_freqs[key]
                     cfrac = human_host_freqs[key]
                     low.write(str(frac) + '\t' + str(cfrac) + '\n')

   


