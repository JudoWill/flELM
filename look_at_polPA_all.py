from collections import defaultdict
import os, utils

def get_host_freqs(afile):
    freqs = {}
    len_counts = defaultdict(utils.init_zero)
    with open(afile) as f:
        for line in f:
            elm_len, seq, count, freq = line.strip().split('\t')
            if 'LIG_MAPK_1' in elm_len:
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
            if protein == 'polymerase PA':
                if 'LIG_MAPK_1' in elm:
                    seq = elm.split(':')[1] + ':' + str(len(elm.split(':')[1]))
                    d[protein][seq] = True
    return d

def get_uniq(uniq, ls1, ls2):
    """Of the conserved ELM seqs, which ones are unique?"""

    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[protein + ':' + elm] = True

human_elmseqs_file = 'working/Jul7/mammal_elms'
bird_elmseqs_file = 'working/Jul7/bird_elms'

human_elmseqs = get_annotations(human_elmseqs_file)
bird_elmseqs = get_annotations(bird_elmseqs_file)

uniq_human = {}
uniq_bird = {}
get_uniq(uniq_human, human_elmseqs, bird_elmseqs)
get_uniq(uniq_bird, bird_elmseqs, human_elmseqs)

print len(uniq_human), len(uniq_bird)

dir = 'working/Jul7'
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
                                               'elms.conservation')))
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
                                               'elms.conservation')))
            if os.path.exists(file):
                total_human += 1
                get_freqs(file, seq_percents_human)

human_host_file = 'working/Jul7/elmdict_H_sapiens.RWlenInit'
chicken_host_file = 'working/Jul7/elmdict_Gallus_gallus.RWlenInit'
human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)

print 'BIRD'
sum = float(0)
for elmseq in uniq_bird:
    seq = elmseq.split(':')[1] + ':' + elmseq.split(':')[2]
    elmseq = 'LIG_MAPK_1' + ':' + seq
    diff = chicken_host_freqs[elmseq]-human_host_freqs[elmseq]
    print elmseq, diff
    if diff > float(0):
        sum += 1
print sum

print 'MAMMAL'
s = float(0)
for elmseq in uniq_human:
    seq = elmseq.split(':')[1] + ':' + elmseq.split(':')[2]
    elmseq = 'LIG_MAPK_1' + ':' + seq
    diff = human_host_freqs[elmseq]-chicken_host_freqs[elmseq]
    print elmseq, diff
    if diff > float(0):
        s += 1
print s
