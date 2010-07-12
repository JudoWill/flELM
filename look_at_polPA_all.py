from collections import defaultdict
import os, utils

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
            if protein == 'nonstructural protein 1':
                seq = elm + ':' + str(len(elm.split(':')[1]))
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

notUniq_human = {}
for protein in human_elmseqs:
    for elm in human_elmseqs[protein]:
        key = protein + ':' + elm
        if key not in uniq_human:
            notUniq_human[key] = True

notUniq_bird = {}
for protein in bird_elmseqs:
    for elm in bird_elmseqs[protein]:
        key = protein + ':' + elm
        if key not in uniq_bird:
            notUniq_bird[key] = True

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
hosts = ('swine',)
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
sums = defaultdict(utils.init_zero)
totals = defaultdict(utils.init_zero)

for elmseq in uniq_bird:
    protein, elm, seq, length = elmseq.split(':')
    key = elm + ':' + seq + ':' + length
    if key in chicken_host_freqs and key in human_host_freqs:
        diff = chicken_host_freqs[key]-human_host_freqs[key]
        if diff > float(0):
            sums[elm] += 1
        totals[elm] += 1

full_sum = 0
full_total = 0
for elm in totals:
    #print elm, float(100)*float(sums[elm])/float(totals[elm])
    full_sum += sums[elm]
    full_total += totals[elm]
print full_sum, full_total-full_sum, float(100)*float(full_sum)/float(full_total)

print 'NOT UNIQ BIRD'
sums = defaultdict(utils.init_zero)
totals = defaultdict(utils.init_zero)

for elmseq in notUniq_bird:
    protein, elm, seq, length = elmseq.split(':')
    key = elm + ':' + seq + ':' + length
    if key in chicken_host_freqs and key in human_host_freqs:
        diff = chicken_host_freqs[key]-human_host_freqs[key]
        if diff > float(0):
            sums[elm] += 1
        totals[elm] += 1

full_sum = 0
full_total = 0
for elm in totals:
    #print elm, float(100)*float(sums[elm])/float(totals[elm])
    full_sum += sums[elm]
    full_total += totals[elm]
print full_sum, full_total-full_sum, float(100)*float(full_sum)/float(full_total)

print 'MAMMAL'
s = defaultdict(utils.init_zero)
t = defaultdict(utils.init_zero)
for elmseq in uniq_human:
    protein, elm, seq, length = elmseq.split(':')
    key = elm + ':' + seq + ':' + length
    if key in chicken_host_freqs and key in human_host_freqs:
        diff = human_host_freqs[key]-chicken_host_freqs[key]
        if diff > float(0):
            s[elm] += 1
        t[elm] += 1

full_sum = 0
full_total = 0
for elm in t:
    #print elm, float(100)*float(s[elm])/float(t[elm])
    full_sum += s[elm]
    full_total += t[elm]
print full_sum, full_total-full_sum, float(100) * float(full_sum)  / float(full_total)

print 'NOT UNIQ MAMMAL'
s = defaultdict(utils.init_zero)
t = defaultdict(utils.init_zero)
for elmseq in notUniq_human:
    protein, elm, seq, length = elmseq.split(':')
    key = elm + ':' + seq + ':' + length
    if key in chicken_host_freqs and key in human_host_freqs:
        diff = human_host_freqs[key]-chicken_host_freqs[key]
        if diff > float(0):
            s[elm] += 1
        t[elm] += 1

full_sum = 0
full_total = 0
for elm in t:
    #print elm, float(100)*float(s[elm])/float(t[elm])
    full_sum += s[elm]
    full_total += t[elm]
print full_sum, full_total-full_sum, float(100) * float(full_sum)  / float(full_total)
