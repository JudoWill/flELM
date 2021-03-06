# it is miss named now
from collections import defaultdict
import utils, utils_stats

def get_annotations(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            protein, elm = line.strip().split('\t')
            d[protein][elm] = True
    return d

def get_uniq(uniq, ls1, ls2):
    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[protein + ':' + elm] = True

def get_freqs(afile):
    freqs = {}
    with open(afile) as f:
        for line in f:
            protein, elmseq, freq = line.strip().split('\t')
            key = ':'.join([protein, elmseq])
            freqs[key] = float(freq)
    return freqs

def get_host_freqs(afile):
    freqs = {}
    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            freqs[elm + ':' + seq] = float(freq)
    return freqs

human_elmseqs_file = 'working/Jul1_year/mammal_elmseqs'
bird_elmseqs_file = 'working/Jul1_year/bird_elmseqs'

human_elmseqs = get_annotations(human_elmseqs_file)
bird_elmseqs = get_annotations(bird_elmseqs_file)

uniq = {}
get_uniq(uniq, human_elmseqs, bird_elmseqs)
get_uniq(uniq, bird_elmseqs, human_elmseqs)
print len(uniq)

# print H5N1 frequencies
freq_file_human = 'working/Jul1_year/human.H5N1.2009.elmseqs.conservation'
freq_file_chicken = 'working/Jul1_year/human.H1N1.2009.elmseqs.conservation'

human_host_file = 'working/Jun29/elmdict_H_sapiens.init'
chicken_host_file = 'working/Jun29/elmdict_Gallus_gallus.init'
human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)

human_freqs = get_freqs(freq_file_human)
chicken_freqs = get_freqs(freq_file_chicken)
v1 = []
v2 = []
for key in human_freqs:
    if key in uniq and key in chicken_freqs:
        hf = human_freqs[key]
        cf = chicken_freqs[key]
        if abs(hf-cf) > 5:
            sp = key.split(':')
            elmseq = sp[1] + ':' + sp[2]
            hhf = 0
            chf = 0
            if elmseq in human_host_freqs:
                hhf = human_host_freqs[elmseq]
            if elmseq in chicken_host_freqs:
                chf = chicken_host_freqs[elmseq]
            if chf != 0 and hhf != 0:
                v1.append(hf/cf)
                v2.append(float(hhf)/float(chf))
                print key
print utils_stats.pcc(v1,v2)
with open('working/vals', 'w') as f:
    for a,b in zip(v1,v2):
        f.write(str(a) + '\t' + str(b) + '\n')



