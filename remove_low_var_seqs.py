""" Remake elmdicts using sequences
    with variance across hosts
    that meets some threshold.
"""
import sys, global_settings
import numpy
from collections import defaultdict

threshold = float(sys.argv[1])

def get_seqs(afile, elm_seqs):
    d = defaultdict(dict)
    counts = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            elm, seq, count, frac = line.strip().split('\t')
            d[elm][seq] = float(frac)
            counts[elm][seq] = int(count)
            elm_seqs[elm][seq] = True
    return (d, counts)

def renorm(elm, seq_counts, use_seqs, afile):
    total = 0
    for seq in use_seqs:
        if seq in seq_counts:
            total += seq_counts[seq]
    for seq in use_seqs:
        if seq in seq_counts:
            c = seq_counts[seq]
            afile.write('%s\t%s\t%d\t%.10f\n' %
                        (elm, seq, c, float(c)/float(total)))

freqs = {}
counts = {}
elm_seqs = defaultdict(dict)
for g in global_settings.GENOMES:
     freq_d, count_d = get_seqs('results/elmdict_'
                                + g + '.redo',
                                elm_seqs)
     freqs[g] = freq_d
     counts[g] = count_d

use_seqs = defaultdict(dict)
for elm in elm_seqs:
    for seq in elm_seqs[elm]:
        freq_v = []
        for g in freqs:
            if elm in freqs[g]:
                if seq in freqs[g][elm]:
                    freq_v.append(freqs[g][elm][seq])
                else:
                    freq_v.append(float(0))
            else:
                freq_v.append(float(0))
        v = numpy.var(numpy.array(freq_v))
        if v >= threshold:
            print('%s\t%s\t%.10f' %
                  (elm, seq, v))
            use_seqs[elm][seq] = True

for g in counts:
    with open('results/elmdict_'
              + g + '.varRedo', 'w') as f:
        for elm in counts[g]:
            renorm(elm, counts[g][elm], use_seqs[elm], f)
