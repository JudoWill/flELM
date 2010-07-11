""" The elmdicts for the simplified
    residue alphabet are formated:
    ELM:seq elm ...,
    and they need to be ELM seq

    This reformats them.

    Not done."""
import sys, utils
from collections import defaultdict

elmdict_file = sys.argv[1]
elm2seq2count = defaultdict(utils.init_zero)
elm2count = defaultdict(utils.init_zero)

with open(elmdict_file) as f:
    for line in f:
        (elmseq, seq, count, frq) = line.strip().split('\t')
        elm = elmseq.split(':')[0] + ':' + str(len(seq))
        elm2seq2count[elmseq + ':' + str(len(seq))] += int(count)
        elm2count[elm] += int(count)

for elmSeq in elm2seq2count:   
    elm, seq, length = elmSeq.split(':')
    count = elm2seq2count[elmSeq]
    freq = float(count)/float(elm2count[elm + ':' + length])
    print('%s\t%s\t%d\t%.10f' %
          (elm + ':' + length, 
           seq, count, freq))


    
