""" The elmdicts for the simplified
    residue alphabet are formated:
    ELM:seq elm ...,
    and they need to be ELM seq

    This reformats them."""
import sys, utils
from collections import defaultdict

elmdict_file = sys.argv[1]
elm2seq2count = defaultdict(utils.init_zero)
elm2count = defaultdict(utils.init_zero)

with open(elmdict_file) as f:
    for line in f:
        (elmseq, seq, count, frq) = line.strip().split('\t')
        elm = elmseq.split(':')[0]
        elm2seq2count[elmseq] += int(count)
        elm2count[elm] += int(count)

for elmSeq in elm2seq2count:   
    elm, seq = elmSeq.split(':')
    count = elm2seq2count[elmSeq]
    freq = float(count)/float(elm2count[elm])
    print('%s\t%s\t%d\t%.10f' %
          (elm, seq, count, freq))


    
