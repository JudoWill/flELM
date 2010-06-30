"""There is little overlap between
   host and flu sequences. To solve
   this problem, I'm substiting some
   residues with symbols for their 
   properties"""
import sys, global_settings, utils
from collections import defaultdict

elmdict_file = sys.argv[1]
elm2seq2count = defaultdict(utils.init_zero)
with open(elmdict_file) as f:
    for line in f:
        (elm, seq, count, frq) = line.strip().split('\t')
        if 'X' not in seq:
            new_seq = utils.mk_sub(seq)
            elm2seq2count[elm+':'+new_seq] += int(count)
for elmSeq in elm2seq2count:   
    elm, seq = elmSeq.split(':')
    count = elm2seq2count[elmSeq]
    print('%s\t%s\t%d\t0' %
          (elm, seq, count))
    
