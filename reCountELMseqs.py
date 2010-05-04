""" Only count ELM seqs based on
    the differences for specificed positions.
"""
import sys, utils
from collections import defaultdict

def getNewSeq(positions, seq):
    """make a sequence that only has
       positions specified in the ELM
       regular exp"""
    new_seq = ''
    for a_chr, pos in zip(seq, positions):
        if pos == '.':
            new_seq += '.'
        else:
            new_seq += a_chr
    return new_seq

# get positions that matter
elm2pos = {}
with open('exp_use') as f:
    for line in f:
        elm, pos = line.strip().split('\t')
        elm2pos[elm] = pos

elm2seq2count = {}
elm2count = defaultdict(utils.init_zero)
with open(sys.argv[1]) as f:
    for line in f:
        elm, seq, count, frac = line.strip().split('\t')
        if elm in elm2pos:
            new_seq = getNewSeq(elm2pos[elm], seq)
            if not elm in elm2seq2count:
                elm2seq2count[elm] = {}
            if not new_seq in elm2seq2count[elm]:
                elm2seq2count[elm][new_seq] = 0
            elm2seq2count[elm][new_seq] += int(count)
            elm2count[elm] += int(count)

for elm in elm2seq2count:
    for seq in elm2seq2count[elm]:
        count = elm2seq2count[elm][seq]
        print(elm + '\t' + seq + '\t' 
              + str(count) + '\t'
              + str(float(count)/float(elm2count[elm])))
        
