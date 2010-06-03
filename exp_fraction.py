""" How many of the ELM hits are expected by chance /
    residues in proteome.
"""
import sys, math, utils
from collections import defaultdict
from global_settings import *

def prob_of_seq(seq, aa2prob):
    p = float(1)
    for c in seq:
        p *= aa2prob[c]
    return p

def prob_of_seq_trans(seq, aa2prob, trans_prob):
    p = aa2prob[seq[0]]
    for i in xrange(len(seq)-1):
        p *= trans_prob[seq[i]][seq[i+1]]
    return p

aa2prob = {}
with open(sys.argv[1]) as f:
    for line in f:
        [aa, prob] = line.strip().split('\t')
        aa2prob[aa] = float(prob)

trans_prob = defaultdict(dict)
with open(sys.argv[4]) as f:
    for line in f:
        [aa, prob] = line.strip().split('\t')
        c1 = aa[0]
        c2 = aa[1]
        trans_prob[c1][c2] = float(prob)

total_aa = 0
if sys.argv[2].find('fa') != -1:
    with open(sys.argv[2]) as f:
        for line in f:
            if line[0] != '>':
                total_aa += len(line.strip())
else:
     flu_gen = utils.GetFluSeqs(organism = FLU_NAMES[sys.argv[2]])
     for seq in flu_gen:
         total_aa += len(seq)

elm_expected = defaultdict(utils.init_zero)
elm_obs = defaultdict(utils.init_zero)
real = defaultdict(utils.init_zero)
with open(sys.argv[3]) as f:
    for line in f:
        [elm, seq, count, frac] = line.strip().split('\t')
        if seq.find('X') == -1 and seq.find('U') == -1 and seq.find('B') == -1 and seq.find('J') == -1 and seq.find('Z') == -1:
            p = int(math.floor(prob_of_seq(seq, aa2prob)*float(total_aa)))
            p_trans = int(math.floor(prob_of_seq_trans(seq, aa2prob, trans_prob)*float(total_aa)))
            elm_expected[elm] += p
            #print p, p_trans
            if int(count) > p:
                real[elm] += int(count)-p
        else:
            real[elm] += int(count)
        elm_obs[elm] += int(count)

for elm in elm_obs:
#    if elm_obs[elm] > elm_expected[elm]:
    if real[elm]:
        val = float(real[elm])/float(elm_obs[elm])
        print('%s\t%.10f' %
          (elm, val))
#    else:
#        val = float(0)
    
                
                
        
