""" How many of the ELM hits (ignore what sequence was matched)
    are not expected by chance?  Report # as fraction of 'real' hits.
"""
import sys, math, utils
from collections import defaultdict
from global_settings import *

def prob_of_seq(seq, aa2prob):
    p = float(1)
    for chr in seq:
        p *= aa2prob[chr]
    return p

aa2prob = {}
with open(sys.argv[1]) as f:
    for line in f:
        [aa, prob] = line.strip().split('\t')
        aa2prob[aa] = float(prob)

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
with open(sys.argv[3]) as f:
    for line in f:
        [elm, seq, count, frac] = line.strip().split('\t')
        if seq.find('X') == -1 and seq.find('U') == -1 and seq.find('B') == -1 and seq.find('J') == -1 and seq.find('Z') == -1:
            p = int(math.ceil(prob_of_seq(seq, aa2prob)*float(total_aa)))
            elm_expected[elm] += p
            elm_obs[elm] += int(count)
for elm in elm_obs:
    if elm_obs[elm] > elm_expected[elm]:
        val = float(elm_obs[elm]-elm_expected[elm])/float(elm_obs[elm])
        print('%s\t%.10f' %
              (elm, val))
                
                
        
