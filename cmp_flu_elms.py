import os
from collections import defaultdict

def get_elm_seqs(uniq_file, use_protein):
    elm2seq = defaultdict(dict)
    with open(uniq_file) as f:
        for line in f:
            protein, elmseq, elm, seq = line.strip().split('\t')
            if protein == use_protein:
                elm2seq[elm][seq] = True
    return elm2seq

protein = 'matrix protein 1'
dir = 'working/Jul12/'

mammal = get_elm_seqs(os.path.join(dir, 'mammal_uniq'),
                      protein)
bird = get_elm_seqs(os.path.join(dir, 'bird_uniq'),
                      protein)

for elm in mammal:
    if elm in bird:
        for mseq in mammal[elm]:
            for bseq in bird[elm]:
                print elm + '\t' + mseq + '\t' + bseq

