""" Get flu sequences """
import sys, os
from utils import GetFluSeqs, init_zero
from global_settings import *
from local_settings import *
from collections import defaultdict

short_name = sys.argv[1]
strain = sys.argv[2]

for short_name in sys.argv[1:]:
    flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name],
                         strain=['influenza a virus'],
                         source=['h5n1'])
    seqs = {}
    for seq, protein, gb in flu_gen:
        if protein in VIRUS_PROTEINS:
            seqs[gb + '.' + protein.replace(' ','')] = seq

    with open('results/' + short_name + '.H5N1.fa', 'w') as f:
        for protein in seqs:
            f.write('>' + protein + '\n')
            f.write(seqs[protein] + '\n')
