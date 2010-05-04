""" Get flu sequences """
import sys, os
from utils import GetFluSeqs, init_zero
from global_settings import *
from local_settings import *
from collections import defaultdict

short_name = sys.argv[1]
strain = sys.argv[2]

flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name],
                     strain=['influenza a virus'],
                     source=[strain.lower()])
seqs = {}
for seq, protein, gb in flu_gen:
    if protein in VIRUS_PROTEINS:
        seqs[gb + '.' + protein] = seq

with open('results/' + short_name + '.' + strain + '.fa', 'w') as f:
    for protein in seqs:
        f.write('>' + protein + '\n')
        f.write(seqs[protein] + '\n')
with open('results/' + short_name + '.' + strain + '.fa.proteins', 'w') as f:
    for protein in seqs:
        f.write(protein + '\n')
