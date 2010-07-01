""" Get flu sequences based on strain, host species, and year"""
import sys, os
from utils import GetFluSeqs, init_zero
from global_settings import *
from local_settings import *
from collections import defaultdict

short_name = sys.argv[1]
strain = sys.argv[2]
use_year = sys.argv[3] # enter NA to take all years
working_dir = sys.argv[4]

if use_year == 'NA':
    flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name],
                         strain=['influenza a virus'],
                         source=[strain.lower()])
else:
    flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name],
                         year=[use_year],
                         strain=['influenza a virus'],
                         source=[strain.lower()])

seqs = {}
for seq, protein, gb in flu_gen:
    if protein in FLU_PROTEINS:
        seqs[gb + '.' + protein] = seq

fname = os.path.join(working_dir, short_name + '.' 
                     + strain + '.' + use_year + '.fa')
with open(fname, 'w') as f:
    for protein in seqs:
        f.write('>' + protein + '\n')
        f.write(seqs[protein] + '\n')
fname = os.path.join(working_dir, short_name + '.' 
                     + strain + '.' + use_year + '.fa.proteins')
with open(fname, 'w') as f:
    for protein in seqs:
        f.write(protein + '\n')
