""" Find the residue probabilities
    for each flu genome.
"""
import sys, os
from utils import GetFluSeqs, init_zero
from global_settings import *
from local_settings import *
from collections import defaultdict

for short_name in sys.argv[1:]:
    flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name])
    d = defaultdict(init_zero)
    total = 0
    for seq in flu_gen:
        for aa in seq:
            if aa != 'X' and aa != 'B' and aa != 'J' and aa != 'Z':
                d[aa] += 1
                total += 1
    outfile = os.path.join(RESULTSDIR, 'flu.' + short_name + '.aa_freq')
    with open(outfile, 'w') as f:
        for aa in d:
            f.write(aa + '\t' + str(float(d[aa])/float(total)) + '\n')


    

