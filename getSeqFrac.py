import utils, os
from global_settings import *
from local_settings import *

species2dict = {}
for g in GENOMES:
    with open(os.path.join(RESULTSDIR, g + '.elm_seq_frac'), 'w') as f:
        with open('results/elmdict_'
                  + g + '.txt') as f2:
            for line in f2:
                [elm, seq, count, frac] = line.strip().split('\t')
                f.write(elm+':'+seq + '\t' + frac + '\n')

for flu in FLU_NAMES:
    with open(os.path.join(RESULTSDIR, 'flu_' + flu + '.elm_seq_frac'), 'w') as f:
        with open('results/flu_elmdict_'
                  + flu) as f2:
            for line in f2:
                [elm, seq, count, frac] = line.strip().split('\t')
                f.write(elm+':'+seq + '\t' + frac + '\n')
