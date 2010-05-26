import utils, os
from global_settings import *
from local_settings import *

species2dict = {}
for g in GENOMES:
    d = utils.calc_elm_frequency('results/elmdict_'
                                 + g + '.redo')
    with open(os.path.join(RESULTSDIR, g + '.redo.elm_freq'), 'w') as f:
        for elm in d:
            f.write(elm + '\t' + str(d[elm]) + '\n')

for flu in FLU_NAMES:
    d = utils.calc_elm_frequency('results/flu_elmdict_'
                                 + flu)
    with open(os.path.join(RESULTSDIR, 'flu_' + flu + '.elm_freq'), 'w') as f:
        for elm in d:
            f.write(elm + '\t' + str(d[elm]) + '\n')
