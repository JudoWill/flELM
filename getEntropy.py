import utils, os
from global_settings import *
from local_settings import *

species2dict = {}
for g in GENOMES:
    d = utils.get_seq2count_dict('results/elmdict_'
                                 + g + '.txt',
                                 float(0))
    entropy = utils.get_species_entropy(d)
    with open(os.path.join(RESULTSDIR, g + '.elm_entropy'), 'w') as f:
        for elm in entropy:
            f.write(elm + '\t' + str(entropy[elm]) + '\n')
for flu in FLU_NAMES:
    d = utils.get_seq2count_dict('results/flu_elmdict_'
                                 + flu,
                                 float(0))
    entropy = utils.get_species_entropy(d)
    with open(os.path.join(RESULTSDIR, 'flu_' + flu + '.elm_entropy'), 'w') as f:
        for elm in entropy:
            f.write(elm + '\t' + str(entropy[elm]) + '\n')
