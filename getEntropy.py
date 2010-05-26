import utils, os, sys
from global_settings import *
from local_settings import *

def no_args():
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

def one_arg(afile):
    d = utils.get_seq2count_dict(afile, float(0))
    entropy = utils.get_species_entropy(d)
    for elm in entropy:
        print elm + '\t' + str(entropy[elm])

def main():
    if len(sys.argv) == 1:
        no_args()
    else:
        one_arg(sys.argv[1])

if __name__ == '__main__': main()
