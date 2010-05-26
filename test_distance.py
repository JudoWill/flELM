import elm_hists, utils_plot, utils
import utils_distance
import itertools
from global_settings import *

# human = elm_hists.get_seq2count_dict('results/elmdict_H_sapiens.txt', 
#                                     float(.05))
# mouse = elm_hists.get_seq2count_dict('results/elmdict_M_musculus.txt', 
#                                     float(.05))
# monkey = elm_hists.get_seq2count_dict('results/elmdict_Macaca_mulatta.txt',
#                                      float(.05))
# print utils_distance.distance_species(human,
#                                       monkey)
# print utils_distance.distance_species(human,
#                                       mouse)
# print utils_distance.distance_species(monkey,
#                                       mouse)

species2dict = {}
virus2dict = {}

virus2dict['swineFlu'] = utils.get_seq2count_dict('results/flu_elmdict_swine',
                                                  float(.4))
virus2dict['chickenFlu'] = utils.get_seq2count_dict('results/flu_elmdict_chicken',
                                                    float(.4))
virus2dict['humanFlu'] = utils.get_seq2count_dict('results/flu_elmdict_human',
                                                  float(.4))
for g in ('H_sapiens', 'Gallus_gallus', 'Sus_scrofa'):
    species2dict[g] = utils.get_seq2count_dict_for_seqs('results/elmdict_'
                                                        + g + '.txt',
                                                        float(0),
                                                        virus2dict)
for v in virus2dict:
    species2dict[v] = virus2dict[v]

d = utils_distance.distance_matrix(species2dict)
elm_d = utils_distance.elm_distance_matrix(species2dict)
#for elm in elm_d:
#    for species_pair in elm_d[elm]:
#        print elm + '\t' + species_pair + '\t' + str(elm_d[elm][species_pair])
for s1, s2 in itertools.combinations(d.keys(), 2):
    print s1 + '\t' + s2 + '\t' + str(d[s1][s2])
utils_plot.distance_heatmap(elm_d, 'test.png')

