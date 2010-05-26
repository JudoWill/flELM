""" I need to find ELM sequences that differ across
    species hosts and viruses.
"""
import sys, utils, utils_distance, itertools, utils_plot
from global_settings import *

virus2dict = {}
virus2dict['chicken'] = utils.get_seq2count_dict('results/flu_elmdict_swine',
                                                  float(0.05))
virus2dict['swine'] = utils.get_seq2count_dict('results/flu_elmdict_chicken',
                                                    float(0.05))
virus2dict['human'] = utils.get_seq2count_dict('results/flu_elmdict_human',
                                                  float(0.05))

genomes = ('H_sapiens', 'Gallus_gallus', 'Sus_scrofa')
species2dict = {}
for g in genomes:
    species2dict[g] = utils.get_seq2count_dict('results/elmdict_'
                                               + g + '.txt',
                                               float(0.01))

use_elms = {}
for elm in virus2dict['human']:
    if elm in virus2dict['chicken'] and elm in virus2dict['swine']:
        distance_is_0 = False
        for v1, v2 in itertools.combinations(virus2dict.keys(), 2):
            distance = utils_distance.distance_elms(virus2dict[v1][elm],
                                                    virus2dict[v2][elm])
            if distance == float(0):
                distance_is_0 = True
        if not distance_is_0:
            use_elms[elm] = True

host_elms = {}
for elm in species2dict['H_sapiens']:
    if elm in species2dict['Gallus_gallus'] and elm in species2dict['Sus_scrofa']:
        distance_is_0 = False
        for h1, h2 in itertools.combinations(genomes, 2):
            distance = utils_distance.distance_elms(species2dict[h1][elm],
                                                    species2dict[h2][elm])
            if distance == float(0):
                distance_is_0 = True
        if not distance_is_0:
            host_elms[elm] = True

right_count = {'chicken':0, 'human':0, 'swine':0}
# wrong_count = {'chicken':0, 'human':0, 'swine':0}
# eq_count = {'chicken':0, 'human':0, 'swine':0}
v2h = {'chicken':'Gallus_gallus',
       'human':'H_sapiens',
       'swine':'Sus_scrofa'}
for v in virus2dict:
    for elm in use_elms:
        if elm in host_elms:
            dis = {}
            for g in genomes:
                dis[g] = utils_distance.distance_elms(species2dict[g][elm],
                                                      virus2dict[v][elm])
            is_less = True
            is_gtr = True
            for g in genomes:
                if g != v2h[v]:
                    if dis[v2h[v]] > dis[g]:#it is not less than at least one
                        is_less = False
            if is_less:
                right_count[v] += 1
            new_d = {}
            for h in species2dict:
                new_d[h] = species2dict[h]
            for v in virus2dict:
                new_d[v] = virus2dict[v]
            utils_plot.elm_host_barplot(new_d, 
                                        elm, elm + '.png')


for v in right_count:
    print v, right_count[v]

            
#for elm in use_elms:
#    utils_plot.elm_host_barplot(virus2dict, elm, elm + '.png')

# v = 'human'
# for elm in use_elms:
#     if elm in host_elms:
#         dis = {}
#         for g in genomes:
#             dis[g] = utils_distance.distance_elms(species2dict[g][elm],
#                                                   virus2dict[v][elm])
#         if dis['H_sapiens'] >= dis['Gallus_gallus'] and dis['H_sapiens'] >= dis['Sus_scrofa']:
#             print elm

# for k in virus_different:
#     [v1, v2, elm] = k.split(':')
#     seqs = {}
#     for v in virus2dict.keys():
#         for seq in virus2dict[v][elm]: seqs[seq] = True
#     new_species_dict = {}
#     for seq in seqs:
#         for g in genomes:
#             if elm in species2dict[g]:
#                 if seq in species2dict[g][elm]:
#                     if not g in new_species_dict:
#                         new_species_dict[g] = {}
#                     if not elm in new_species_dict[g]:
#                         new_species_dict[g][elm] = {}
#                     new_species_dict[g][elm][seq] = species2dict[g][elm][seq]
#     if len(new_species_dict.keys()) != 0:
#         utils_plot.elm_host_barplot(new_species_dict, elm, elm + '.host.png')
#         utils_plot.elm_host_barplot(virus2dict, elm, elm + '.png')
#         sys.exit(0)

# v = 'chickenFlu'
# for elm in virus2dict[v]:
#     for seq in virus2dict[v][elm]:
#         for g in genomes:
#             if elm in species2dict[g]:
#                 if seq in species2dict[g][elm]:
#                     print g + '\t' + elm + '\t' + seq + '\t' + str(species2dict[g][elm][seq]) + '\t' + str(virus2dict[v][elm][seq])



