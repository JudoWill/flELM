""" Each human ELM/seq pair has a normalized fraction 
    (normalized to that ELM). I'll split the ELM/seq 
    pairs into virus and non-virus, and use the 
    one-sided wilcoxon test to see if virus 
    ELM/seq pairs have lower fractions.
"""
import utils, os, utils_stats
from global_settings import *
from local_settings import *

species2dict = {}
flu2dict = {}

for g in ['H_sapiens']:#GENOMES:
    species2dict[g] = utils.get_seq2count_dict('results/elmdict_'
                                               + g + '.redo',
                                               float(0))

for flu in ['HIV']:#FLU_NAMES:
    flu2dict[flu] = utils.get_seq2count_dict('results/hiv_freq',
                                             float(0.1))

virus_like = []
non_virus = []
host = 'H_sapiens'
virus = 'HIV'
for elm in species2dict[host]:
    if elm in flu2dict[virus]:
        for seq in species2dict[host][elm]:
            if seq in flu2dict[virus][elm]:
                #if flu2dict[virus][elm][seq] > float(.05):
#                virus_like.append([elm+':'+seq,species2dict[host][elm][seq]])
                virus_like.append(species2dict[host][elm][seq])
                #else:
                #    non_virus.append(species2dict[host][elm][seq])
            else:
                non_virus.append(species2dict[host][elm][seq])
    else:
        for seq in species2dict[host][elm]:
            non_virus.append(species2dict[host][elm][seq])
print utils_stats.wilcox_gtr(virus_like, non_virus)
with open('virus', 'w') as f:
    for item in virus_like:
        f.write('blank\t' + str(item) + '\n')
with open('nonvirus', 'w') as f:
    for item in non_virus:
        f.write('blank\t' + str(item) + '\n')
print len(virus_like), len(non_virus)

    
