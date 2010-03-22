""" Assume that ELM sequences in a virus
    will have higher frequencies in the 
    correct host than other hosts.

    Using horse, human, chicken, and pig,
    compare the correct hosts to all other
    species and find the # ELM/sequence 
    combinations for which the correct host has
    a higher ELM/sequence frequency.    
"""
import utils, os
from global_settings import *
from local_settings import *

def get_dis(species2dict, flu2dict, species, flu):
    dis = {}
    for elm in flu2dict[flu]:
        for seq in flu2dict[flu][elm]:
            if seq in species2dict[species][elm]:
                dis[elm + ':' + seq] = species2dict[species][elm][seq]
    return dis

species2dict = {}
flu2dict = {}

for g in GENOMES:
    species2dict[g] = utils.get_seq2count_dict('results/elmdict_'
                                               + g + '.redo',
                                               float(0))

for flu in ['human', 'equine', 'chicken', 'swine']:
    flu2dict[flu] = utils.get_seq2count_dict(os.path.join(RESULTSDIR, 
                                                          'flu_elmdict_' + flu + '.redo'),
                                             float(0))

for host, flu in [ ['Equus_caballus', 'equine'],
                   ['Gallus_gallus', 'chicken'],
                   ['H_sapiens', 'human'],
                   ['Sus_scrofa', 'swine'] ]:
    this_dis = get_dis(species2dict, flu2dict, host, flu)
    
    for g in GENOMES:
        if g != host:
            better = 0
            match = 0
            dis = get_dis(species2dict, flu2dict, g, flu)
            for elm in this_dis:
                if elm in dis:
                    if this_dis[elm] > dis[elm]:
                        better += 1
                    match += 1
                else:
                    better += 1
            total = float(len(this_dis.keys()))
            #total = float(match)
            print('HOST:' + host + '\t' + g 
                  + '\t' + str(better) + '\t' 
                  + str(float(better)/total))
    
