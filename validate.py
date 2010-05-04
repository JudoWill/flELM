""" Assume that ELM sequences in a virus
    will have higher frequencies in the 
    correct host than other hosts.

    Using horse, human, chicken, and pig,
    compare the correct hosts to all other
    species and find the # ELM/sequence 
    combinations for which the correct host has
    a higher ELM/sequence frequency.    
"""
import utils, os, random, utils_distance
from global_settings import *
from local_settings import *

conserved_hiv = {}
with open('results/HIV1.clean.elms.90.conserved') as f:
    for line in f:
        conserved_hiv[line.split('\t')[3]] = True

def get_dis(species2dict, flu2dict, species, flu):
    """ grab the sequence fraction """
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
                                                          flu + '.elms.90.freq.redo'),
                                             float(0.01))

tmp_input = 'tmp' + str(random.randint(0, 100))
with open(tmp_input,'w') as f:
    f.write('Host\tSpecies\tPercent\n')
    for host, flu in [ ['Equus_caballus', 'equine'],
                       ['Gallus_gallus', 'chicken'],
                       ['H_sapiens', 'human'],
                       ['Sus_scrofa', 'swine'] ]:
        this_dis = get_dis(species2dict, flu2dict, host, flu)
        
        for g in GENOMES:
            if g != host:
                better = 0
                dis = get_dis(species2dict, flu2dict, g, flu)
                total_elms = utils_distance.get_elements(this_dis, dis)
                diff = float(0)
                for elm_seq in total_elms:
                    elm, seq = elm_seq.split(':')
                    if elm in conserved_hiv:                        
                        if elm_seq in this_dis:
                            if elm_seq in dis:
                                diff += (this_dis[elm_seq] - dis[elm_seq])
                            else:
                                diff += this_dis[elm_seq]
                        else:
                            diff -= dis[elm_seq]
                f.write(host + '\t' + g + '\t'
                        + str(diff) + '\n')
    
