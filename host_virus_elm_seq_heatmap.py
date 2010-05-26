""" Compare the ELM sequence utilization for hosts
    to that of viruses.
"""
import sys, os, local_settings, utils, utils_motif, utils_graph, random
from collections import defaultdict

suffix = sys.argv[1]

def getConservedELMs(virus, subtypes):
    ls = [utils_motif.annotation2protein(os.path.join(local_settings.RESULTSDIR,
                                                      virus + '.' + subtype 
                                                      + '.elms.70.controled'),
                                         {'ELM':True}) 
          for subtype in subtypes[virus]]
    return utils_graph.intersectLists(ls)

hosts = {'H_sapiens':'M', 'Sus_scrofa':'P', 'Gallus_gallus':'C', 
         'Equus_caballus':'H', 'Taeniopygia_guttata':'F'}
viruses = {'human':'m', 'swine':'p', 'equine':'h',
           'duck':'d', 'chicken':'c'}
subtypes = {'human':('H1N1', 'H3N2', 'H5N1'),
            'swine':('H1N1', 'H3N2'),
            'equine':('H3N8',),
            'chicken':('H5N1', 'H9N2'),
            'duck':('H5N1', 'H9N2')}

virus2conservedELMs = {}
all_elms = {}
for virus in viruses:
    virus2conservedELMs[virus] = getConservedELMs(virus, subtypes)
    for elm in virus2conservedELMs[virus]:
        all_elms[elm] = True

# load ELM seq fractions
host2elmFreqs = {}
virus2elmFreqs = {}
host2elm2seqs = defaultdict(dict)
virus2elm2seqs = defaultdict(dict)
for host in hosts:
    host2elmFreqs[host] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                'elmdict_' + host + suffix),
                                                   float(0.001))
    for elm in host2elmFreqs[host]:
        for seq in host2elmFreqs[host][elm]:
            host2elm2seqs[elm][seq] = True

for virus in viruses:
        virus2elmFreqs[virus] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                      'flu_elmdict_' + virus), float(0.1))
        for elm in virus2elmFreqs[virus]:
            for seq in virus2elmFreqs[virus][elm]:
                virus2elm2seqs[elm][seq] = True

virus_total = 0
virus_match = 0
for elm in virus2elm2seqs:
    for seq in virus2elm2seqs[elm]:
        virus_total += 1
        if seq in host2elm2seqs[elm]:
            virus_match += 1
print('matching virus ELM sequences')
print('%.10f\t%d%d' %
      (float(virus_match)/float(virus_total),
       virus_match, virus_total))

host_total = 0
host_match = 0
for elm in host2elm2seqs:
    if elm in virus2elm2seqs:
        for seq in host2elm2seqs[elm]:
            host_total += 1
            if seq in virus2elm2seqs[elm]:
                host_match += 1
print('matching host ELM sequences')
print('%.10f\t%d%d' %
      (float(host_match)/float(host_total),
       host_match, host_total))
            
                
                      
