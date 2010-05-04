""" For conserved virus ELMs, find the cosine distance
    between sequence fraction vectors of virus and host.

    Updated to use HIV and HCV conserved ELMs.
"""
from scipy.spatial import distance
from numpy import *
import utils, local_settings, os, utils_motif, utils_graph, random, os

def getConservedELMs(virus, subtypes):
    ls = [utils_motif.annotation2protein(os.path.join(local_settings.RESULTSDIR,
                                                      virus + '.' + subtype 
                                                      + '.elms.70.controled'),
                                         {'ELM':True}) 
          for subtype in subtypes[virus]]
    return utils_graph.intersectLists(ls)

def getDistance(virus_d, host_d):
    seqs = utils_graph.unionLists([virus_d,
                                   host_d])
    host_v = []
    virus_v = []
    for seq in seqs:
        for v,d in ( (host_v, host_d),
                     (virus_v, virus_d) ):
            if seq in d:
                v.append(d[seq])
            else:
                v.append(0)
    return distance.cosine(host_v, virus_v)

hosts = {'H_sapiens':'M', 'Sus_scrofa':'P', 'Gallus_gallus':'C', 
         'Equus_caballus':'H', 'Taeniopygia_guttata':'F'}
viruses = {'human':'m', 'swine':'p', 'equine':'h',
           'duck':'d', 'chicken':'c', 'HIV':'a', 'HCV':'l'}
subtypes = {'human':('H1N1', 'H3N2', 'H5N1'),
            'swine':('H1N1', 'H3N2'),
            'equine':('H3N8',),
            'chicken':('H5N1', 'H9N2'),
            'duck':('H5N1', 'H9N2'),
            'HIV':('all',),
            'HCV':('all',)}

virus2conservedELMs = {}
for virus in viruses:
    virus2conservedELMs[virus] = getConservedELMs(virus, subtypes)

# load ELM seq fractions
host2elmFreqs = {}
virus2elmFreqs = {}
for host in hosts:
    host2elmFreqs[host] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                'elmdict_' + host + '.txt'),
                                                   float(0))
#tmp_input = 'tmp_input' + str(random.randint(0,100))
tmp_input = 'plots/for_aydin/cos_host_virus.tab'
with open(tmp_input, 'w') as f:
    f.write('Virus_Host\tELM\tDistance\n')
    for virus in viruses:
        virus2elmFreqs[virus] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                      'flu_elmdict_' + virus), float(0))
        for elm in virus2conservedELMs[virus]:
            if 'FAIL' not in elm:
                for host in hosts:
                    dis = getDistance(virus2elmFreqs[virus][elm], 
                                      host2elmFreqs[host][elm])
                    f.write('%s\t%s\t%.10f\n'
                            % (viruses[virus] + hosts[host], elm, dis))
out_file = 'plots/for_aydin/cos_dis_heatmap.png'
tmp_r = 'tmp_r' + str(random.randint(0,100))
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d,aes(Virus_Host,ELM)) + opts(axis.text.y = theme_blank()) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='red',high='steelblue')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')

