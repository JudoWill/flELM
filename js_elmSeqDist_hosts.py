"""Use Jensen-Shannon divergence to make a dendrogram for eukaryotic hosts"""
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

species = global_settings.TEST_GENOMES
short_names = global_settings.ALIASES

results = 'working/runs/Jun24/'

# count elm:seq occurence
counts = {}
all_elmSeqs = {}
for host in species:
    counts[host] = defaultdict(utils.init_zero)
    with open(results + 'elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            all_elmSeqs[elmSeq] = True
            counts[host][elmSeq] += int(count)

host_vecs = utils.mk_count_vecs(counts, all_elmSeqs)
host_dists = utils.mk_count_dists(host_vecs)

tmp_input = 'tmp_data'
tmp_r = 'tmp_r' + str(random.randint(0,100))
tmp_labels = 'labels' + str(random.randint(0,100))
out_file = results + 'js_hosts_elmSeq_dendrogram.png'

js_distances = defaultdict(dict)
for host1, host2 in itertools.combinations(species, 2):
    js_dis = utils.jensen_shannon_dists(host_dists[host1],
                                        host_dists[host2])
    js_distances[host1][host2] = js_dis
    js_distances[host2][host1] = js_dis

with open(tmp_input, 'w') as f:    
    for host1 in species:
        line = ''
        for host2 in species:
            if host1 == host2:
                line += '0\t'
            else:
                line += str(js_distances[host1][host2]) + '\t'
        f.write(line.strip('\t') + '\n')

with open(tmp_labels, 'w') as f:
    f.write('\t'.join(species) + '\n')

with open(tmp_r, 'w') as f:
    f.write("source('funcs.R')\n")
    f.write("library('MASS')\n")
    f.write("d<-read.delim('"
            + tmp_input
            + "',header=FALSE,sep='\\t')\n")
    f.write('dist.r<-as.dist(d)\n')
    f.write("labels.d<-read.delim('"
            + tmp_labels
            + "',header=FALSE,sep='\\t')\n")
    f.write('labels<-as.matrix(labels.d)\n')
    f.write("h<-hclust(dist.r,method='average')\n")
    f.write("png('" + out_file + "')\n")
    f.write("plot(h,hang=-1,labels=labels[1,],main='Species Dendrogram')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')
os.system('rm ' + tmp_r + ' ' + tmp_labels + ' ' + tmp_input)
