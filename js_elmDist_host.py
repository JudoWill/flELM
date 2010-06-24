"""Count ELMs (not ELM sequences) on hosts
   and make phylogeny use Jensen-Shannon
   divergence."""
import itertools, utils, global_settings, os, random
from collections import defaultdict

def get_host_counts(ls_of_hosts):
    """Count total ELM hits for each
       host in ls_of_hosts.
       Return counts & ELMs found."""

    host_elmCounts = {}
    seen_elms = defaultdict(dict)
    for host in ls_of_hosts:
        host_elmCounts[host] = defaultdict(utils.init_zero)
        with open('results/roundup_all/elmdict_' + host + '.init') as f:
            for line in f:
                (elm, seq, count, fq) = line.strip().split('\t')
                host_elmCounts[host][elm] += int(count)
                seen_elms[elm][host] = True
    use_elms = {}
    for elm in seen_elms:
        if len(seen_elms[elm]) == len(ls_of_hosts):
            use_elms[elm] = True
    print len(use_elms)
    return (host_elmCounts, use_elms)

def mk_vec(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for this host's counts"""
    
    vec = []
    for elmseq in all_elmSeqs:
        if elmseq in counts:
            vec.append(counts[elmseq])
        else:
            vec.append(float(0))
    return vec

def mk_count_vecs(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for all hosts"""

    vecs = {}
    for host in counts:
        vecs[host] = mk_vec(counts[host],
                            all_elmSeqs)
    return vecs

def mk_count_dists(vecs):
    """change count vectors into distributions"""

    dists = {}
    for host in vecs:
        dists[host] = utils.getDistFromCount(vecs[host])
    return dists

hosts = global_settings.TEST_GENOMES#('H_sapiens', 'Macaca_mulatta', 'Pan_troglodytes', 'R_norvegicus', 'M_musculus')
host_elmCounts, elms = get_host_counts(hosts)
host_vecs = mk_count_vecs(host_elmCounts, elms)
host_dists = mk_count_dists(host_vecs)

tmp_input = 'tmp_data'
tmp_r = 'tmp_r' + str(random.randint(0,100))
tmp_labels = 'labels' + str(random.randint(0,100))
out_file = 'working/js_elm_host_phylogeny.png'

js_distances = defaultdict(dict)
for host1, host2 in itertools.combinations(host_dists, 2):
    js_dis = utils.jensen_shannon_dists(host_dists[host1],
                                        host_dists[host2])
    js_distances[host1][host2] = js_dis
    js_distances[host2][host1] = js_dis

with open(tmp_input, 'w') as f:    
    for host1 in host_dists:
        line = ''
        for host2 in host_dists:
            if host1 == host2:
                line += '0\t'
            else:
                line += str(js_distances[host1][host2]) + '\t'
        f.write(line.strip('\t') + '\n')

with open(tmp_labels, 'w') as f:
    f.write('\t'.join(host_dists) + '\n')

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
