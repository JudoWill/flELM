"""Use Jensen-Shannon divergence 
   to make a dendrogram for eukaryotic hosts.
   Choose to cluster the ELM sequences before
   calculating JS divergence.
   To skip clusteirng, enter NA as the first
   argument. Otherwise, enter a closest flu 
   distance file computed by flu_project_host_flu_closest.py.
"""
import itertools, sys, os, utils, random, global_settings, numpy, utils_plot
from collections import defaultdict

distance_file = sys.argv[1] # closest_dis
results_dir = sys.argv[2] # working/runs/Jun24/
out_file = sys.argv[3] # js_hosts_elmSeq_dendrogram_new.png
dis_cutoff_init = int(sys.argv[4]) # 2
dis_cutoff_meta = float(sys.argv[5]) # 3

do_clustering = True
if distance_file == 'NA':
    do_clustering = False

if do_clustering:
    dis_file = os.path.join(results_dir, distance_file)
    mapping = utils.get_clusters(dis_file, dis_cutoff_init,
                                 dis_cutoff_meta)
    
# count elm:seq occurence
counts = {}
all_elmSeqs = {}
for host in global_settings.TEST_GENOMES:
    counts[host] = defaultdict(utils.init_zero)
    elm_file = os.path.join(results_dir, 
                            'elmdict_' + host + '.init')
    with open(elm_file) as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if do_clustering:
                if elmSeq in mapping[elm]:
                    key = mapping[elm][elmSeq]
                    all_elmSeqs[key] = True
                    counts[host][key] += int(count)
            else:
                all_elmSeqs[elmSeq] = True
                counts[host][elmSeq] += int(count)

host_vecs = utils.mk_count_vecs(counts, all_elmSeqs)
host_dists = utils.mk_count_dists(host_vecs)
utils_plot.phylogeny_js(os.path.join(results_dir,
                                     out_file), host_dists)
