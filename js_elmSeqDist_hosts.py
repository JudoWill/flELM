"""Use Jensen-Shannon divergence 
   to make a dendrogram for eukaryotic hosts.
   Choose to cluster the ELM sequences before
   calculating JS divergence.
   To skip clusteirng, enter NA as the first
   argument. Otherwise, enter a closest flu 
   distance file computed by flu_project_host_flu_closest.py.
"""
import itertools, sys, os, utils, random, global_settings, numpy, utils_plot, utils_graph
from collections import defaultdict

distance_file = sys.argv[1] # closest_dis
results_dir = sys.argv[2] # working/runs/Jun24/
out_file = sys.argv[3] # js_hosts_elmSeq_dendrogram_new.png
dis_cutoff_init = int(sys.argv[4]) # 2
dis_cutoff_meta = float(sys.argv[5]) # 3
use_elms_file = sys.argv[6]
suffix = sys.argv[7]

use_elms = {}
with open(use_elms_file) as f:
    for line in f:
        (elm, stuff) = line.strip().split('\t')
        use_elms[elm] = True

do_clustering = True
if distance_file == 'NA':
    do_clustering = False

if do_clustering:
    dis_file = os.path.join(results_dir, distance_file)
    mapping = utils.get_clusters(dis_file, dis_cutoff_init,
                                 dis_cutoff_meta)
else:
    mapping = {}
    
counts = utils.count_host_elmSeqs(global_settings.TEST_GENOMES,
                                  do_clustering, mapping,
                                  results_dir, use_elms, suffix)

ls = []
for host in counts:
    ls.append(counts[host])
all_elmSeqs = {}
#all_elmSeqs = utils_graph.intersectLists(ls)
for host in counts:
    for elmSeq in counts[host]:
        all_elmSeqs[elmSeq] = True

host_vecs = utils.mk_count_vecs(counts, all_elmSeqs)
host_dists = utils.mk_count_dists(host_vecs)
utils_plot.phylogeny_js(os.path.join(results_dir,
                                     out_file), host_dists)
