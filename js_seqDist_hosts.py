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

results_dir = sys.argv[1] # working/runs/Jun24/
out_file = sys.argv[2]

f = os.path.join(results_dir, 'test_host_seqs')
use_seqs = utils_graph.getNodes(f)
    
counts = utils.count_host_seqs(global_settings.PLT_GENOMES,
                               results_dir, use_seqs, '.init')

ls = []
for host in counts:
    ls.append(counts[host])
all_elmSeqs = {}
for host in counts:
    for elmSeq in counts[host]:
        all_elmSeqs[elmSeq] = True

host_vecs = utils.mk_count_vecs(counts, all_elmSeqs)
host_dists = utils.mk_count_dists(host_vecs)
utils_plot.phylogeny_js(os.path.join(results_dir,
                                     out_file), host_dists)
