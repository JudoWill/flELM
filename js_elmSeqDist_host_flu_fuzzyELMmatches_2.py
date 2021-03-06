"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import sys, os, utils, global_settings, utils_graph
from collections import defaultdict

cluster_distance_file = sys.argv[1] # NA for skip
elm_count_dir = sys.argv[2] # results/roundup_all/

do_clustering = True
if cluster_distance_file == 'NA':
    do_clustering = False
# this comes from my scratch experiments
#human_distance_file = '../../scratch/human_flu_distances'
#chicken_distance_file = '../../scratch/chicken_flu_distances'
#both_distance_file = 'working/runs/Jun24/closest_dis'
if do_clustering:
    f = os.path.join(elm_count_dir, cluster_distance_file)
    mapping = utils.get_clusters(f, 2.5, float(2.5))    
else:
    mapping = {}
hosts = global_settings.TEST_GENOMES
#all_elmSeqs = {}
flus = ('human',)
flu_counts = {}
seen_seqs = {}
seen_seqs_ls = []

for flu in flus:
    flu_elm_file = os.path.join('results/',
                                flu + '.H5N1.elms')
    utils.count_flu_sampled(flu, flu_elm_file, flu_counts,
                            seen_seqs, mapping, do_clustering)
    seen_seqs_ls.append(seen_seqs[flu])
if len(seen_seqs_ls) > 1:
    all_elmSeqs = utils_graph.unionLists(seen_seqs_ls)
else:
    all_elmSeqs = seen_seqs_ls[0]

host_counts = utils.count_host_elmSeqs(hosts, do_clustering, mapping,
                                       elm_count_dir)
                                                   
host_found_seqs = utils_graph.unionLists([host_counts['H_sapiens'],
                                          host_counts['Gallus_gallus']])
use_seqs = utils_graph.unionLists([all_elmSeqs, host_found_seqs])

flu_vecs = utils.mk_count_vecs(flu_counts, use_seqs)                   
host_vecs = utils.mk_count_vecs(host_counts, use_seqs)
host_dists = utils.mk_count_dists(host_vecs)
flu_dists = utils.mk_count_dists(flu_vecs)

js_distances = defaultdict(dict)
for host in ('H_sapiens', 'Gallus_gallus'):
    for flu in flus:
        js_dis = utils.jensen_shannon_dists(host_dists[host],
                                            flu_dists[flu])
        js_distances[host][flu] = js_dis
        print host, flu, js_dis

# utils.print_it('chicken_flu', flu_vecs['chicken'])
utils.print_it('human flu', flu_vecs['human'])
utils.print_it('H sapiens', host_vecs['H_sapiens'])
utils.print_it('Gallus gallus', host_vecs['Gallus_gallus'])

