"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import sys, os, utils, random, global_settings, utils_graph
from collections import defaultdict

# this comes from my scratch experiments
human_distance_file = '../../scratch/human_flu_distances'
chicken_distance_file = '../../scratch/chicken_flu_distances'
both_distance_file = 'working/runs/Jun24/closest_dis'

def print_it(name, vec):
    print name, float(count_0s(vec))/float(len(vec))

def print_results(elm, clusters, overlap):
    """Print out clustering of ELM sequences.
       Ignore sequences in overlap that cannot
       be assigned."""

    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                print('%s\t%d\t%s' %
                      (elm, cluster, seq))

mapping = utils.get_clusters(both_distance_file)    
hosts = global_settings.TEST_GENOMES
#all_elmSeqs = {}
flus = ('human',)
flu_counts = {}
seen_seqs = {}
seen_seqs_ls = []
for flu in flus:
     pre = utils.get_flu_counts('results/' + flu + '.H5N1.elms', 
                                global_settings.FLU_PROTEINS)
     flu_counts_sampled = {}
     flu_proteins_sampled = {}
     protein_counts = []
     for protein in pre:
         protein_counts.append(len(pre[protein]))
     m = min(protein_counts)
     for protein in pre:
         if m == len(pre[protein]):
             flu_proteins_sampled[protein] = pre[protein].keys()
         else:
             flu_proteins_sampled[protein] = random.sample(pre[protein], m)
     for protein in flu_proteins_sampled:
         flu_counts_sampled[protein] = {}
         for sampled_protein in flu_proteins_sampled[protein]:
             flu_counts_sampled[protein][sampled_protein] = pre[protein][sampled_protein]
     seen_seqs[flu] = {}
     flu_counts[flu] = utils.count_flu(flu_counts_sampled, 
                                       mapping, seen_seqs[flu])
     seen_seqs_ls.append(seen_seqs[flu])
if len(seen_seqs_ls) > 1:
    all_elmSeqs = utils_graph.intersectLists(seen_seqs_ls)
else:
    all_elmSeqs = seen_seqs_ls[0]

host_counts = {}
found_seqs = {}
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    found_seqs[host] = {}
    with open('working/runs/Jun24/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm in mapping:
                if elmSeq in mapping[elm]:
                    key = mapping[elm][elmSeq]
                    if key in all_elmSeqs:
                        host_counts[host][key] += int(count)
                        found_seqs[host][key] = True
                    # else:
                    #     host_counts[host][key] += int(count)
                    #     found_seqs[host][key] = True
host_found_seqs = utils_graph.intersectLists([found_seqs['H_sapiens'],
                                              found_seqs['Gallus_gallus']])
use_seqs = utils_graph.intersectLists([all_elmSeqs, host_found_seqs])
        

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

# print_it('chicken_flu', flu_vecs['chicken'])
# print_it('human flu', flu_vecs['human'])
# print_it('H sapiens', host_vecs['H_sapiens'])
# print_it('Gallus gallus', host_vecs['Gallus_gallus'])

