"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import itertools, sys, os, utils, random, global_settings, numpy, math
import Bio.Cluster
from collections import defaultdict

def print_results(elm, clusters, overlap):
    """Print out clustering of ELM sequences.
       Ignore sequences in overlap that cannot
       be assigned."""

    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                print('%s\t%d\t%s' %
                      (elm, cluster, seq))

# this comes from my scratch experiments
human_distance_file = '../../scratch/human_flu_distances'
chicken_distance_file = '../../scratch/chicken_flu_distances'
both_distance_file = 'working/runs/Jun24/closest_dis'

distance_file = both_distance_file

# map each flu elmSeq to a host elmSeq
host_flu_elmSeq_mapping = {}
flu_host_elmSeq_mapping = {}
host_flu_elmSeq_best_dis = {}
flu_host_elmSeq_best_dis = {}
mapped = {}

with open(distance_file) as f:
    for line in f:
        (elm, flu_seq, host_seq, distance) = line.strip().split('\t')
        dis = int(distance)
        host_elmSeq = elm + ':' + host_seq
        flu_elmSeq = elm + ':' + flu_seq
        
        if dis < 5:
            if elm not in flu_host_elmSeq_best_dis:
                flu_host_elmSeq_best_dis[elm] = {}
                flu_host_elmSeq_mapping[elm] = {}

            if flu_elmSeq not in flu_host_elmSeq_best_dis[elm]:
                flu_host_elmSeq_best_dis[elm][flu_elmSeq] = {}
                flu_host_elmSeq_mapping[elm][flu_elmSeq] = {}

            if host_elmSeq not in flu_host_elmSeq_best_dis[elm][flu_elmSeq]:
                flu_host_elmSeq_best_dis[elm][flu_elmSeq][host_elmSeq] = dis
                flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

            # if dis < flu_host_elmSeq_best_dis[elm][flu_elmSeq][host_elmSeq]:
            #     # replace
            #     flu_host_elmSeq_best_dis[elm][flu_elmSeq][host_elmSeq] = dis
            #     flu_host_elmSeq_mapping[elm][flu_elmSeq] = {}
            #     flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True
            # elif dis == flu_host_elmSeq_best_dis[elm][flu_elmSeq][host_elmSeq]:
                # append
            flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

for elm in flu_host_elmSeq_mapping:
    dis_mat = []
    for seq1,seq2 in itertools.combinations(flu_host_elmSeq_mapping[elm], 2):
        match = float(len(set(seq1) & set(seq2)))
        dis = numpy.average([match/float(x) for x in [len(seq1), len(seq2)]])
        dis_mat.append(dis)
    mat = numpy.array(dis_mat)
    num_clusters = 10
    percent_overlap = float(1)
    while num_clusters > 1 and percent_overlap > float(.1):
        if len(flu_host_elmSeq_mapping[elm]) > num_clusters:
            ans, error, nfound = Bio.Cluster.kmedoids(mat, nclusters=num_clusters, 
                                                      npass=10)
            clusters = defaultdict(dict)
            total = {}
            for flu_seq, cluster_id in zip(flu_host_elmSeq_mapping[elm],
                                           ans):
                clusters[cluster_id][flu_seq] = True
                total[flu_seq] = True
                for host_seq in flu_host_elmSeq_mapping[elm][flu_seq]:
                    clusters[cluster_id][host_seq] = True  
                    total[host_seq] = True
            overlap = {}
            for cluster1, cluster2 in itertools.combinations(clusters, 2):
                for gene in set(clusters[cluster1]) & set(clusters[cluster2]):
                    overlap[gene] = True
            percent_overlap = float(len(overlap))/float(len(total))
            if percent_overlap < float(.1):
                print_results(elm, clusters, overlap)
                break
            else:
                num_clusters -= 1
        else:
            break
sys.exit(0)

for flu_elmSeq in flu_host_elmSeq_mapping:
    for host_elmSeq in flu_host_elmSeq_mapping:
        if not host_elmSeq in host_flu_elmSeq_mapping:
            host_flu_elmSeq_mapping[host_elmSeq] = {}
        host_flu_elmSeq_mapping[host_elmSeq][flu_elmSeq] = True
        mapped[flu_elmSeq] = True

def count_0s(ls):
    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

hosts = ('H_sapiens', 'Gallus_gallus')
flus = ('human', 'chicken')

# count elm:seq occurence
flu_counts = {}
pre_flu_counts = {}
host_counts = {}
all_elmSeqs = {}
for flu in flus:
    pre_flu_counts[flu] = utils.get_flu_counts('results/' + flu + '.H5N1.elms', 
                                               global_settings.FLU_PROTEINS)

# sample protein sequenes from chicken
new_chicken_counts = {}
new_chicken_proteins = {}
for protein in proteins:
    c = len(pre_flu_counts['human'][protein].keys())
    new_chicken_proteins[protein] = random.sample(pre_flu_counts['chicken'][protein], c)

for protein in new_chicken_proteins:
    new_chicken_counts[protein] = {}
    for seq in new_chicken_proteins[protein]:
        new_chicken_counts[protein][seq] = pre_flu_counts['chicken'][protein][seq]

flu_counts['human'] = utils.count_flu(pre_flu_counts['human'], all_elmSeqs)
flu_counts['chicken'] = utils.count_flu(new_chicken_counts, all_elmSeqs)

for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('working/runs/Jun24/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elmSeq in host_flu_elmSeq_mapping:
                 for a_flu_elmSeq in host_flu_elmSeq_mapping[elmSeq]:
                     host_counts[host][a_flu_elmSeq] += int(count)
            # else:
            #     host_counts[host][elmSeq] += int(count)
            #     all_elmSeqs[elmSeq] = True

flu_vecs = utils.mk_count_vecs(flu_counts, mapped)
flu_dists = utils.mk_count_dists(flu_vecs)
host_vecs = utils.mk_count_vecs(host_counts, mapped)
host_dists = utils.mk_count_dists(host_vecs)

js_distances = defaultdict(dict)
for host in hosts:
    for flu in flus:
        js_dis = utils.jensen_shannon_dists(host_dists[host],
                                            flu_dists[flu])
        js_distances[host][flu] = js_dis
        print host, flu, js_dis

def print_it(name, vec):
    print name, float(count_0s(vec))/float(len(vec))

print_it('chicken_flu', flu_vecs['chicken'])
print_it('human flu', flu_vecs['human'])
print_it('H sapiens', host_vecs['H_sapiens'])
print_it('Gallus gallus', host_vecs['Gallus_gallus'])
