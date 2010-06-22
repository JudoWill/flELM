"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import itertools, sys, os, utils, random, global_settings, numpy, math
import Bio.Cluster, Levenshtein, utils_graph
from collections import defaultdict

# this comes from my scratch experiments
human_distance_file = '../../scratch/human_flu_distances'
chicken_distance_file = '../../scratch/chicken_flu_distances'
both_distance_file = '../../scratch/test_genomes_distances'

def print_results(elm, clusters, overlap):
    """Print out clustering of ELM sequences.
       Ignore sequences in overlap that cannot
       be assigned."""

    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                print('%s\t%d\t%s' %
                      (elm, cluster, seq))

def fix_overlap(elm, mapping, overlap, new_clusters):
    """Assign sequences in overlap
       based on edit distance. Assume
       there are no ties to break."""

    for elmSeq in overlap:
        dis_cluster = []
        for cluster in new_clusters:
            dis = numpy.average([Levenshtein.distance(a_elmSeq.split(':')[1],
                                                      elmSeq.split(':')[1])
                                 for a_elmSeq in new_clusters[cluster]])
            dis_cluster.append([dis, cluster])
        dis_cluster.sort()
        best_cluster = dis_cluster[0]
        #print dis_cluster[0], dis_cluster[1]
        mapping[elm][elmSeq] = elm + ':' + str(best_cluster[1])

def mk_mapping(elm, clusters, overlap, mapping):
    """Update a map of sequences to cluster.
       Ignore sequences in overlap that cannot
       be assigned."""

    new_clusters = defaultdict(dict)
    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                mapping[elm][seq] = elm + ':' + str(cluster)
                new_clusters[cluster][seq] = True
    fix_overlap(elm, mapping, overlap, new_clusters)

def get_initial_clusters(distance_file):
    """Make a cluster for each flu sequence.
       Place in the cluster any host sequence
       below some threshold."""

    # map each flu elmSeq to a host elmSeq
    flu_host_elmSeq_mapping = {}

    with open(distance_file) as f:
        for line in f:
            (elm, flu_seq, host_seq, distance) = line.strip().split('\t')
            dis = int(distance)
            host_elmSeq = elm + ':' + host_seq
            flu_elmSeq = elm + ':' + flu_seq

            if dis < 3:
                if elm not in flu_host_elmSeq_mapping:
                    flu_host_elmSeq_mapping[elm] = {}

                if flu_elmSeq not in flu_host_elmSeq_mapping[elm]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq] = {}

                if host_elmSeq not in  flu_host_elmSeq_mapping[elm][flu_elmSeq]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

                flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

    return flu_host_elmSeq_mapping

def get_meta_clusters(flu_host_elmSeq_mapping):
    """Cluster flu sequence clusters based on 
       overlapping host sequences."""

    mapping = defaultdict(dict)
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
                ans, error, nfound = Bio.Cluster.kmedoids(mat, 
                                                          nclusters=num_clusters, 
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
                    #print_results(elm, clusters, overlap)
                    mk_mapping(elm, clusters, overlap, mapping)
                    break
                else:
                    num_clusters -= 1
            else:
                break
    return mapping

def get_clusters():
    """Return map of sequences to cluster name"""

    distance_file = both_distance_file
    flu_host_elmSeq_mapping = get_initial_clusters(distance_file)
    mapping = get_meta_clusters(flu_host_elmSeq_mapping)
    return mapping

def count_0s(ls):
    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

def count_flu(protein2counts, mapping, all_elmSeqs):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                elm = elmSeq.split(':')[0]
                if elm in mapping:
                    if elmSeq in mapping[elm]:
                        key = mapping[elm][elmSeq]
                        counts[key] += protein2counts[protein][seq][elmSeq]
                        all_elmSeqs[key] = True
    return counts

def get_flu_counts(afile, proteins):
    """Make protein_name -> seq_name -> elm_seq_counts"""

    counts = {}
    with open(afile) as f:
        for line in f:
            (protein, st, stp,
             elm, seq, junk) = line.strip().split('\t')
            name = protein.split('.')[-1]
            elmSeq = elm + ':' + seq
            if name in proteins:
                if name not in counts:
                    counts[name] = {}
                if protein not in counts[name]:
                    counts[name][protein] = {}
                if elmSeq not in counts:
                    counts[name][protein][elmSeq] = 0
                counts[name][protein][elmSeq] += 1
    return counts

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

mapping = get_clusters()    
hosts = global_settings.TEST_GENOMES
#all_elmSeqs = {}
flus = ('human', 'chicken')
proteins = ('hemagglutinin', 'neuraminidase', 'nucleocapsid protein',
            'matrix protein 1', 'nonstructural protein 1', 'matrix protein 2',
            'nonstructural protein 2', 'polymerase PA', 'polymerase PB2',
            'polymerase PB1', 'PB1-F2 protein')
flu_counts = {}
seen_seqs = {}
for flu in flus:
     pre = get_flu_counts('results/' + flu + '.H5N1.elms', 
                          proteins)
     seen_seqs[flu] = {}
     flu_counts[flu] = count_flu(pre, mapping, seen_seqs[flu])

all_elmSeqs = utils_graph.intersectLists([seen_seqs['human'],
                                          seen_seqs['chicken']])

host_counts = {}
found_seqs = {}
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    found_seqs[host] = True
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm in mapping:
                if elmSeq in mapping[elm]:
                    key = mapping[elm][elmSeq]
                    if key in all_elmSeqs:
                        host_counts[host][key] += int(count)
                        found_seqs[host][key] = True
host_found_seqs = utils_graph.intersectLists([found_seqs['H_sapiens'],
                                              found_seqs['Gallus_gallus']])
use_seqs = utils_graph.intersectLists([all_elmSeqs, host_found_seqs])
        

flu_vecs = mk_count_vecs(flu_counts, use_seqs)                   
host_vecs = mk_count_vecs(host_counts, use_seqs)
host_dists = mk_count_dists(host_vecs)
flu_dists = mk_count_dists(flu_vecs)

js_distances = defaultdict(dict)
for host in ('H_sapiens', 'Gallus_gallus'):
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

