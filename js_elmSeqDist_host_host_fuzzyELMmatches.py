"""Use Jensen-Shannon divergence to find distances between chicken/human
   host/flu.  Since host/flu share few ELM sequences, only consider
   flu ELMs & find host ELMs that are closest in edit distance."""
import itertools, sys, os, utils, random, global_settings, numpy, math
import Bio.Cluster
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

def mk_mapping(elm, clusters, overlap, mapping):
    """Update a map of sequences to cluster.
       Ignore sequences in overlap that cannot
       be assigned."""

    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                mapping[elm][seq] = elm + ':' + str(cluster)

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

            if dis < 2:
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
        num_clusters = 5
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

def count_flu(protein2counts, all_elmSeqs):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                counts[elmSeq] += protein2counts[protein][seq][elmSeq]
                all_elmSeqs[elmSeq] = True
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
all_elmSeqs = {}
host_counts = {}
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm in mapping:
                if elmSeq in mapping[elm]:
                    key = mapping[elm][elmSeq]
                    all_elmSeqs[key] = True
                    host_counts[host][key] += int(count)
                    
host_vecs = mk_count_vecs(host_counts, all_elmSeqs)
host_dists = mk_count_dists(host_vecs)

tmp_input = 'tmp_data'
tmp_r = 'tmp_r' + str(random.randint(0,100))
tmp_labels = 'labels' + str(random.randint(0,100))
out_file = 'plots/for_aydin_2/roundup_all/js.dendrogram.fuzzy.png'

js_distances = defaultdict(dict)
for host1, host2 in itertools.combinations(hosts, 2):
    js_dis = utils.jensen_shannon_dists(host_dists[host1],
                                        host_dists[host2])
    js_distances[host1][host2] = js_dis
    js_distances[host2][host1] = js_dis

with open(tmp_input, 'w') as f:    
    for host1 in hosts:
        line = ''
        for host2 in hosts:
            if host1 == host2:
                line += '0\t'
            else:
                line += str(js_distances[host1][host2]) + '\t'
        f.write(line.strip('\t') + '\n')

with open(tmp_labels, 'w') as f:
    f.write('\t'.join(hosts) + '\n')

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

