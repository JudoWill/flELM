"""Host and flu share few ELM sequences,
   so I'm making clusters of ELM sequences
   to be able to make comparisons between
   host and flu.

   I'll start each flu ELM sequence in a
   different cluster.  I'll add to these
   clusters any host sequence less than X
   edit distances away. Then I'll cluster
   these clusters choosing k as large as
   possible to limit the overlap between
   clusters to X%.
   
   Then I will add back the overlap
   sequences to clusters keeping the
   average distance below some threshold.

   I will check the clustering to make sure
   that:
   1) hosts produce the correct phylogeny
   2) at least X percent of the flu vectors
      are filled
"""
import Bio.Cluster, Levenshtein, utils_graph
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

distance_file = '../../scratch/test_genomes_distances'

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
        if best_cluster[0] < float(2):
            mapping[elm][elmSeq] = elm + ':' + str(best_cluster[1])
        else: # make new cluster
                mapping[elm][elmSeq] = elm + ':' + str(len(new_clusters)+1)
                new_clusters[len(new_clusters)+1][elmSeq] = True

def mk_mapping(elm, clusters, overlap, mapping, unmapped):
    """Update a map of sequences to cluster."""

    new_clusters = defaultdict(dict)
    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                mapping[elm][seq] = elm + ':' + str(cluster)
                new_clusters[cluster][seq] = True
    fix_overlap(elm, mapping, overlap, new_clusters)
    unmapped_formatted = {}
    for u in unmapped:
        unmapped_formatted[elm + ':' + u] = True
    fix_overlap(elm, mapping, unmapped_formatted, new_clusters)

def get_initial_clusters(distance_file):
    """Make a cluster for each flu sequence.
       Place in the cluster any host sequence
       below some threshold."""

    # map each flu elmSeq to a host elmSeq
    flu_host_elmSeq_mapping = {}
    total = defaultdict(dict)
    mapped = defaultdict(dict)
    with open(distance_file) as f:
        for line in f:
            (elm, flu_seq, host_seq, distance) = line.strip().split('\t')
            dis = int(distance)
            host_elmSeq = elm + ':' + host_seq
            flu_elmSeq = elm + ':' + flu_seq

            total[elm][host_seq] = True
            total[elm][flu_seq] = True

            if dis < 2:
                if elm not in flu_host_elmSeq_mapping:
                    flu_host_elmSeq_mapping[elm] = {}

                if flu_elmSeq not in flu_host_elmSeq_mapping[elm]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq] = {}

                if host_elmSeq not in  flu_host_elmSeq_mapping[elm][flu_elmSeq]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

                flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True
                mapped[elm][host_seq] = True
                mapped[elm][flu_seq] = True

    unmapped = {}
    for elm in total:
        if elm in mapped:
            unmapped[elm] = set(total[elm]) - set(mapped[elm])
        else:
            unmapped[elm] = total[elm].keys()
                                
    return (flu_host_elmSeq_mapping, unmapped)

def get_meta_clusters(flu_host_elmSeq_mapping, unmapped):
    """Cluster flu sequence clusters based on 
       overlapping host sequences.
       Placed unmapped sequences into clusters
       as well."""

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
                                                          npass=20)
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
                    mk_mapping(elm, clusters, overlap, mapping, unmapped[elm])
                    break
                else:
                    num_clusters -= 1
            else:
                break
    return mapping

def get_clusters():
    """Return map of sequences to cluster name"""

    (flu_host_elmSeq_mapping, unmapped) = get_initial_clusters(distance_file)
    mapping = get_meta_clusters(flu_host_elmSeq_mapping, unmapped)
    return mapping

def count_0s(ls):
    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

def mk_dis_mat(strings):
    """Make distance matrix"""

    dis_mat = []
    for s1 in strings:
        for s2 in strings:
            if s1 == s2:
                break
            else:
                dis = Levenshtein.distance(s1, s2)
                dis_mat.append(dis)
    return numpy.array(dis_mat)

def check_clusters(cluster2seq, host_counts, flu_counts, elm):
    """make sure flu sequences are not in clusters alone"""

    host_seqs = {}
    for host in host_counts:
        for elmSeq in host_counts[host]:
            helm,hseq = elmSeq.split(':')
            if helm == elm:
                host_seqs[hseq] = True
    flu_seqs = {}
    for flu in flu_counts:
        for elmSeq in flu_counts[flu]:
            helm,hseq = elmSeq.split(':')
            if helm == elm:
                flu_seqs[hseq] = True

    for cluster in cluster2seq:
        found_host = False
        found_flu = False
        for seq in cluster2seq[cluster]:

            if seq in host_seqs:
                found_host = True
            if seq in flu_seqs:
                found_flu = True
        if found_flu and not found_host:
            return False
    return True

def count_flu(protein2counts, all_elmSeqs):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                elm, sequence = elmSeq.split(':')
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
host_counts = {}
found_seqs = []
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    found_seqs.append({})
    with open('results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            if elm in mapping:
                if elmSeq in mapping[elm]:
                    key = mapping[elm][elmSeq]
                    host_counts[host][key] += int(count)
                    found_seqs[-1][key] = True
               #  else:
            #         host_counts[host][elmSeq] += int(count)
            #         found_seqs[-1][elmSeq] = True
            # else:
            #     host_counts[host][elmSeq] += int(count)
            #     found_seqs[-1][elmSeq] = True

use_seqs = utils_graph.intersectLists(found_seqs)
host_vecs = mk_count_vecs(host_counts, use_seqs)
host_dists = mk_count_dists(host_vecs)

tmp_input = 'tmp_data'
tmp_r = 'tmp_r' + str(random.randint(0,100))
tmp_labels = 'labels' + str(random.randint(0,100))
out_file = 'working/try_clusters.png'

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
os.system('rm ' + ' '.join([tmp_labels, tmp_input, tmp_r]))
