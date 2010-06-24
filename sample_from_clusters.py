"""I need to sample sequences from
   host clusters so that a single
   cluster does not dominate the
   host feature vectors.

   This script does that.
"""
import sys, random, utils
from collections import defaultdict

name2label = {"Homo sapiens": 'H_sapiens',
              "Mus musculus": 'M_musculus',
              "Rattus norvegicus": 'R_norvegicus',
              "Pan troglodytes": 'Pan_troglodytes',
              "Bos taurus": 'Bos_taurus',
              "Taeniopygia guttata": 'Taeniopygia_guttata',
              "Gallus gallus": 'Gallus_gallus',
              "Canis familiaris": 'Canis_familiaris',
              'Macaca mulatta': 'Macaca_mulatta',
              'Danio rerio': 'D_rerio'}

random.seed()
roundup_file = sys.argv[1]
fasta_dir = sys.argv[2]

seqs_w_fasta = defaultdict(dict)
clusters = {}
hosts = {}
# load roundup file
with open(roundup_file) as f:
    for line in f:
        (cluster, host, seq_id) = line.strip().split('\t')
        if cluster not in clusters:
            clusters[cluster] = {}
        if host not in clusters[cluster]:
            clusters[cluster][host] = {}
        clusters[cluster][host][seq_id] = True
        hosts[host] = True
start_clusters = len(clusters)
# load fasta file
for host in hosts:
    with open(fasta_dir
              + name2label[host] + '.fa') as f:
        for line in f:
            if line[0] == '>':
                seqs_w_fasta[host][line.split('|')[1].strip()] = True
# remove seq_ids from roundup that
# are not in fasta
for cluster in clusters:
    for host in clusters[cluster]:
        to_remove = {}
        for seq_id in clusters[cluster][host]:
            if seq_id not in seqs_w_fasta[host]:
                to_remove[seq_id] = True
            else:
                print seq_id, seqs_w_fasta[host].keys()[0]
        for rm in to_remove:
            del clusters[cluster][host][rm]
# remove hosts with no seq_ids
for cluster in clusters:
    to_remove = {}
    for host in clusters[cluster]:
        if not len(clusters[cluster][host]):
            to_remove[host] = True
    for rm in to_remove:
        del clusters[cluster][rm]
# remove clusters missing hosts
to_remove = {}
for cluster in clusters:
    if len(hosts) != len(clusters[cluster]):
        to_remove[cluster] = True
for rm in to_remove:
    del clusters[rm]
print start_clusters, len(clusters)
