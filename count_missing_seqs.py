"""I need to know what sequences
   correspond to seq_ids in the
   roundup files.

   When I grab files from NCBI,
   I get different GIs than what
   I started with.

   This script counts the seq_ids
   for which I do not have sequences.

   I found that just fish & macaque were
   giving me problems.
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

roundup_file = sys.argv[1]
fasta_dir = sys.argv[2]

seqs_w_fasta = defaultdict(dict)
seqs_in_clusters = defaultdict(dict)
hosts = {}
# load roundup file
with open(roundup_file) as f:
    for line in f:
        (cluster, host, seq_id) = line.strip().split('\t')
        seqs_in_clusters[host][seq_id] = True
        hosts[host] = True
# load fasta file
for host in hosts:
    with open(fasta_dir
              + name2label[host] + '.fa') as f:
        for line in f:
            if line[0] == '>':
                ID = line.split('|')[1].strip()
                seqs_w_fasta[host][ID] = True

# count missing seqs
for host in seqs_in_clusters:
    print host, len(seqs_in_clusters[host]), len(set(seqs_in_clusters[host]) &
                                                 set(seqs_w_fasta[host]))
