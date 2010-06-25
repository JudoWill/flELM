"""I need to limit myself to clusters
   that have only one sequence from
   each organism, but have FASTA seqs
   for all organisms

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
              "Canis familiaris": 'Canis_familiaris'}

def getID(line):
    return line.split('|')[1].strip()

random.seed()
roundup_file = sys.argv[1] # results/Homo_Mus_Pan_Rat_Bos_Can_Gal_Tae_Dan_Mac.roundup.parsed
fasta_dir = sys.argv[2] # data/roundup_all/
outdir = sys.argv[3]

seqs_w_fasta = defaultdict(dict)
clusters = {}
hosts = {}
# load roundup file
with open(roundup_file) as f:
    for line in f:
        (cluster, host, seq_id) = line.strip().split('\t')
        if host in name2label:
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
                seqs_w_fasta[host][getID(line)] = True
# remove seq_ids from roundup that
# are not in fasta
for cluster in clusters:
    for host in clusters[cluster]:
        to_remove = {}
        for seq_id in clusters[cluster][host]:
            if seq_id not in seqs_w_fasta[host]:
                to_remove[seq_id] = True
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
    if len(set(hosts) & set(clusters[cluster])) != len(hosts):
        to_remove[cluster] = True
for rm in to_remove:
    del clusters[rm]
# remove clusters with more than
# one sequence for a host
to_remove = {}
for cluster in clusters:
    for host in clusters[cluster]:
        if len(clusters[cluster][host]) > 1:
            to_remove[cluster] = True
            break
for rm in to_remove:
    del clusters[rm]

# make new fasta from
# resulting clusters
seqs = defaultdict(dict)
for cluster in clusters:
    for host in clusters[cluster]:
        seq_id = clusters[cluster][host].keys()[0]
        seqs[host][seq_id] = True

for host in hosts:
    fasta_itr = utils.fasta_iter(fasta_dir 
                                 + name2label[host] + '.fa',
                                 getID)
    with open(outdir + name2label[host] + '.fa', 'w') as outf:
        for ID, seq in fasta_itr:
            if ID in seqs[host]:
                outf.write('>' + ID + '\n')
                outf.write(seq + '\n')
print start_clusters, len(clusters)
                
