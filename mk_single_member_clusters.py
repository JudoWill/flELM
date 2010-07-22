"""I need to choose sequences from
   host clusters so that a single
   cluster does not dominate the
   host feature vectors.

   This script does that. Sequences
   are chosen to be similar in length
   to other host sequences in the cluster.
"""
import sys, random, utils, numpy
from collections import defaultdict

name2label = {"Homo sapiens": 'H_sapiens',
              "Gallus gallus": 'Gallus_gallus'}
              # "Mus musculus": 'M_musculus',
              # "Pan troglodytes": 'Pan_troglodytes',
              # "Sus scrofa": 'Sus_scrofa',
              # "Taeniopygia guttata": 'Taeniopygia_guttata',
              
              # "Equus caballus": 'Equus_caballus'}

def getID_local(line):
    return line.split('|')[1].strip()

def get_best_seq(cluster, host, seqs_w_fasta):
    """ For this host, find the seq that is
        closest in length to all other single
        member seqs in the cluster """

    lens = []
    for ahost in cluster:
        if len(cluster[ahost]) == 1:
            ID = cluster[ahost].keys()[0]
            lens.append(seqs_w_fasta[ahost][ID])
    if len(lens) == 0:
        return 'NA'

    dis_seq = []
    for candidate_seq in cluster[host]:
        candidate_len = seqs_w_fasta[host][candidate_seq]
        dis = numpy.average([abs(candidate_len-x) for x in lens])
        dis_seq.append([dis, candidate_seq])
    dis_seq.sort()
    if dis_seq[0][0] == dis_seq[1][0]:
        # sort alphabetically
        # and choose first one
        pair = [dis_seq[0][1], dis_seq[1][1]]
        pair.sort()
        ret_seq = pair[0]
    else:
        ret_seq = dis_seq[0][1]
    return ret_seq

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
# and count seq lengths
for host in hosts:
    afile = fasta_dir + name2label[host] + '.fa'
    for ID, seq in utils.fasta_iter(afile, getID_local):
        seqs_w_fasta[host][ID] = len(seq)
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

# make new fasta file by 
# sampling one species
# from each cluster

seqs = {}
for host in hosts:
    seqs[host] = {}

ignore_clusters = 0
for cluster in clusters:
    tmp_seqs = {}
    do_not_use = False
    for host in clusters[cluster]:
        if len(clusters[cluster][host]) > 1:
            seq = get_best_seq(clusters[cluster], 
                               host, seqs_w_fasta)
            if seq == 'NA':
                do_not_use = True
                break
        else:
            seq = clusters[cluster][host].keys()[0]
        tmp_seqs[host] = seq
    if not do_not_use:
        for host in tmp_seqs:
            seqs[host][tmp_seqs[host]] = True
    else:
        ignore_clusters += 1
print ignore_clusters, 'ignore clusters'
print len(seqs['Homo sapiens']), 'sequences per host'
print start_clusters, len(clusters), 'used clusters'

# writed sampled seq_id fasta
for host in hosts:
    fasta_itr = utils.fasta_iter(fasta_dir 
                                 + name2label[host] + '.fa',
                                 getID_local)
    with open(outdir + name2label[host] + '.fa', 'w') as outf:
        for ID, seq in fasta_itr:
            if ID in seqs[host]:
                outf.write('>' + ID + '\n')
                outf.write(seq + '\n')

                
