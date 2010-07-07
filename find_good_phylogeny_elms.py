"""Which ELMs correctly separate birds
   and mammals?"""
import sys, utils, itertools, global_settings
from collections import defaultdict

def check_phylogeny(dis):
    """make sure that chicken & finch
       are closer to each other than
       any mammal"""

    birds = ('Gallus_gallus', 'Taeniopygia_guttata')
    bird_dis = dis['Gallus_gallus']['Taeniopygia_guttata']
    for host in dis:
        if host not in birds:
            for bird in birds:
                new_dis = dis[bird][host]
                if new_dis <= bird_dis:
                    return False
    return True

results_dir = sys.argv[1]
suffix = sys.argv[2]
elmfile = sys.argv[3]

elms = {}
with open(elmfile) as f:
    for line in f:
        elm, exp = line.strip().split('\t')
        elms[elm] = True

for elm in elms:
    counts = utils.count_host_elmSeqs(global_settings.TEST_GENOMES,
                                      False, {},
                                      results_dir, {elm:True}, suffix)
    all_elmSeqs = {}
    for host in counts:
        for elmSeq in counts[host]:
            all_elmSeqs[elmSeq] = True
    host_vecs = utils.mk_count_vecs(counts, all_elmSeqs)
    host_dists = utils.mk_count_dists(host_vecs)
    js = defaultdict(dict)
    for host1, host2 in itertools.combinations(host_dists, 2):
        dis = utils.jensen_shannon_dists(host_dists[host1],
                                         host_dists[host2])
        js[host1][host2] = dis
        js[host2][host1] = dis
    if check_phylogeny(js):
        print elm
