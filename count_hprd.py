import sys, utils_motif, utils
from collections import defaultdict
elm2proteins = utils_motif.annotation2protein(sys.argv[1], {'ELM':True})
for elm in elm2proteins:
    seqs = defaultdict(utils.init_zero)
    total = 0
    for protein in elm2proteins[elm]:
        for [st, stp, seq] in elm2proteins[elm][protein]:
            seqs[seq] += 1
            total += 1
    for seq in seqs:
        print elm + '\t' + seq + '\t' + str(seqs[seq]) + '\t' + str(float(seqs[seq])/float(total))
