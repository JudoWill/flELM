import utils_fasta, sys, utils, global_settings
from collections import defaultdict

aa = 0
for g in global_settings.GENOMES:
    fasta = utils_fasta.loadFASTA('data/' + g + '.fa')
    for p in fasta:
        aa += len(fasta[p])

elm2count = defaultdict(utils.init_zero)
for g in global_settings.GENOMES:
    with open('results/elmdict_' + g + '.txt') as f:
        for line in f:
            [elm, seq, count_st, freq] = line.strip().split('\t')
            elm2count[elm] += int(count_st)

with open('results/all_elm_aa_freq', 'w') as f:
    for elm in elm2count:
        v = float(elm2count[elm])/float(aa)
        f.write(elm + '\t' + str(v) + '\n')
