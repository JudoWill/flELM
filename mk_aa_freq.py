import utils_fasta, sys, utils
from collections import defaultdict

fasta = utils_fasta.loadFASTA(sys.argv[1])

aa = 0
for p in fasta:
    aa += len(fasta[p])

elm2count = defaultdict(utils.init_zero)
with open(sys.argv[2]) as f:
    for line in f:
        [elm, seq, count_st, freq] = line.strip().split('\t')
        elm2count[elm] += int(count_st)

with open(sys.argv[3], 'w') as f:
    for elm in elm2count:
        v = float(elm2count[elm])/float(aa)
        f.write(elm + '\t' + str(v) + '\n')
