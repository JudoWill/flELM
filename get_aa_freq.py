import sys, utils
from collections import defaultdict
counts = defaultdict(utils.init_zero)
total = 0
with open(sys.argv[1]) as f:
    for line in f:
        if line[0] != '>':
            for chr in line.strip():
                if chr != 'X' and chr != 'U':
                    counts[chr] += 1
                    total += 1
for aa in counts:
    print aa + '\t' + str(float(counts[aa])/float(total))
