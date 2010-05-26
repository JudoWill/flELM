from collections import defaultdict
import utils, sys
elms = defaultdict(utils.init_zero)
total = 0
trg_elm = 'MOD_SUMO'
with open(sys.argv[1]) as f:
    for line in f:
        [elm, seq, count, frac] = line.strip().split('\t')
        if elm == trg_elm:
            c = seq[0]
            elms[c] += int(count)
            total += int(count)
for c in elms:
    print c + '\t' + str(float(elms[c])/float(total))
