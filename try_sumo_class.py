from collections import defaultdict
import utils, sys
elms = defaultdict(utils.init_zero)
total = 0
#trg_elm = ('MOD_SUMO', 'LIG_KEPE_1', 'LIG_KEPE_2', 'LIG_KEPE_3')
#trg_elm = ('LIG_SH3_1', 'LIG_SH3_2', 'LIG_SH3_3', 'LIG_SH3_4', 'LIG_SH3_5')
with open(sys.argv[1]) as f:
    for line in f:
        [elm, seq, count, frac] = line.strip().split('\t')
        if elm.find('LIG_SH2_') != -1:
            elms[elm] += int(count)
            total += int(count)
for c in elms:
    print c + '\t' + str(float(elms[c])/float(total))
