"""When simplifying the residue alphabet,
   I made so many ELM patterns that I had
   to split up the file. This combines the
   results.

   Enter results directory & maximum suffix """
import sys, os, global_settings, utils
from collections import defaultdict

dir = sys.argv[1] # working/Jul12/
end = int(sys.argv[2]) # 10

for g in ('H_sapiens', 'Gallus_gallus'):
    counts = defaultdict(utils.init_zero)
    for x in xrange(end+1):
        file = dir + 'elmdict_' + g + '.simple' + str(x)
        with open(file) as f:
            for line in f:
                str1, str2, count, freq = line.strip().split('\t')
                counts[str2] += int(count)
    new_file = dir + 'elmdict_' + g + '.simple'
    with open(new_file, 'w') as f:
        for seq in counts:
            f.write('%s\t%s\t%d\t%s\t\n' %
                    (seq, seq, counts[seq], '1'))
