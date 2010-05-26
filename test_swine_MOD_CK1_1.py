import utils

def get_counts(afile):
    d = {}
    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            if elm == 'MOD_CK1_1':
                d[seq] = float(freq)
    return d

flu_counts = get_counts('results/flu_elmdict_human')
host_counts = get_counts('results/elmdict_Sus_scrofa.redo')
print utils.klDistance(flu_counts, host_counts)
