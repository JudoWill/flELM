import utils_stats
elm_ls = [[1, 'YIIK'], [1, 'YIVK'], [1, 'YLDK'], [1, 'YTIR'], [2, 'YLMA'], [3, 'YLLV'], [5, 'YIEG'], [10, 'YVNT'], [15, 'YTID'], [28, 'YVSM'], [129, 'YLLA'], [143, 'YLLT'], [260, 'YTLD'], [266, 'YINT'], [271, 'YVRT'], [273, 'YCVL'], [274, 'YLEK'], [275, 'YFTA'], [275, 'YIMK'], [277, 'YVDG']]
cut = 200

virus = []
nonvirus = []
found_seqs = {}
with open('results/elmdict_Gallus_gallus.txt') as f:
    for line in f:
        elm, seq, count, frac_st = line.strip().split('\t')
        if elm == 'LIG_SH2_STAT5':
            appended = False
            for elm_count, elm_seq in elm_ls:
                if seq == elm_seq:
                    found_seqs[seq] = True
                    if elm_count > cut:
                        virus.append(float(frac_st))
                    else:
                        nonvirus.append(float(frac_st))
                    appended = True
                    break
            if not appended:
                nonvirus.append(float(frac_st))
for count, seq in elm_ls:
    if not seq in found_seqs and count > cut:
        virus.append(float(0))
print utils_stats.wilcox_gtr(virus, nonvirus)
