import utils_stats
elm_ls = [[1, 'YFDD'], [1, 'YFTT'], [1, 'YIST'], [1, 'YIXK'], [1, 'YLAK'], [1, 'YLPD'], [1, 'YLTA'], [1, 'YLXA'], [1, 'YLYV'], [1, 'YTID'], [1, 'YTMK'], [1, 'YTTG'], [1, 'YVDR'], [1, 'YVEE'], [1, 'YVEG'], [1, 'YVRP'], [2, 'YCRA'], [2, 'YIDG'], [2, 'YIKY'], [2, 'YILD'], [2, 'YIRT'], [2, 'YTIR'], [4, 'YCIL'], [6, 'YIVK'], [7, 'YIAS'], [8, 'YIIK'], [11, 'YLLA'], [12, 'YLMS'], [18, 'YFTS'], [81, 'YVVK'], [859, 'YIEG'], [860, 'YLLS'], [1318, 'YLMA'], [2182, 'YVDG'], [2183, 'YCVL'], [2184, 'YVRT'], [2185, 'YTLD'], [2186, 'YLEK'], [2191, 'YFTA'], [2192, 'YIMK'], [2193, 'YINT']]

virus = []
nonvirus = []
found_seqs = {}
with open('results/elmdict_H_sapiens.txt') as f:
    for line in f:
        elm, seq, count, frac_st = line.strip().split('\t')
        if elm == 'LIG_SH2_STAT5':
            appended = False
            for elm_count, elm_seq in elm_ls:
                if seq == elm_seq:
                    found_seqs[seq] = True
                    if elm_count > 2000:
                        virus.append(float(frac_st))
                    else:
                        nonvirus.append(float(frac_st))
                    appended = True
                    break
            if not appended:
                nonvirus.append(float(frac_st))
for count, seq in elm_ls:
    if not seq in found_seqs and count > 2000:
        virus.append(float(0))
print utils_stats.wilcox_gtr(virus, nonvirus)
