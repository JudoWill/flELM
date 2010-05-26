import utils_stats
elm_ls = [[1, 'YFSA'], [1, 'YFSS'], [1, 'YIAS'], [1, 'YINA'], [1, 'YLLV'], [1, 'YVKT'], [1, 'YVVK'], [2, 'YLLS'], [4, 'YTIR'], [6, 'YLMA'], [7, 'YVSM'], [14, 'YLLT'], [20, 'YLPD'], [66, 'YCVL'], [100, 'YCIL'], [108, 'YIEG'], [155, 'YLLA'], [164, 'YLEK'], [164, 'YTLD'], [164, 'YVRT'], [166, 'YFTA'], [167, 'YINT'], [167, 'YVDG'], [168, 'YIMK']]
cut = 125

virus = []
nonvirus = []
found_seqs = {}
with open('results/elmdict_Sus_scrofa.txt') as f:
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
