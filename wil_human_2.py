import utils_stats
elm_ls = [[1, 'YILD'], [2, 'YCIL'], [2, 'YLMA'], [3, 'YCRA'], [11, 'YIMR'], [18, 'YLLT'], [43, 'YTLD'], [72, 'YTID'], [88, 'YVSM'], [96, 'YLLA'], [101, 'YLEK'], [104, 'YIMK'], [114, 'YCVL'], [115, 'YFTA'], [116, 'YINT'], [116, 'YVDG'], [116, 'YVRT']]
cut = 90

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
