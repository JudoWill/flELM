import utils_motif, sys

flu = sys.argv[1]
species = sys.argv[2]
strain = sys.argv[3]

human = utils_motif.protein2annotation('human.' + strain + '.elms',
                                       {'ELM':True})
human_conserved = utils_motif.protein2annotation('human.' + strain + '.elms.90',
                                                 {'ELM':True})
swine = utils_motif.protein2annotation(flu + '.' + strain + '.elms',
                                       {'ELM':True})
swine_conserved = utils_motif.protein2annotation(flu + '.' + strain + '.elms.90',
                                                 {'ELM':True})

def get_freq(afile):
    freq = {}
    with open(afile) as f:
        for line in f:
            #[elm_seq, freq_st] = line.strip().split('\t')
            #elm, seq = elm_seq.split(':')
            elm, seq, num, freq_st = line.strip().split('\t')
            if not elm in freq:
                freq[elm] = {}
            freq[elm][seq] = float(freq_st)
    return freq

def get_best_seq(seqs):
    ls = []
    for seq in seqs:
        ls.append([seqs[seq],seq])
    ls.sort()
    if len(ls) > 1:
        if ls[-1][0] == ls[-2][0]:
            print 'TIE'
    print ls
    return ls[-1]

def getProtein2elm2seq(elms):
    protein2elm2seq = {}
    for protein_id in elms:
        protein = protein_id.split('.')[1]
        if protein in human_conserved and protein in swine_conserved:
            for elm in elms[protein_id]:
                if elm in human_conserved[protein] and elm in swine_conserved[protein]:
                    if not protein in protein2elm2seq:
                        protein2elm2seq[protein] = {}
                    if not elm in protein2elm2seq[protein]:
                        protein2elm2seq[protein][elm] = {}
                    for [st, stp, seq] in elms[protein_id][elm]:
                        if not seq in protein2elm2seq[protein][elm]:
                            protein2elm2seq[protein][elm][seq] = 0
                        protein2elm2seq[protein][elm][seq] += 1
    return protein2elm2seq

human_d = getProtein2elm2seq(human)
swine_d = getProtein2elm2seq(swine)
human_freq = get_freq('results/elmdict_H_sapiens.redo')
swine_freq = get_freq('results/elmdict_' + species + '.redo')

for protein in ['polymerase PA']:#human_d:
    if protein in swine_d:
        for elm in ['LIG_SH2_STAT5']:#human_d[protein]:
            if elm in swine_d[protein]:
                human_c, best_human = get_best_seq(human_d[protein][elm])
                swine_c, best_swine = get_best_seq(swine_d[protein][elm])
                if best_human != best_swine:
                    human_score = {best_human:float(0),best_swine:float(0)}
                    swine_score = {best_human:float(0),best_swine:float(0)}

                    if best_human in human_freq[elm]:
                        human_score[best_human] = human_freq[elm][best_human]
                    if best_swine in human_freq[elm]:
                        human_score[best_swine] = human_freq[elm][best_swine]

                    if best_human in swine_freq[elm]:
                        swine_score[best_human] = swine_freq[elm][best_human]
                    if best_swine in swine_freq[elm]:
                        swine_score[best_swine] = swine_freq[elm][best_swine]
                    
                    result = 'NOT_COR'
                    if swine_score[best_swine] > swine_score[best_human] and human_score[best_human] > human_score[best_swine]:
                        result = 'POS_COR'
#                    else:
#                        print protein, elm, swine_score[best_swine], swine_score[best_human], human_score[best_human], human_score[best_swine]
                    print protein, elm, best_human, best_swine, result
        
