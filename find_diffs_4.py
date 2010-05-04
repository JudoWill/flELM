import utils_motif

human = utils_motif.protein2annotation('human.H1N1.elms',
                                       {'ELM':True})
human_conserved = utils_motif.protein2annotation('human.H1N1.elms.90',
                                                 {'ELM':True})
swine = utils_motif.protein2annotation('swine.H1N1.elms',
                                       {'ELM':True})
swine_conserved = utils_motif.protein2annotation('swine.H1N1.elms.90',
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
    #if len(ls) > 1:
    #    print ls[0], ls[1]
    
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
swine_freq = get_freq('results/elmdict_Sus_scrofa.redo')

for protein in human_d:
    if protein in swine_d:
        for elm in human_d[protein]:
            if elm in swine_d[protein]:
                human_c, best_human = get_best_seq(human_d[protein][elm])
                swine_c, best_swine = get_best_seq(swine_d[protein][elm])
                if best_human != best_swine:
                    if elm in human_freq:
                        if not best_human in human_freq[elm]:
                            print 'swine DIFF', protein, elm, 'ZERO'
                        else:
                            print 'swine DIFF', protein, elm, 'PRESENT'
                    else:
                        print 'swine DIFF', protein, elm, 'ZERO'
                else:
                     if elm in human_freq:
                         if not best_human in human_freq[elm]:
                             print 'swine SAME', protein, elm, 'ZERO'
                         else:
                             print 'swine SAME', protein, elm, 'PRESENT'
                     else:
                         print 'swine SAME', protein, elm, 'ZERO'
