import utils_motif, utils_graph

use_elms = utils_graph.getNodes('use_elms')

human = utils_motif.protein2annotation('human.H1N1.elms',
                                       {'ELM':True})
human_conserved = utils_motif.protein2annotation('human.H1N1.elms.90',
                                                 {'ELM':True})
swine = utils_motif.protein2annotation('swine.H1N1.elms',
                                       {'ELM':True})
swine_conserved = utils_motif.protein2annotation('swine.H1N1.elms.90',
                                                 {'ELM':True})

def get_entropy(afile):
    entropy = {}
    with open(afile) as f:
        for line in f:
            [elm, entropy_st] = line.strip().split('\t')
            if not elm in entropy:
                entropy[elm] = {}
            entropy[elm] = float(entropy_st)
    return entropy

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
human_entropy = get_entropy('results/H_sapiens.elm_entropy')
swine_entropy = get_entropy('results/Sus_scrofa.elm_entropy')

for protein in human_d:
    if protein in swine_d:
        for elm in human_d[protein]:
            if elm in swine_d[protein]:
                human_c, best_human = get_best_seq(human_d[protein][elm])
                swine_c, best_swine = get_best_seq(swine_d[protein][elm])
                h_entropy = human_entropy[elm]
                s_entropy = swine_entropy[elm]
                if elm in use_elms:
                    if best_human != best_swine and abs(human_c-swine_c):
                        print 'DIFF\t' + protein + '\t' + elm + '\t' + str(h_entropy) + '\t' + str(s_entropy)
                    else:
                        print 'SAME\t' + protein + '\t' + elm + '\t' + str(h_entropy) + '\t' + str(s_entropy)
                    
        
