import utils_motif, sys

conserved_file = sys.argv[1]
elm_file = sys.argv[2]
conserved = utils_motif.protein2annotation(conserved_file,
                                           {'ELM':True})
elms_pre = utils_motif.protein2annotation(elm_file,
                                          {'ELM':True})
elms = {}
for protein in elms_pre:
    vp = protein.split('.')[-1]
    if not vp in elms:
        elms[vp] = {}
    for elm in elms_pre[protein]:
        if not elm in elms[vp]:
            elms[vp][elm] = []
        for tri in elms_pre[protein][elm]:
            elms[vp][elm].append(tri)

elm2seq2count = {}
for vp in conserved:
    for elm in conserved[vp]:
        for [st, stp, seq] in elms[vp][elm]:
            if not elm in elm2seq2count:
                elm2seq2count[elm] = {}
            if not seq in elm2seq2count[elm]:
                elm2seq2count[elm][seq] = 0
            elm2seq2count[elm][seq] += 1
for elm in elm2seq2count:
    total = 0
    for seq in elm2seq2count[elm]:
        total += elm2seq2count[elm][seq]
    for seq in elm2seq2count[elm]:
        count = elm2seq2count[elm][seq]
        print(elm + '\t' + seq + '\t' + str(count) + '\t' + str(float(count)/float(total)))
