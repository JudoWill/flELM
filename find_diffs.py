import utils

d = {}
for g in ('human', 'chicken', 'swine', 'equine'):
    d[g] = utils.get_seq2count_dict('results/' + g + '.elms.90.freq.redo', float(0))

top_seqs = {}
for g in d:
    for elm in d[g]:
        ls = []
        for seq in d[g][elm]:
            ls.append([d[g][elm][seq], seq])
        ls.sort()
        if not g in top_seqs:
            top_seqs[g] = {}
        top_seqs[g][elm] = ls[-1][1]
for elm in top_seqs['human']:
    if elm in top_seqs['swine']:
        if top_seqs['human'][elm] != top_seqs['swine'][elm] and top_seqs['human'][elm] != top_seqs['chicken'][elm] and top_seqs['human'][elm] != top_seqs['equine'][elm]:
            print elm + '\t' + top_seqs['human'][elm] + '\t' + top_seqs['swine'][elm] + '\t' + top_seqs['chicken'][elm] + '\t' + top_seqs['equine'][elm]
        
