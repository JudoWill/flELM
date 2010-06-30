""" What are the sequence differences for the
    given ELM for chicken & human H5N1?"""
import sys, utils, os, utils_graph

use_elm = sys.argv[1]

flu_counts = {}
seen_seqs = {}
seen_seqs_ls = []
flus = ('human', 'chicken')
for flu in flus:
    flu_elm_file = os.path.join('results/',
                                flu + '.H5N1.elms')
    utils.count_flu_sampled(flu, flu_elm_file, flu_counts,
                            seen_seqs, {}, False)
    use_seqs = {}
    for elmSeq in seen_seqs[flu]:
        elm, seq = elmSeq.split(':')
        if elm == use_elm:
            if flu_counts[flu][elmSeq] > 50:
                use_seqs[elmSeq] = True
    seen_seqs_ls.append(use_seqs)

# remove seqs seen less than 10x
#for flu in flu_counts:
#    for elmSeq in 

use_seqs_pre = utils_graph.intersectLists(seen_seqs_ls)

counts = utils.count_host_elmSeqs(('H_sapiens',),
                                  False, {},
                                  'working/Jun29/', {use_elm:True})
use_seqs = utils_graph.intersectLists([use_seqs_pre, counts['H_sapiens']])
host_vecs = utils.mk_count_vecs(counts, use_seqs)
host_dists = utils.mk_count_dists(host_vecs)
flu_vecs = utils.mk_count_vecs(flu_counts, use_seqs)  
flu_dists = utils.mk_count_dists(flu_vecs)

# for flu in flu_dists:
#     print sum(flu_dists[flu])
# sys.exit(0)

for seq, human, chicken in zip(use_seqs, flu_dists['human'],
                               flu_dists['chicken']):
    print seq + '\t' + str(human) + '\t' + str(chicken)

# for flu in flu_counts:
#     for elmSeq in flu_counts[flu]:
#         elm, seq = elmSeq.split(':')
#         if elm == use_elm:
#             print flu, seq, flu_counts[flu][elmSeq]
