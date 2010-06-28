""" What is the Jensen-shannon distance
    between 2 flu groups?
"""
import utils, sys, utils_graph

flu_group_1_file = sys.argv[1]
flu_group_2_file = sys.argv[2]

flu_counts = {}
seen_elmSeqs = {}
seen_seqs_ls = []
for name, file in (('g1', flu_group_1_file),
                   ('g2', flu_group_2_file)):
    utils.count_flu_sampled(name, file,
                            flu_counts,
                            seen_elmSeqs, {}, False)
    seen_seqs_ls.append(seen_elmSeqs[name])
use_elmSeqs = utils_graph.unionLists(seen_seqs_ls)
flu_vecs = utils.mk_count_vecs(flu_counts, use_elmSeqs)
flu_dists = utils.mk_count_dists(flu_vecs)
js_dis = utils.jensen_shannon_dists(flu_dists['g1'],
                                    flu_dists['g2'])
print js_dis
    

