""" What are the sequence differences for the
    given ELM for chicken & human H5N1?"""
import sys, utils, os, utils_graph

use_elm = sys.argv[1]

flu_counts = {}
seen_seqs = {}
seen_seqs_ls = []
flus = ('chicken','human')
for flu in flus:
    if flu == 'human':
        flu_elm_file = os.path.join('working/Jul1',
                                    flu + '.H1N1.simpleELMs')
    else:
        flu_elm_file = os.path.join('working/Jul1',
                                    flu + '.H5N1.simpleELMs')
    # flu_elm_file = os.path.join('working/Jun30/',
    #                             flu + '.H5N1.simpleELMs')
    utils.count_flu_sampled(flu, flu_elm_file, flu_counts,
                            seen_seqs, {}, False)
    use_seqs = {}
    for elmSeq in seen_seqs[flu]:
        elm, seq = elmSeq.split(':')
        if elm == use_elm:
            if flu_counts[flu][elmSeq] > 0:
                use_seqs[elmSeq] = True
    seen_seqs_ls.append(use_seqs)

# remove seqs seen less than 10x
#for flu in flu_counts:
#    for elmSeq in 

use_seqs_pre = utils_graph.unionLists(seen_seqs_ls)

counts = utils.count_host_elmSeqs(('Gallus_gallus','H_sapiens'),
                                  False, {},
                                  'working/Jul1/', {use_elm:True},
                                  '.simple')

use_seqs = utils_graph.intersectLists([use_seqs_pre, counts['Gallus_gallus'],
                                   counts['H_sapiens']])
host_vecs = utils.mk_count_vecs(counts, use_seqs)
host_dists = utils.mk_count_dists(host_vecs)
flu_vecs = utils.mk_count_vecs(flu_counts, use_seqs)  
flu_dists = utils.mk_count_dists(flu_vecs)

# for flu in flu_dists:
#     print sum(flu_dists[flu])
# sys.exit(0)

for seq, hhuman, hchicken, hf, cf in zip(use_seqs, host_dists['Gallus_gallus'],
                                         host_dists['H_sapiens'],
                                         flu_dists['human'],
                                         flu_dists['chicken']):
    print seq + '\t' + str(hhuman) + '\t' + str(hchicken) + '\t' + str(hf) + '\t' + str(cf)
for host in host_dists:
    for flu in flu_dists:
        print host, flu, utils.jensen_shannon_dists(host_dists[host],
                                                    flu_dists[flu])
# for flu in flu_counts:
#     for elmSeq in flu_counts[flu]:
#         elm, seq = elmSeq.split(':')
#         if elm == use_elm:
#             print flu, seq, flu_counts[flu][elmSeq]
