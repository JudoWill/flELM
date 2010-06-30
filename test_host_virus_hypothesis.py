""" For the input ELMs, take the JS
    divergence for chicken/human H5N1
    to human & chicken. For which ELMs
    does the hypothesis holds."""
import sys, utils, os, utils_graph
from collections import defaultdict

elm_file = sys.argv[1]

working_elms = utils_graph.getNodes(elm_file)

flu_counts = {}
seen_seqs = {}
seen_seqs_ls = []
elm2seqs = defaultdict(dict)
flus = ('human', 'chicken')
for flu in flus:
    # flu_elm_file = os.path.join('results',
    #                             flu + '.H5N1.elms')
    flu_elm_file = os.path.join('working/Jun30/',
                                flu + '.H5N1.simpleELMs')
    utils.count_flu_sampled(flu, flu_elm_file, flu_counts,
                            seen_seqs, {}, False)
    for elmseq in seen_seqs[flu]:
        elm, seq = elmseq.split(':')
        elm2seqs[elm][elmseq] = True

counts = utils.count_host_elmSeqs(('Gallus_gallus','H_sapiens'),
                                  False, {},
                                  'working/Jun30/', working_elms,
                                  '.simple')

for elm in working_elms:
    use_seqs = elm2seqs[elm]
    host_vecs = utils.mk_count_vecs(counts, use_seqs)
    host_dists = utils.mk_count_dists(host_vecs)
    flu_vecs = utils.mk_count_vecs(flu_counts, use_seqs)  
    flu_dists = utils.mk_count_dists(flu_vecs)
    
    flu = flu_dists['human']
    human_score_H =  utils.jensen_shannon_dists(host_dists['H_sapiens'],
                                              flu)
    chicken_score_H =  utils.jensen_shannon_dists(host_dists['Gallus_gallus'],
                                                flu)
    
    flu = flu_dists['chicken']
    human_score_C =  utils.jensen_shannon_dists(host_dists['H_sapiens'],
                                              flu)
    chicken_score_C =  utils.jensen_shannon_dists(host_dists['Gallus_gallus'],
                                                flu)

    if human_score_C > chicken_score_C and human_score_H < chicken_score_H:
        print elm
