""" Find ELMs that are not conserve by change.
    Enter subtype & cutoff for # trial an ELM for a 
    protein can be found by chance.
"""
import sys, utils_motif, utils_graph
sys.path.append('../')
import utils
from collections import defaultdict
elm_file = sys.argv[1]
cutoff = int(sys.argv[2])

real = utils_motif.protein2annotation('../results/' + elm_file,
                                      {'ELM':True})

protein2elms = {}
for x in xrange(10):
    protein2annotation = utils_motif.protein2annotation(str(x) + '/'
                                                        + elm_file,
                                                        {'ELM':True})
    for protein in protein2annotation:
        for elm in protein2annotation[protein]:
            if not protein in protein2elms:
                protein2elms[protein] = {}
            if not elm in protein2elms[protein]:
                protein2elms[protein][elm] = 0
            protein2elms[protein][elm] += 1

for protein in real:
    if protein in protein2elms:
        for elm in real[protein]:
            if elm in protein2elms[protein]:
                if protein2elms[protein][elm] < cutoff:
                    print protein + '\t' + elm
                else:
                    print "FAIL\t" + protein + '\t' + elm
            else:
                print protein + '\t' + elm
    else:
        print protein + '\tno protein control hits'
