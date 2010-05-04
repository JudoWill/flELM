""" Find ELMs conserved on some protein
    for X % of all subtypes w/ at least
    50 sequences.

    Made for HIV and HCV data.
"""
import utils_motif, sys

def get_annotations(protein2elm):
    type_counts = {}
    elm_counts = {}
    for p in protein2elm:
        h = p.split('.')[-1]
        t = p.split('.')[0]
        if h not in type_counts:
            type_counts[h] = {}
            elm_counts[h] = {}
        if t not in type_counts[h]:
            type_counts[h][t] = {}
            elm_counts[h][t] = {}
        type_counts[h][t][p] = True
        for elm in protein2elm[p]:
            if elm not in elm_counts[h][t]:
                elm_counts[h][t][elm] = {}
            for start, stop, seq in protein2elm[p][elm]:
                elm_counts[h][t][elm][p] = True
    elm2type2percent = {}
    for h in type_counts:
        for t in type_counts[h]:
            seqs = len(type_counts[h][t].keys())
            if seqs >= 50:
                if h not in elm2type2percent:
                    elm2type2percent[h] = {}
                for elm in elm_counts[h][t]:
                    if elm not in elm2type2percent[h]:
                        elm2type2percent[h][elm] = {}
                    elm2type2percent[h][elm][t] = float(len(elm_counts[h][t][elm].keys()))/float(seqs)
    return elm2type2percent

protein2elm = utils_motif.protein2annotation(sys.argv[1],
                                             {'ELM':True})
cut = float(sys.argv[2])/float(100)

protein2elm2type2percent = get_annotations(protein2elm)
for p in protein2elm2type2percent:
    for elm in protein2elm2type2percent[p]:
        all_pass = True
        for t in protein2elm2type2percent[p][elm]:
            if protein2elm2type2percent[p][elm][t] >= cut:
                print('%s\t0\t0\t%s\tNA\tELM' %
                      (p, elm))
        
