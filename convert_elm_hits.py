""" I need to add HIV and HCV to the project,
    but I must first convert the ELMs hits to
    the frequencies used for this project.
"""
import utils_motif, sys

def printELMusage(elm, seq2count):
    total = 0
    for seq in seq2count:
        total += seq2count[seq]
    for seq in seq2count:
        print('%s\t%s\t%d\t%.10f' %
              (elm, seq, seq2count[seq],
               float(seq2count[seq]/float(total))))

protein2elm = utils_motif.protein2annotation(sys.argv[1],
                                             {'ELM':True})
elm2seq2count = {}
for p in protein2elm:
    for elm in protein2elm[p]:
        if elm not in elm2seq2count:
            elm2seq2count[elm] = {}
        for [st, stp, seq] in protein2elm[p][elm]:
            if seq not in elm2seq2count[elm]:
                elm2seq2count[elm][seq] = 0
            elm2seq2count[elm][seq] += 1
for elm in elm2seq2count:
    printELMusage(elm, elm2seq2count[elm])
