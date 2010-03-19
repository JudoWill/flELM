import utils_stats, itertools, utils_distance

# when using entropy
# maybe you should only inlude
# elms common to both species
#measure = '.elm_entropy'
#measure = '.elm_freq'
measure = '.elm_seq_frac'

def loadEntropy(afile):
    e = {}
    with open(afile) as f:
        for line in f:
            [elm, ent] = line.strip().split('\t')
            e[elm] = float(ent)
    return e

def cor(ent1, ent2):
    ls1 = []
    ls2 = []
    for elm in utils_distance.get_elements(ent1, ent2):
        if elm in ent2 and elm in ent1:
            ls1.append(ent1[elm])
            ls2.append(ent2[elm])
        #elif elm in ent1:
        #    ls1.append(ent1[elm])
        #    ls2.append(float(0))
        #else:
        #    ls1.append(float(0))
        #    ls2.append(ent2[elm])
    return utils_stats.pcc(ls1, ls2)

ents = {}
gs = ('H_sapiens', 'Sus_scrofa', 'Gallus_gallus')
flus = ('human', 'swine', 'chicken')
for g  in gs:
    ents[g] = loadEntropy('results/' + g + measure)
for g  in flus:
    ents[g] = loadEntropy('results/flu_' + g + measure)

for g1, g2 in itertools.combinations(gs, 2):
    print g1, g2, cor(ents[g1], ents[g2])

for g in gs:
    for flu in flus:
        print g, flu, cor(ents[g], ents[flu])

for g1, g2 in itertools.combinations(flus, 2):
    print g1, g2, cor(ents[g1], ents[g2])
