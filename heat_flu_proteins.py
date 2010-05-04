import utils_plot
from collections import defaultdict

def get_conserved(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            [protein, j1, j2, 
             elm, seq, jelm] = line.strip().split('\t')
            d[protein][elm] = float(1)
    return d

def get_all_conserved(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            [protein, elm, cons] = line.strip().split('\t')
            d[protein][elm] = float(cons)/float(100)
    return d

genomes = ('human', 'swine', 'chicken', 'equine')
conserved = {}
all_conserved = {}
for g in genomes:
    conserved[g] = get_conserved('results/' + g + '.elms.90')
    all_conserved[g] = get_all_conserved('results/' + g + '.elms.conservation')
for protein in conserved['human']:
    d = defaultdict(dict)
    elms = {}
    for g in genomes:
        for elm in conserved[g][protein]:
            elms[elm] = True
    for g in genomes:
        for elm in elms:
            if elm in conserved[g][protein]:
                d[elm][g] = float(1)
            elif elm in all_conserved[g][protein]:
                d[elm][g] = all_conserved[g][protein][elm]
            else:
                d[elm][g] = float(0)
    utils_plot.distance_heatmap(d, protein + '.all.png')
