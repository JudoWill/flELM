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

#genomes = ('H5N1', 'H9N2')
#species = 'chicken'

genomes = ('H1N1', 'H3N2')
species = 'swine'
conserved = {}
all_conserved = {}
for g in genomes:
    conserved[g] = get_conserved(species + '.' + g + '.elms.90')
    all_conserved[g] = get_all_conserved(species + '.' + g + '.elms.conservation')
for protein in conserved['H1N1']:
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
    utils_plot.distance_heatmap(d, species + '.' + protein + '.png')
