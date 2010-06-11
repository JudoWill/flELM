"""Split up genomes for blasting on blades"""

import itertools

def chunks(ls, n):
    """Yield n-sized chunks for ls"""
    for i in xrange(0, len(ls), n):
        yield ls[i:i+n]

GENOMES = ('M_musculus', 'Bos_taurus','Canis_familiaris','D_rerio',
           'Equus_caballus', 'Gallus_gallus', 'H_sapiens',
           'Macaca_mulatta', 'R_norvegicus', 'Sus_scrofa',
           'Taeniopygia_guttata', 'Pan_troglodytes')
pairs = []
for g1,g2 in itertools.combinations(GENOMES,2):
    pairs.append((g1,g2))
    pairs.append((g2,g1))
for count,c in enumerate(chunks(pairs, 33)):
    f = open('chunks/roundup' + str(count), 'w')
    for g1,g2 in c:
        f.write(g1 + '\t' + g2 + '\n')
    f.close()
