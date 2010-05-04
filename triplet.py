import utils_motif

def get(afile):
    d = {}
    with open(afile) as f:
        for line in f:
            a_class, elm, freq = line.strip().split('\t')
            d[elm] = freq
    return d

human = get('H_sapiens.redo.ratio')
swine = get('Sus_scrofa.redo.ratio')
chicken = get('Gallus_gallus.redo.ratio')
finch = get('Taeniopygia_guttata.redo.ratio')

all_elms = {}
for d in (human, swine, chicken):
    for elm in d:
        all_elms[elm] = True
print 'elm\thuman\tswine\tchicken\tfinch'
for elm in all_elms:
    ls = []
    for d in (human, swine, chicken, finch):
        if elm in d:
            ls.append(d[elm])
        else:
            ls.append('0')
    #if not (ls[0] == float(1) and ls[1] == float(1) and ls[2] == float(1)):
    print elm + '\t' + '\t'.join(ls)
