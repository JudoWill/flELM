ent1 = {}
with open('junk.entropy') as f:
    for line in f:
        [elm, ent] = line.strip().split('\t')
        ent1[elm] = ent
with open('results/H_sapiens.elm_entropy') as f:
    for line in f:
        [elm, ent] = line.strip().split('\t')
        if elm in ent1:
            if float(ent1[elm]) < float(ent):
                print 'yes'
#            print elm + '\t' + ent + '\t' + ent1[elm]
