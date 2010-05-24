import utils_motif, sys, utils, random, os
from collections import defaultdict

elms = {}
with open('mammal_bird.different.70.test') as f:
    for line in f:
        elms[line.strip().split()[0]] = True

d = {'ELM':True}
swine_H1N1_elms = utils_motif.protein2annotation('results/swine.H1N1.elms', d)
swine_H3N2_elms = utils_motif.protein2annotation('results/swine.H3N2.elms', d)
swine_H1N1_cons = utils_motif.protein2annotation('results/swine.H1N1.elms.70', d)
swine_H3N2_cons = utils_motif.protein2annotation('results/swine.H3N2.elms.70', d)
hits = {'pig.H1N1':swine_H1N1_elms, 'pig.H3N2':swine_H3N2_elms}
hits_conserved = {'pig.H1N1':swine_H1N1_cons, 'pig.H3N2':swine_H3N2_cons}

human_H1N1_elms = utils_motif.protein2annotation('results/human.H1N1.elms', d)
human_H3N2_elms = utils_motif.protein2annotation('results/human.H3N2.elms', d)
human_H5N1_elms = utils_motif.protein2annotation('results/human.H5N1.elms', d)
human_H1N1_cons = utils_motif.protein2annotation('results/human.H1N1.elms.70', d)
human_H3N2_cons = utils_motif.protein2annotation('results/human.H3N2.elms.70', d)
human_H5N1_cons = utils_motif.protein2annotation('results/human.H5N1.elms.70', d)
hits['human.H1N1'] = human_H1N1_elms
hits['human.H3N2'] = human_H3N2_elms
hits['human.H5N1'] = human_H5N1_elms
hits_conserved['human.H1N1'] = human_H1N1_cons
hits_conserved['human.H3N2'] = human_H3N2_cons
hits_conserved['human.H5N1'] = human_H5N1_cons

# horse_H3N8_elms = utils_motif.protein2annotation('results/equine.H3N8.elms', d)
# horse = [horse_H3N8_elms]

chicken_H5N1_elms = utils_motif.protein2annotation('results/chicken.H5N1.elms', d)
chicken_H9N2_elms = utils_motif.protein2annotation('results/chicken.H9N2.elms', d)
chicken_H5N1_cons = utils_motif.protein2annotation('results/chicken.H5N1.elms', d)
chicken_H9N2_cons = utils_motif.protein2annotation('results/chicken.H9N2.elms', d)
hits['chicken.H5N1'] = chicken_H5N1_elms
hits['chicken.H9N2'] = chicken_H9N2_elms
hits_conserved['chicken.H5N1'] = chicken_H5N1_cons
hits_conserved['chicken.H9N2'] = chicken_H9N2_cons

# duck_H5N1_elms = utils_motif.protein2annotation('results/duck.H5N1.elms', d)
# duck_H9N2_elms = utils_motif.protein2annotation('results/duck.H9N2.elms', d)
# duck = [duck_H5N1_elms, duck_H9N2_elms]

use_elms = {}
for species in hits_conserved:
    for protein in hits_conserved[species]:
        for elm in hits_conserved[species][protein]:
            use_elms[elm] = True

counts = {}
for pig in hits:    
    for protein in hits[pig]:
        protein_class = protein.split('.')[1] 
        for elm in hits[pig][protein]:
            if elm not in counts:
                counts[elm] = {}
            if pig not in counts[elm]:
                counts[elm][pig] = {}
            if protein_class not in counts[elm][pig]:
                counts[elm][pig][protein_class] = {}
            for [st, stp, seq] in hits[pig][protein][elm]:
                if seq not in counts[elm][pig][protein_class]:
                    counts[elm][pig][protein_class][seq] = {}
                counts[elm][pig][protein_class][seq][protein] = True
for elm in elms:
    line = ''
    for pig in counts[elm]:
        for protein in counts[elm][pig]:
            if protein in hits_conserved[pig]:
                if elm in hits_conserved[pig][protein]:
                    for seq in counts[elm][pig][protein]:
                        line += '%s\t%s\t%s\t%d\n' % (pig.split('.')[0], protein + '.' + pig.split('.')[1], 
                                                      seq, len(counts[elm][pig][protein][seq]))
    if line:
        tmp_input = 'tmp_i' + str(random.randint(0,100))
        with open(tmp_input, 'w') as f:
            f.write('Species\tProtein\tSeq\tCount\n')
            f.write(line)
        r_file = 'tmp_r' + str(random.randint(0,100))
        with open(r_file, 'w') as f:
            f.write('library(ggplot2)\n')
            f.write("d<-read.delim('"
                    + tmp_input + "', header=T, sep='\\t')\n")
            f.write("png('plots/swine_chick_flu_elms/" + elm + ".png')\n")
            line = "ggplot(d) + aes(x=Seq,y=Count) + geom_bar(aes(fill=Species)) + facet_grid(Protein~.) + opts(title='" + elm + "')\n"
            f.write(line)
            f.write('dev.off()\n')
        os.system('R < ' + r_file + ' --no-save')
        os.system('rm ' + r_file + ' ' + tmp_input)
