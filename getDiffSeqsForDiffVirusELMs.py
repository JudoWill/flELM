import utils_motif, sys, utils, random, os
from collections import defaultdict

elms = {}
with open('mammal_bird.different.70.test') as f:
    for line in f:
        elms[line.strip().split()[0]] = True

d = {'ELM':True}
swine_H1N1_elms = utils_motif.protein2annotation('results/swine.H1N1.elms', d)
swine_H3N2_elms = utils_motif.protein2annotation('results/swine.H3N2.elms', d)
swine_H1N1_cons = utils_motif.protein2annotation('results/swine.H1N1.elms.90', d)
swine_H3N2_cons = utils_motif.protein2annotation('results/swine.H3N2.elms.90', d)
swine = {'H1N1':swine_H1N1_elms, 'H3N2':swine_H3N2_elms}
swine_conserved = {'H1N1':swine_H1N1_cons, 'H3N2':swine_H3N2_cons}

# human_H1N1_elms = utils_motif.protein2annotation('results/human.H1N1.elms', d)
# human_H3N2_elms = utils_motif.protein2annotation('results/human.H3N2.elms', d)
# human_H5N1_elms = utils_motif.protein2annotation('results/human.H5N1.elms', d)
# human = [human_H1N1_elms, human_H3N2_elms, human_H5N1_elms]

# horse_H3N8_elms = utils_motif.protein2annotation('results/equine.H3N8.elms', d)
# horse = [horse_H3N8_elms]

# chicken_H5N1_elms = utils_motif.protein2annotation('results/chicken.H5N1.elms', d)
# chicken_H9N2_elms = utils_motif.protein2annotation('results/chicken.H9N2.elms', d)
# chicken = [chicken_H5N1_elms, chicken_H9N2_elms]

# duck_H5N1_elms = utils_motif.protein2annotation('results/duck.H5N1.elms', d)
# duck_H9N2_elms = utils_motif.protein2annotation('results/duck.H9N2.elms', d)
# duck = [duck_H5N1_elms, duck_H9N2_elms]

for pig in swine:    
    counts = {}
    for protein in swine[pig]:
        protein_class = protein.split('.')[1] 
        for elm in swine[pig][protein]:
            if elm not in counts:
                counts[elm] = {}
            if protein_class not in counts[elm]:
                counts[elm][protein_class] = {}
            for [st, stp, seq] in swine[pig][protein][elm]:
                if seq not in counts[elm][protein_class]:
                    counts[elm][protein_class][seq] = {}
                counts[elm][protein_class][seq][protein] = True
    for elm in counts:
        line = ''
        for protein in counts[elm]:
            if elm in swine_conserved[pig][protein]:
                for seq in counts[elm][protein]:
                    line += '%s\t%s\t%d\n' % (protein, seq, len(counts[elm][protein][seq]))
        if line:
            tmp_input = 'tmp_i' + str(random.randint(0,100))
            with open(tmp_input, 'w') as f:
                f.write('Protein\tSeq\tCount\n')
                f.write(line)
            r_file = 'tmp_r' + str(random.randint(0,100))
            with open(r_file, 'w') as f:
                f.write('library(ggplot2)\n')
                f.write("d<-read.delim('"
                        + tmp_input + "', header=T, sep='\\t')\n")
                f.write("png('plots/swine_flu_elms/" + pig + '.' + elm + ".png')\n")
                f.write("ggplot(d) + aes(x=Seq,y=Count) + geom_bar(aes(fill=Protein)) + facet_grid(Protein~.) + opts(legend.position='none') + opts(title='" + elm + "')\n")
                f.write('dev.off()\n')
            os.system('R < ' + r_file + ' --no-save')
            os.system('rm ' + r_file + ' ' + tmp_input)
