import itertools, sys, os, utils, random
from collections import defaultdict

suffix = sys.argv[1]

# 'Macaca_mulatta 'R_norvegicus',  'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus', 'Taeniopygia_guttata'
species = ('H_sapiens', 'Macaca_mulatta', 'R_norvegicus',
           'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus',
           'Sus_scrofa', 'Equus_caballus',
           'Gallus_gallus', 'Taeniopygia_guttata')
short_names = {'H_sapiens':'Us',
               'Macaca_mulatta':'Chimp',
               'M_musculus':'Mice',
               'R_norvegicus':'Rat',
               'Sus_scrofa':'Pig',
               'Equus_caballus':'Hrse',
               'Canis_familiaris':'Cat',
               'Bos_taurus':'Cow',
               'Gallus_gallus':'Chick',
               'D_rerio':'Fish',
               'Taeniopygia_guttata':'Fnch'}

#'.elm_aa_freq'
freqs = defaultdict(dict)
elms = {}
for host in species:
    with open('results/' + host + suffix + '.elm_aa_freq') as f:
        for line in f:
            (elm, fq) = line.strip().split('\t')
            elms[elm] = True
            freqs[host][elm] = float(fq)

tmp_input = 'tmp_i' + str(random.randint(0,100))
#tmp_input = 'plots/for_aydin_2/elm_freq_dis' + suffix + '.tab'
tmp_r = 'tmp_r' + str(random.randint(0,100))
tmp_labels = 'labels' + str(random.randint(0,100))
out_file = 'plots/for_aydin_2/cos_elmFreq_host_dis.dendrogram' + suffix + '.png'
species_lines = {}
for s in species:
    species_lines[s] = ''

for elm in elms:
    for s in species:
        if elm in freqs[s]:
            species_lines[s] += str(freqs[s][elm]) + '\t'
        else:
            species_lines[s] += '0\t'
with open(tmp_input, 'w') as f:    
    for s in species:
        f.write(species_lines[s].strip('\t') + '\n')
with open(tmp_labels, 'w') as f:
    f.write('\t'.join(species) + '\n')
with open(tmp_r, 'w') as f:
    f.write("d<-read.delim('"
            + tmp_input
            + "',header=FALSE,sep='\t')\n")
    f.write("labels.d<-read.delim('"
            + tmp_labels
            + "',header=FALSE,sep='\t')\n")
    f.write('labels<-as.matrix(labels.d)\n')
    f.write('dist.r<-dist(d)\n')
    f.write("h<-hclust(dist.r,method='complete')\n")
    f.write("png('" + out_file + "')\n")
    f.write('plot(h,hang=-1,labels=labels[1,])\n')
    f.write('dev.off()\n')
