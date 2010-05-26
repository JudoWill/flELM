import sys, os, random
from collections import defaultdict
#'Sus_scrofa',
species = ('H_sapiens', 'Macaca_mulatta', 'R_norvegicus',
           'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus',
            'Equus_caballus',
           'Gallus_gallus', 'Taeniopygia_guttata')
short_names = {'H_sapiens':'a_Us',
               'Macaca_mulatta':'a_Chimp',
               'M_musculus':'b_Mice',
               'R_norvegicus':'b_Rat',
               'Sus_scrofa':'d_Pig',
               'Equus_caballus':'d_Hrse',
               'Canis_familiaris':'d_Cat',
               'Bos_taurus':'d_Cow',
               'Gallus_gallus':'c_Chick',
               'D_rerio':'f_Fish',
               'Taeniopygia_guttata':'c_Fnch'}

r_input = str(random.randint(0,100))
r = 'r' + str(random.randint(0,100))
out = 'plots/for_aydin_2/realFraction_hist.png'

host2elms = {}
for host in species:
    elms = {}
    with open('results/elmdict_' + host + '.realFraction') as finput:
        for line in finput:
            (elm, fq) = line.strip().split('\t')
            elms[elm] = float(fq)
    host2elms[host] = elms

# plot fraction of sequences counting species:seq as a unique item     
with open(r_input, 'w') as f:
    f.write('Species\tELM\tRealFraction\n')
    for host in host2elms:
        for elm in host2elms[host]:
            f.write('%s\t%s\t%.10f\n' %
                    (short_names[host],
                     elm, host2elms[host][elm]))
            
with open(r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + r_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out + "')\n")
    f.write("ggplot(d) + aes(x=ELM,y=RealFraction) + geom_bar(aes(fill=ELM)) + facet_grid(Species~.) + opts(legend.position='none') + opts(title='ELM real hits')\n")
    f.write('dev.off()\n')
os.system('R < ' + r + ' --no-save')

