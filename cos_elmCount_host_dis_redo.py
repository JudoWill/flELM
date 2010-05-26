""" Find cosine distance between ELM counts that occur
    more than expected.  This is different from using
    ELM frequencies.  I'm not using frequencies b/c
    the counts have already been corrected in a sense.
"""
import itertools, sys, os, utils, random
from collections import defaultdict

# 'Macaca_mulatta 'R_norvegicus',  'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus', 'Taeniopygia_guttata', 'Sus_scrofa',
species = ('H_sapiens', 'Macaca_mulatta', 'R_norvegicus',
           'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus',
            'Equus_caballus',
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
freqs = {}
for host in species:
    freqs_single = defaultdict(utils.init_zero)
    with open('results/elmdict_' + host + '.realFraction') as f:
        for line in f:
            (elm, fq) = line.strip().split('\t')
            freqs_single[elm] = float(fq)
    freqs[host] = freqs_single

tmp_input = 'tmp_i' + str(random.randint(0,100))
tmp_r = 'tmp_r' + str(random.randint(0,100))
out_file = 'plots/for_aydin_2/elm_count_dis_redo.png'
with open(tmp_input, 'w') as f:
    f.write('Host1\tHost2\tDistance\n')
    for i in xrange(len(species)):
        for j in xrange(len(species)):
            if i != j:
                host1 = species[i]
                host2 = species[j]
                f.write('%s\t%s\t%.10f\n'
                        % (short_names[host1], 
                           short_names[host2], 
                           utils.getDistance(freqs[host1], 
                                             freqs[host2])))
with open(tmp_r, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("d<-read.delim('"
                + tmp_input + "', header=T, sep='\\t')\n")
        f.write("png('" + out_file + "')\n")
        f.write("ggplot(d,aes(Host1,Host2)) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='white',high='steelblue')\n")
        f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')
#os.system('rm ' + tmp_r + ' ' + tmp_input)

    
