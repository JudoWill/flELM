import sys, os, random
from collections import defaultdict

suffix = sys.argv[1]

species = ('H_sapiens', 'Macaca_mulatta', 'R_norvegicus',
           'Canis_familiaris', 'Bos_taurus', 'D_rerio', 'M_musculus',
           'Sus_scrofa', 'Equus_caballus',
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
out = 'plots/for_aydin_2/elm_seq_hits_hist_speciesSeq' + suffix + '.png'

host2seq = {}
elm2seqs = defaultdict(dict)
for host in species:
    seqs = defaultdict(dict)
    with open('results/elmdict_' + host + suffix) as finput:
        for line in finput:
            (elm, seq, count, fq) = line.strip().split('\t')
            seqs[elm][seq] = True
            elm2seqs[elm][seq] = True
    host2seq[host] = seqs

# plot fraction of sequences counting species:seq as a unique item     
with open(r_input, 'w') as f:
    f.write('Species\tELM\tSeqs\n')
    for elm in elm2seqs:
        total = 0
        for host in species:
            total += len(host2seq[host][elm].keys())
        for host in species:
            frac = float(len(host2seq[host][elm].keys()))/float(total)
            f.write('%s\t%s\t%.10f\n' %
                    (short_names[host],
                     elm, frac))
            
with open(r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + r_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out + "')\n")
    f.write("ggplot(d) + aes(x=ELM,y=Seqs) + geom_bar(aes(fill=ELM)) + facet_grid(Species~.) + opts(legend.position='none') + opts(title='ELM species:seq')\n")
    f.write('dev.off()\n')
os.system('R < ' + r + ' --no-save')

# plot fraction of sequences counting seq as a unique item  
out = 'plots/for_aydin_2/elm_seq_hits_hist_seq' + suffix + '.png'     
with open(r_input, 'w') as f:
    f.write('Species\tELM\tSeqs\n')
    for elm in elm2seqs:
        total = float(len(elm2seqs[elm].keys()))
        for host in species:
            frac = float(len(host2seq[host][elm].keys()))/total
            f.write('%s\t%s\t%.10f\n' %
                    (short_names[host],
                     elm, frac))
            
with open(r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + r_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out + "')\n")
    f.write("ggplot(d) + aes(x=ELM,y=Seqs) + geom_bar(aes(fill=ELM)) + facet_grid(Species~.) + opts(legend.position='none') + opts(title='ELM seq')\n")
    f.write('dev.off()\n')
os.system('R < ' + r + ' --no-save')
os.system('rm ' + r + ' ' + r_input)
