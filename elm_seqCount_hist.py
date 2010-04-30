""" Make a histogram that counts the
    # of strings an ELM can match.
"""
import global_settings, local_settings, os, random
from collections import defaultdict

def get_seqs(afile, d):
    with open(afile) as f:
        for line in f:
            elm, seq, count, frac = line.strip().split('\t')
            d[elm][seq] = True

counts = defaultdict(dict)
for genome in global_settings.GENOMES:
    get_seqs(os.path.join(local_settings.RESULTSDIR,
                          'elmdict_' + genome + '.txt'),
             counts)
ls = []
for elm in counts:
    ls.append( (len(counts[elm].keys()),elm) )
ls.sort()
ls.reverse()
tmp_input = 'tmp_input' + str(random.randint(0,100))
with open(tmp_input, 'w') as f:
    f.write('ELM\tSeqCount\n')
    for count,elm in ls:
        f.write('%s\t%d\n' 
                % (elm, count))

tmp_r = 'tmp_r' + str(random.randint(0,100))
out_file = 'plots/for_aydin/elm_seqCount_hist.png'
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d) + aes(x=ELM,y=SeqCount) + geom_bar() + opts(legend.position='none', axis.text.x = theme_blank()) + opts(title='Sequence Counts Per ELM')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')
os.system('rm ' + tmp_r + ' ' + tmp_input)
