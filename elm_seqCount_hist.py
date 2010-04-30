""" Make a histogram that counts the
    # of strings an ELM can match.

    Updated to use log10(seqCount) and entropy.
"""
import global_settings, local_settings, os, random, math, utils
from collections import defaultdict

def get_seqs(afile, d):
    with open(afile) as f:
        for line in f:
            elm, seq, count, frac = line.strip().split('\t')
            if elm not in d:
                d[elm] = {}
            if seq not in d[elm]:
                d[elm][seq] = 0
            d[elm][seq] += int(count)

def get_entropy(counts):
    elm2entropy = {}
    for elm in counts:
        ls = []
        total = 0
        for seq in counts[elm]:
            c = counts[elm][seq]
            ls.append(c)
            total += c
        entropy = float(0)
        for item in ls:
            prob = float(item)/float(total)
            entropy =- prob * math.log(prob, 2)
        elm2entropy[elm] = entropy
    return elm2entropy

counts = {}
for genome in global_settings.GENOMES:
    get_seqs(os.path.join(local_settings.RESULTSDIR,
                          'elmdict_' + genome + '.txt'),
             counts)
entropy = get_entropy(counts)
tmp_input = 'tmp_input' + str(random.randint(0,100))
with open(tmp_input, 'w') as f:
    f.write('ELM\tMeasure\tVal\n')
    for elm in counts:
        f.write('%s\tLog10SeqCount\t%.10f\n' 
                % (elm, math.log(len(counts[elm].keys()),10)))
    for elm in entropy:
        f.write('%s\tEntropy\t%.10f\n' 
                % (elm, entropy[elm]))

tmp_r = 'tmp_r' + str(random.randint(0,100))
out_file = 'plots/for_aydin/elm_seqCount_hist.png'
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d) + aes(x=ELM,y=Val) + geom_bar() + facet_grid(Measure~., scales='free') + opts(legend.position='none', axis.text.x = theme_blank()) + opts(title='Metrics Per ELM')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')
#os.system('rm ' + tmp_r + ' ' + tmp_input)
