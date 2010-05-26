""" For conserved virus ELMs, find the cosine distance
    between sequence fraction vectors of virus and host.

    Updated to use HIV and HCV conserved ELMs.
"""
import utils, local_settings, os, utils_motif, utils_graph, random, os, sys, numpy
from collections import defaultdict

suffix = sys.argv[1]

def getConservedELMs(virus, subtypes):
    ls = [utils_motif.annotation2protein(os.path.join(local_settings.RESULTSDIR,
                                                      virus + '.' + subtype 
                                                      + '.elms.70.controled'),
                                         {'ELM':True}) 
          for subtype in subtypes[virus]]
    return utils_graph.intersectLists(ls)

hosts = {'H_sapiens':'M', 'Sus_scrofa':'P', 'Gallus_gallus':'C', 
         'Equus_caballus':'H', 'Taeniopygia_guttata':'F'}
viruses = {'human':'m', 'swine':'p', 'equine':'h',
           'duck':'d', 'chicken':'c'}#, 'HIV':'a', 'HCV':'l'}
subtypes = {'human':('H1N1', 'H3N2', 'H5N1'),
            'swine':('H1N1', 'H3N2'),
            'equine':('H3N8',),
            'chicken':('H5N1', 'H9N2'),
            'duck':('H5N1', 'H9N2')}
            #'HIV':('all',),
            #'HCV':('all',)}

virus2conservedELMs = {}
all_elms = {}
for virus in viruses:
    virus2conservedELMs[virus] = getConservedELMs(virus, subtypes)
    for elm in virus2conservedELMs[virus]:
        all_elms[elm] = True

# load ELM seq fractions
host2elmFreqs = {}
virus2elmFreqs = {}
use_seqs = {}
for host in hosts:
    host2elmFreqs[host] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                'elmdict_' + host + suffix),
                                                   float(0))
    for elm in host2elmFreqs[host]:
        if elm not in use_seqs:
            use_seqs[elm] = {}
        for seq in host2elmFreqs[host][elm]:
            if seq not in use_seqs[elm]:
                use_seqs[elm][seq] = 0
            use_seqs[elm][seq] += 1

for elm in use_seqs:
    rm_ls = []
    for seq in use_seqs[elm]:
        if use_seqs[elm][seq] != len(hosts.keys()):
            rm_ls.append(seq)
    for seq in rm_ls:
        del use_seqs[elm][seq]

means = {}
vars = {}
for elm in all_elms:
    means[elm] = {}
    vars[elm] = {}
    seq_vals = defaultdict(list)
    for host in hosts:
        for seq in host2elmFreqs[host][elm]:
            seq_vals[seq].append(host2elmFreqs[host][elm][seq])
    for seq in seq_vals:
        missing = len(hosts.keys())-len(seq_vals[seq])
        for m in xrange(missing):
            seq_vals[seq].append(0)
    
        means[elm][seq] = numpy.average(numpy.array(seq_vals[seq]))
        #inverse variance
        vars[elm][seq] = numpy.var(numpy.array(seq_vals[seq]))

#tmp_input = 'tmp_input' + str(random.randint(0,100))
tmp_input = 'plots/for_aydin/cos_host_virus' + suffix + '.tab'
with open(tmp_input, 'w') as f:
    f.write('Virus_Host\tELM\tDistance\n')
    for virus in viruses:
        virus2elmFreqs[virus] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                      'flu_elmdict_' + virus), float(0))
        for elm in virus2conservedELMs[virus]:
            if 'FAIL' not in elm:
                for host in hosts:
                    if elm in use_seqs:
                        dis = utils.klDistance(virus2elmFreqs[virus][elm], 
                                               host2elmFreqs[host][elm],
                                               use_seqs[elm])
                    else:
                        dis = numpy.NaN
                    f.write('%s\t%s\t%.10f\n'
                            % (viruses[virus] + hosts[host], elm, dis))
out_file = 'plots/for_aydin/cos_dis_heatmap' + suffix + '.png'
tmp_r = 'tmp_r' + str(random.randint(0,100))
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d,aes(Virus_Host,ELM)) + opts(axis.text.y = theme_blank()) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='red',high='steelblue')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')

