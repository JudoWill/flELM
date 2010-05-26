""" For conserved virus ELMs, find the cosine distance
    between sequence fraction vectors for host pairs.
"""
import utils, local_settings, os, utils_motif, utils_graph, random, os, itertools, utils, sys, numpy
from collections import defaultdict

suffix = sys.argv[1]

def getConservedELMs(virus, subtypes):
    ls = [utils_motif.annotation2protein(os.path.join(local_settings.RESULTSDIR,
                                                      virus + '.' + subtype + '.elms.70'),
                                         {'ELM':True}) 
          for subtype in subtypes[virus]]
    return utils_graph.intersectLists(ls)

hosts = {'H_sapiens':'A', 'Sus_scrofa':'P', 'Gallus_gallus':'C', 
         'Equus_caballus':'H', 'Taeniopygia_guttata':'F',
         'M_musculus':'M','R_norvegicus':'R','D_rerio':'Z',
         'Bos_taurus':'B','Macaca_mulatta':'L','Canis_familiaris':'G'}
viruses = {'human':'m', 'swine':'p', 'equine':'h',
           'duck':'d', 'chicken':'c'}
subtypes = {'human':('H1N1', 'H3N2', 'H5N1'),
            'swine':('H1N1', 'H3N2'),
            'equine':('H3N8',),
            'chicken':('H5N1', 'H9N2'),
            'duck':('H5N1', 'H9N2'),
            'HIV':('all',),
            'HCV':('all',)}

virus2conservedELMs = {}
all_elms = {}
for virus in viruses:
    virus2conservedELMs[virus] = getConservedELMs(virus, subtypes)
    for elm in virus2conservedELMs[virus]:
        all_elms[elm] = True
# load ELM seq fractions
host2elmFreqs = {}
virus2elmFreqs = {}
for host in hosts:
    host2elmFreqs[host] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                'elmdict_' + host + suffix),
                                                   float(0))
# means = {}
# vars = {}
# for elm in all_elms:
#     means[elm] = {}
#     vars[elm] = {}
#     seq_vals = defaultdict(list)
#     for host in hosts:
#         for seq in host2elmFreqs[host][elm]:
#             seq_vals[seq].append(host2elmFreqs[host][elm][seq])
#     for seq in seq_vals:
#         missing = len(hosts.keys())-len(seq_vals[seq])
#         for m in xrange(missing):
#             seq_vals[seq].append(0)
    
#         means[elm][seq] = numpy.average(numpy.array(seq_vals[seq]))
#         #inverse variance
#         vars[elm][seq] = numpy.var(numpy.array(seq_vals[seq]))
#print means
#tmp_input = 'tmp_input' + str(random.randint(0,100))
tmp_input = 'plots/for_aydin/cos_host_host' + suffix + '.tab'
with open(tmp_input, 'w') as f:
    f.write('Host_Host\tELM\tDistance\n')
    for elm in all_elms:
        counter = 0
        for host1,host2 in (('H_sapiens', 'Macaca_mulatta'),
                            ('H_sapiens', 'M_musculus'),
                            ('H_sapiens', 'R_norvegicus'),
                            ('H_sapiens', 'Sus_scrofa'),
                            ('H_sapiens', 'Equus_caballus'),
                            ('H_sapiens', 'Canis_familiaris'),
                            ('H_sapiens', 'Bos_taurus'),
                            ('H_sapiens', 'Gallus_gallus'),
                            ('H_sapiens', 'Taeniopygia_guttata'),
                            ('H_sapiens', 'D_rerio')
                            ):#itertools.combinations(hosts.keys(),2):
            use_dict_1 = {}
            use_dict_2 = {}
            if elm in host2elmFreqs[host1]:
                use_dict_1 = host2elmFreqs[host1][elm]
            if elm in host2elmFreqs[host2]:
                use_dict_2 = host2elmFreqs[host2][elm]
            if len(use_dict_1.keys()) == 0 and len(use_dict_2.keys()) == 0:
                pass
            #    f.write('%s\t%s\t%.10f\n'
            #            % (str(10), elm, float(0)))
            else:
                dis = utils.getDistance(use_dict_2, 
                                        use_dict_1)
                f.write('%s\t%s\t%.10f\n'
                        % (str(counter), elm, dis))
            counter += 1

out_file = 'plots/for_aydin_2/cos_dis_host' + suffix + '.png'
tmp_r = 'tmp_r' + str(random.randint(0,100))
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d,aes(Host_Host,ELM)) + opts(axis.text.x = theme_blank()) + opts(axis.text.y = theme_blank()) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='green',high='steelblue')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')


