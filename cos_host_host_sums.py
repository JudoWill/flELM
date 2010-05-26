""" For conserved virus ELMs, find the cosine distance
    between sequence fraction vectors for host pairs.
"""
import utils, local_settings, os, utils_motif, utils_graph, random, os, itertools, utils, sys, numpy
from collections import defaultdict

suffix = sys.argv[1]

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

# load ELM seq fractions
host2elmFreqs = {}
for host in species:
    host2elmFreqs[host] = utils.get_seq2count_dict(os.path.join(local_settings.RESULTSDIR,
                                                                'elmdict_' + host + suffix),
                                                   float(0))

tmp_input = 'plots/for_aydin_2/cos_host_host' + suffix + '.tab'
with open(tmp_input, 'w') as f:
    f.write('Host1\tHost2\tDistance\n')
    for i in xrange(len(species)):
        for j in xrange(len(species)):
            if i != j:
                host1 = species[i]
                host2 = species[j]
                sum = float(0)
                elms_compared = 0
                for elm in utils_graph.unionLists([host2elmFreqs[host1],
                                                   host2elmFreqs[host2]]):
                    for h in (host1, host2):
                        if elm not in host2elmFreqs[h]:
                            host2elmFreqs[h][elm] = {}
                    if len(host2elmFreqs[host1][elm].keys()) != 0 and len(host2elmFreqs[host2][elm].keys()) != 0:
                        sum += utils.getDistance(host2elmFreqs[host1][elm],
                                                 host2elmFreqs[host2][elm])
                        elms_compared += 1
                f.write('%s\t%s\t%.10f\n'
                        % (short_names[host1], 
                           short_names[host2], 
                           float(sum)/float(elms_compared)))
           
out_file = 'plots/for_aydin_2/cos_dis_sum_host' + suffix + '.png'
tmp_r = 'tmp_r' + str(random.randint(0,100))
with open(tmp_r, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + tmp_input + "', header=T, sep='\\t')\n")
    f.write("png('" + out_file + "')\n")
    f.write("ggplot(d,aes(Host1,Host2)) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='green',high='steelblue')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')
os.system('rm ' + tmp_r)


