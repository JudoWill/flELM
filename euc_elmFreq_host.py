import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

def get_use_ELMs_permissive(host2elm):
    elms = {}
    for host in host2elm:
        for elm in host2elm[host]:
            elms[elm] = True
    return elms

def get_use_ELMs_restrictive(host2elm):
    hosts = len(host2elm.keys())
    elms = defaultdict(utils.init_zero)
    for host in host2elm:
        for elm in host2elm[host]:
            elms[elm] += 1
    use_elms = []
    for elm in elms:
        if elms[elm] == hosts:
            use_elms.append(elm)
    return use_elms

def norm_freqs(elm2freqs, hosts):
    new_data = {}
    for elm in elm2freqs:
        s = sum(elm2freqs[elm])
        new_data[elm] = [x/s for x in elm2freqs[elm]]
        l = len(new_data[elm])
        while l < len(hosts):
            new_data[elm].append(float(0))
            l += 1
    return new_data 

def threshold_vars(normed_elm_freqs, thresh):
    """variance in numpy corresponds w/ wikipedia"""
    keep_elms = {}
    for elm in normed_elm_freqs:
        v = numpy.var(normed_elm_freqs[elm]) 
        if v < thresh:
            keep_elms[elm] = v
    return keep_elms

def threshold_vars_down(normed_elm_freqs, thresh):
    """variance in numpy corresponds w/ wikipedia"""
    keep_elms = {}
    for elm in normed_elm_freqs:
        v = numpy.var(normed_elm_freqs[elm]) 
        if v > thresh:
            keep_elms[elm] = v
    return keep_elms
    
suffix = '.redoWNeg'
results_dir = sys.argv[1]
out_file = sys.argv[2]

species = global_settings.MAMMALS3
short_names = global_settings.ALIASES

#'.elm_aa_freq'
freqs = defaultdict(dict)
elm2freqs = defaultdict(list)
for host in species:
    with open(results_dir + host + suffix + '.elm_aa_freq') as f:
        for line in f:
            (elm, fq) = line.strip().split('\t')
            freqs[host][elm] = float(fq)
            elm2freqs[elm].append(float(fq))

tmp_input = 'tmp_i' + str(random.randint(0,100))
tmp_r = 'tmp_r' + str(random.randint(0,100))
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
# with open(tmp_r, 'w') as f:
#         f.write('library(ggplot2)\n')
#         f.write("d<-read.delim('"
#                 + tmp_input + "', header=T, sep='\\t')\n")
#         f.write("png('" + out_file + "')\n")
#         f.write("ggplot(d,aes(Host1,Host2)) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='white',high='steelblue')\n")
#         f.write('dev.off()\n')
#         f.write('q()\n')
# os.system('R < ' + tmp_r + ' --no-save')
#os.system('rm ' + tmp_r + ' ' + tmp_input)

tmp_labels = 'labels' + str(random.randint(0,100))
species_lines = {}
for s in species:
    species_lines[s] = ''

# don't include ELMs that are not present in all species
nfreqs = norm_freqs(elm2freqs, species)
t = threshold_vars(nfreqs, float(10000))
for elm in t:
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
    f.write("source('funcs.R')\n")
    f.write("library('MASS')\n")
    f.write("d<-read.delim('"
            + tmp_input
            + "',header=FALSE,sep='\\t')\n")
    f.write("labels.d<-read.delim('"
            + tmp_labels
            + "',header=FALSE,sep='\\t')\n")
    f.write('labels<-as.matrix(labels.d)\n')
    #f.write("dist.r<-cos_dist(d)\n")
    f.write('dn<-norm(d)\n')
    #f.write("write.matrix(dn,file='plots/for_aydin_2/norm.tab',sep='\t')\n")
    f.write("dist.r<-dist(d,method='euclidean')\n")
    f.write("h<-hclust(dist.r,method='average')\n")
    f.write("png('" + out_file + "')\n")
    f.write("plot(h,hang=-1,labels=labels[1,],main='Host ELM freq euc Phylogeny')\n")
    f.write('dev.off()\n')
os.system('R < ' + tmp_r + ' --no-save')

# with open('dataf', 'w') as f:
#     for elm in t:
#         f.write(str(t[elm]) + '\n')
# with open('tmpR', 'w') as f:
#     f.write("d<-read.delim('dataf',header=FALSE,sep='\\t')\n")
#     f.write("png('var.png')\n")
#     f.write('plot(density(d[,1]))\n')
#     f.write('dev.off()\n')
# os.system('R < tmpR --no-save')
