""" Functions for making graphs
    for publication & investigation

    R must be installed
"""
import random, os, itertools, utils
from global_settings import *
from collections import defaultdict
random.seed()

def mk_test_data():
    """ make two {}'s of seq counts """

    d1 = {'AAA':5,
          'BBB':3,
          'CCC':7,
          'DDD':1,
          'EEE':15}

    d2 = {'AAA':7,
          'BBB':10,
          'CCC':1,
          'DDD':9,
          'EEE':8}

    return [d1, d2]

def mk_random_test_data():
    """ {} of sequence to
        count
    """

    test_data = {}
    for x in xrange(100):
        test_data[str(x)] = random.randint(0,100)
    return test_data

def get_total_elm_hits(seq2count_dict):
    """ Get total number of ELM sequences
        for this ELM. 
    """

    total = 0
    for seq in seq2count_dict:
        total += seq2count_dict[seq]
    return float(total)

def create_elm_hist_tmp_file(seq2count_dict1, species_name1,
                             seq2count_dict2, species_name2):
    """ Called in elm_freq_histogram """

    tmp_file = 'tmp_input' + str(random.randint(0,100))
    #total_hits1 = get_total_elm_hits(seq2count_dict1)
    #total_hits2 = get_total_elm_hits(seq2count_dict2)
    with open(tmp_file, 'w') as f:
        f.write('Set\tSeq\tNormalizedCount\n')
        for name, a_set, in [[species_name1,
                              seq2count_dict1], 
                             [species_name2,
                              seq2count_dict2]]:
            for k in a_set:
                f.write(name + '\t' 
                        + k + '\t' 
                        + str(a_set[k]) + '\n')

                        #+ str(float(a_set[k])/total_hits) + '\n') already given %s
    return tmp_file

def elm_freq_histogram(seq2count_dict1, species_name1,
                       seq2count_dict2, species_name2,
                       out_file, elm_name):
    """ This is intended to work for 
        one ELM type to compare sequence
        distributions.

        Given a {} of raw sequence counts,
        normalize the counts and plot
        2 histograms
    """
    
    tmp_input = create_elm_hist_tmp_file(seq2count_dict1, species_name1,
                                         seq2count_dict2, species_name2)
    r_file = 'tmp' + str(random.randint(0,100))
    with open(r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("d<-read.delim('"
                + tmp_input + "', header=T, sep='\\t')\n")
        f.write("png('" + out_file + "')\n")
        #f.write("ggplot(d, aes(Seq, Count, group=Set, colour=Set)) + geom_line(size=1) + opts(legend.position='none')\n")
        #f.write("ggplot(d) + aes(x=Seq) + geom_histogram() + facet_grid(Set~.)\n")
        f.write("ggplot(d) + aes(x=Seq,y=NormalizedCount) + geom_bar(aes(fill=Set)) + facet_grid(Set~.) + opts(legend.position='none') + opts(title='" + elm_name + "')\n")
        f.write('dev.off()\n')
    os.system('R < ' + r_file + ' --no-save')
    os.system('rm ' + r_file + ' ' + tmp_input)


def create_elm_barplot_tmp_file(species2elms, elm):
    """ Called in elm_host_barplot """

    tmp_file = 'tmp_input' + str(random.randint(0,100))

    with open(tmp_file, 'w') as f:
        f.write('Set\tSeq\tNormalizedCount\n')
        for species in species2elms:
            if elm in species2elms[species]:
                for seq in species2elms[species][elm]:
                    f.write(ALIASES[species] + '\t' 
                            + seq + '\t' 
                            + str(species2elms[species][elm][seq]) + '\n')

    return tmp_file

def elm_host_barplot(species2elms, elm, out_file):
    """ This is intended to work for 
        one ELM type to compare sequence
        distributions.

        Given a {} species to {} of raw sequence counts,
        normalize the counts and plot
        one histogram per species
    """
    
    tmp_input = create_elm_barplot_tmp_file(species2elms, elm)
    r_file = 'tmp' + str(random.randint(0,100))
    with open(r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("d<-read.delim('"
                + tmp_input + "', header=T, sep='\\t')\n")
        f.write("png('" + out_file + "')\n")
        f.write("ggplot(d) + aes(x=Seq,y=NormalizedCount) + geom_bar(aes(fill=Set)) + facet_grid(Set~.) + opts(legend.position='none') + opts(title='" + elm + "')\n")
        f.write('dev.off()\n')
        f.write('q()\n')
    os.system('R < ' + r_file + ' --no-save')
    #os.system('rm ' + r_file + ' ' + tmp_input)

def distance_heatmap(x_by_y_dict, out_file):
    """ make a heatmap of x vs. y
        using the {} d[x][y] = val

        'steelblue'
    """

    tmp_input_file = 'tmp' + str(random.randint(0,100))
    tmp_r_file = 'rtmp' + str(random.randint(0,100))
    with open(tmp_input_file, 'w') as f:
        f.write('ELM\tSpeciesPair\tDistance\n')
        for x in x_by_y_dict:
            for y in x_by_y_dict[x]:
                f.write(x + '\t' + 
                        y + '\t' +
                        str(x_by_y_dict[x][y]) + '\n')
    with open(tmp_r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("d<-read.delim('"
                + tmp_input_file + "', header=T, sep='\\t')\n")
        f.write("png('" + out_file + "')\n")
        f.write("ggplot(d,aes(SpeciesPair,ELM)) + geom_tile(aes(fill=Distance),colour='white') + scale_fill_gradient(low='red',high='steelblue')\n")
        f.write('dev.off()\n')
    os.system('R < ' + tmp_r_file + ' --no-save')
    os.system('rm ' + tmp_r_file + ' ' + tmp_input_file)

def phylogeny_js(out_file, dists):
    """Make a dendrogram using Jensen-Shannon divergence."""
    
    tmp_input = 'tmp_data'
    tmp_r = 'tmp_r' + str(random.randint(0,100))
    tmp_labels = 'labels' + str(random.randint(0,100))

    js_distances = defaultdict(dict)
    for host1, host2 in itertools.combinations(dists, 2):
        js_dis = utils.jensen_shannon_dists(dists[host1],
                                            dists[host2])
        js_distances[host1][host2] = js_dis
        js_distances[host2][host1] = js_dis

    with open(tmp_input, 'w') as f:    
        for host1 in dists:
            line = ''
            for host2 in dists:
                if host1 == host2:
                    line += '0\t'
                else:
                    line += str(js_distances[host1][host2]) + '\t'
            f.write(line.strip('\t') + '\n')

    with open(tmp_labels, 'w') as f:
        f.write('\t'.join(dists.keys()) + '\n')

    with open(tmp_r, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write("library('MASS')\n")
        f.write("d<-read.delim('"
                + tmp_input
                + "',header=FALSE,sep='\\t')\n")
        f.write('dist.r<-as.dist(d)\n')
        f.write("labels.d<-read.delim('"
                + tmp_labels
                + "',header=FALSE,sep='\\t')\n")
        f.write('labels<-as.matrix(labels.d)\n')
        f.write("h<-hclust(dist.r,method='average')\n")
        f.write("png('" + out_file + "')\n")
        f.write("plot(h,hang=-1,labels=labels[1,],main='Species Dendrogram')\n")
        f.write('dev.off()\n')
    os.system('R < ' + tmp_r + ' --no-save')
    os.system('rm ' + ' '.join((tmp_r, tmp_labels, tmp_input)))

def phylogeny_euc():
    """Make a dendrogram using Euclidean distance."""

    pass
    
