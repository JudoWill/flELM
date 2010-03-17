""" Functions for making graphs
    for publication & investigation

    R must be installed
"""
import random, os
random.seed()

def mk_test_data():
    """ {} of sequence to
        count
    """

    test_data = {}
    for x in xrange(100):
        test_data[str(x)] = random.randint(0,100)
    return test_data


def create_elm_hist_tmp_file(seq2count_dict1, species_name1,
                             seq2count_dict2, species_name2):
    tmp_file = 'tmp_input' + str(random.randint(0,100))
    with open(tmp_file, 'w') as f:
        f.write('Set\tSeq\tCount\n')
        for name, a_set in [[species_name1,
                             seq2count_dict1], 
                            [species_name2,
                             seq2count_dict2]]:
            for k in a_set:
                f.write(name + '\t' 
                        + k + '\t' 
                        + str(a_set[k]) + '\n')
    return tmp_file

def elm_freq_histogram(seq2count_dict1, species_name1,
                       seq2count_dict2, species_name2,
                       out_file):
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
        f.write("ggplot(d, aes(Seq, Count, group=Set, colour=Set)) + geom_line(size=1) + opts(legend.position='none')\n")
        f.write('dev.off()\n')
    os.system('R < ' + r_file + ' --no-save')
