""" Plot comparitive histograms
    of sequence distributions.

    Enter host ELM results,
    flu ELM results,
    and species name for flu file.

    This will output out seq count
    bar graphs for each ELM.
"""
import utils_plot, sys, os
from collections import defaultdict

def get_test_data():
    #test_data1 = utils_plot.mk_test_data()
    #test_data2 = utils_plot.mk_test_data()
    #return [test_data1, test_data2]
    return utils_plot.mk_test_data()

def test_plot_old():
    """ Makes test data """

    [test_data1, test_data2] = get_test_data()
    
    print 'Set\tSeq\tCount'
    for name, a_set in [['set1',
                         test_data1], 
                        ['set2',
                         test_data2]]:
        for k in a_set:
            print name + '\tAAAA' + k + '\t' + str(a_set[k])

def test_plot():
    [test_data1, test_data2] = get_test_data()
    utils_plot.elm_freq_histogram(test_data1, 'human',
                                  test_data2, 'virus',
                                  'test.png', 'LIG_TEST')

def get_seq2count_dict(elm_file, cutoff):
    elm2seq2count = defaultdict(dict)
    with open(elm_file) as f:
        for line in f:
            [elm, seq, count, frac_st] = line.strip().split('\t')
            frac = float(frac_st)
            if frac >= cutoff:
                elm2seq2count[elm][seq] = frac
    return elm2seq2count

def main(args):
    elm_file_1 = sys.argv[1]
    species1 = sys.argv[2]
    elm_file_2 = sys.argv[3]
    species2 = sys.argv[4]
    cutoff = float(sys.argv[5])
    plot_dir = sys.argv[6]

    elms1 = get_seq2count_dict(elm_file_1, cutoff)
    elms2 = get_seq2count_dict(elm_file_2, cutoff)

    elms = {}
    for elm_d in [elms1, elms2]:
        for elm in elm_d:
            elms[elm] = True
    for elm in elms:
        if len(elms1[elm].keys()) == 1 and len(elms2[elm].keys()) == 1:
            pass
        else:
            utils_plot.elm_freq_histogram(elms1[elm], species1,
                                          elms2[elm], species2,
                                          os.path.join(plot_dir,
                                                       elm + '.png'),
                                          elm)

if __name__ == '__main__': main(sys.argv)

