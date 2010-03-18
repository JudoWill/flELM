""" Plot comparitive histograms
    of sequence distributions.

    Enter host ELM results,
    flu ELM results,
    and species name for flu file.

    This will output out seq count
    bar graphs for each ELM.
"""
import utils_plot, sys

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

def main(args):
    pass

if __name__ == '__main__': main(sys.argv)

