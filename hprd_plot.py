import utils, utils_plot, utils_distance

d = {'web':utils.get_seq2count_dict('results/human.website.elm.elmdict',
                                    float(.01)),
     'regex':utils.get_seq2count_dict('results/hprd_new.regex.elms.elmdict',
                                      float(.01))}
elms = utils_distance.get_elements(d['web'], d['regex'])

for elm in elms:
    utils_plot.elm_host_barplot(d, elm, 'plots/hprd/' + elm + '.png')
