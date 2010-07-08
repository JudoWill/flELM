""" Make a plot of host ELM sequence frequencies
    for uniq influenza ELM sequences."""
from collections import defaultdict
import utils, os, utils_graph

good_elms = utils_graph.getNodes('working/Jul7/good_phylogeny_elms')

def write_file(fname, uniq, this_host_freqs, that_host_freqs, this, that):
    used = {}
    with open(fname, 'w') as f:
        f.write('ELM\tSeq\tHost\tFreq\n')
        for elmseq in uniq:
            protein, elm, seq = elmseq.split(':')
            new_seq = utils.mk_sub(seq)
            key = elm + ':' + new_seq
            seq = new_seq
            if key != 'LIG_RGD:RGD' and key != 'LIG_AP2alpha_2:DPF' and key not in used and elm in good_elms:
                
                used[key] = True
                this_val = float(0)
                that_val = float(0)
                if key in this_host_freqs:
                    this_val = this_host_freqs[key]
                if key in that_host_freqs:
                    that_val = that_host_freqs[key]
                diffpos = max([float(0), this_val - that_val])
                diffneg = max([float(0), that_val - this_val])
                if this_val and that_val:
                    # f.write('%s\t%s\t%s\t%.10f\n' %
                    #         (elm, seq, this, this_val))
                    # f.write('%s\t%s\t%s\t%.10f\n' % 
                    #         (elm, seq, that, that_val))
                    f.write('%s\t%s\t%s\t%.10f\n' % 
                            (elm, seq, 'DiffPos', diffpos))
                    f.write('%s\t%s\t%s\t%.10f\n' % 
                            (elm, seq, 'DiffNeg', diffneg))

def get_annotations(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            protein, elm = line.strip().split('\t')
            d[protein][elm] = True
    return d

def get_important(ls1, ls2):
    """Get annotations important on flu"""

    elms = {}
    for ls in (ls1, ls2):
        for protein in ls:
            for elm in ls[protein]:
                elms[protein + ':' + elm] = True
    return elms

def get_uniq(uniq, ls1, ls2):
    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[protein + ':' + elm] = True

def get_freqs(afile, seq_percents):

    with open(afile) as f:
        for line in f:
            protein, elmseq, freq = line.strip().split('\t')
            key = ':'.join([protein, elmseq])
            seq_percents[key].append(float(freq))

def get_host_freqs(afile):
    freqs = {}
    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            freqs[elm + ':' + seq] = float(freq)
    return freqs

dir = 'working/Jul1_year'
years = range(2000,2011,1)

# bird
birds = ('chicken','duck')
strains = ('H9N2', 'H5N1')
seq_percents_bird = defaultdict(list)
total_bird = 0
for bird in birds:
    for strain in strains:
        for year in years:
            file = os.path.join(dir, '.'.join((bird, strain, str(year),
                                               'elmseqs.conservation')))
            if os.path.exists(file):
                total_bird += 1
                get_freqs(file, seq_percents_bird)

# human
hosts = ('human',)
strains = ('H5N1', 'H1N1', 'H3N2', 'H3N8')
seq_percents_human = defaultdict(list)
total_human = 0
for host in hosts:
    for strain in strains:
        for year in years:
            file = os.path.join(dir, '.'.join((host, strain, str(year),
                                               'elmseqs.conservation')))
            if os.path.exists(file):
                total_human += 1
                get_freqs(file, seq_percents_human)

human_elmseqs_file = 'working/Jul1_year/mammal_elmseqs'
bird_elmseqs_file = 'working/Jul1_year/bird_elmseqs'

human_elmseqs = get_annotations(human_elmseqs_file)
bird_elmseqs = get_annotations(bird_elmseqs_file)

uniq_human = {}
uniq_bird = {}
get_uniq(uniq_human, human_elmseqs, bird_elmseqs)
get_uniq(uniq_bird, bird_elmseqs, human_elmseqs)
impt_elms = get_important(human_elmseqs, bird_elmseqs)

# for elmseq in seq_percents_bird:
#     if elmseq in uniq:
#         print elmseq + '\t' + str(sum(seq_percents_bird[elmseq])/float(total_bird)) + '\t' + str(sum(seq_percents_human[elmseq])/float(total_human))

#human_host_file = 'working/Jun29/elmdict_H_sapiens.simple'
#chicken_host_file = 'working/Jun29/elmdict_Gallus_gallus.simple'

human_host_file = 'working/Jul7/elmdict_H_sapiens.RWinit'
chicken_host_file = 'working/Jul7/elmdict_Gallus_gallus.RWinit'

human_host_freqs = get_host_freqs(human_host_file)
chicken_host_freqs = get_host_freqs(chicken_host_file)
              
write_file('working/human_plt_host_freqs_RW', uniq_human, human_host_freqs, 
           chicken_host_freqs, 'Human', 'Chicken')
write_file('working/chicken_plt_host_freqs_RW', uniq_bird, chicken_host_freqs, 
           human_host_freqs, 'Chicken', 'Human')
