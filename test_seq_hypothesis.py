""" Hypothesis: the virus sequences 
    that are most utilized will
    be the most utilized sequences
    for the correct host.
"""
import sys, utils_stats, utils_scripting, utils_motif, utils
from collections import defaultdict

req_args = ['elmdict',
            'virus ELMs',
            'conserved virus ELMs',
            'out file']
examples = ['results/elmdict_H_sapiens.txt',
            'human.H1N1.elms',
            'human.H1N1.elms.90',
            'hypoth_human']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

def get_freq(afile):
    """get host frequency"""
    freq = {}
    with open(afile) as f:
        for line in f:
            elm, seq, num, freq_st = line.strip().split('\t')
            if not elm in freq:
                freq[elm] = {}
            freq[elm][seq] = float(freq_st)
    return freq

def get_best_seq(seqs):
    ls = []
    for seq in seqs:
        ls.append([seqs[seq],seq])
    ls.sort()
    if len(ls) > 1:
        if ls[-1][0] == ls[-2][0]:
            print 'TIE'
    print ls
    return ls[-1]

def getProtein2elm2seq(elms, conserved):
    """get virus ELM seqs for conserved ELMs"""
    protein2elm2seq = {}
    proteinCounts = defaultdict(utils.init_zero)
    for protein_id in elms:
        protein = protein_id.split('.')[1]
        proteinCounts[protein] += 1
        if protein in conserved:
            for elm in elms[protein_id]:
                if elm in conserved[protein]:
                    if not protein in protein2elm2seq:
                        protein2elm2seq[protein] = {}
                    if not elm in protein2elm2seq[protein]:
                        protein2elm2seq[protein][elm] = {}
                    for [st, stp, seq] in elms[protein_id][elm]:
                        if not seq in protein2elm2seq[protein][elm]:
                            protein2elm2seq[protein][elm][seq] = 0
                        protein2elm2seq[protein][elm][seq] += 1
    return (protein2elm2seq, proteinCounts)

host_file = sys.argv[1]
elm_file = sys.argv[2]
cons_file = sys.argv[3]
ofile = sys.argv[4]

elms = utils_motif.protein2annotation(elm_file,
                                      {'ELM':True})
cons_elms = utils_motif.protein2annotation(cons_file,
                                           {'ELM':True})
host_freqs = get_freq(host_file)
(elm_counts, protein_counts) = getProtein2elm2seq(elms, cons_elms)

lines = ''
for protein in elm_counts:
    protein_count = float(protein_counts[protein])
    for elm in elm_counts[protein]:
        if elm in host_freqs:
            virus_freqs = []
            non_virus_freqs = []
            found_seqs = {}
            for seq in elm_counts[protein][elm]:
                if float(elm_counts[protein][elm][seq])/protein_count > float(.9):
                    if seq in host_freqs[elm]:
                        virus_freqs.append(host_freqs[elm][seq])
                        found_seqs[seq] = True
                    else:
                        virus_freqs.append(float(0))
                else:
                     if seq in host_freqs[elm]:
                         non_virus_freqs.append(host_freqs[elm][seq])
                         found_seqs[seq] = True
                     else:
                         non_virus_freqs.append(float(0))
            for seq in host_freqs[elm]:
                if not seq in found_seqs:
                    non_virus_freqs.append(host_freqs[elm][seq])
            #line = ''
            if len(virus_freqs) > 2 and len(non_virus_freqs) > 2:
                lines += protein + '\t'+ elm + '\t'+ str(utils_stats.wilcox_gtr(virus_freqs, non_virus_freqs)) + '\t' + str(utils_stats.wilcox_less(virus_freqs, non_virus_freqs)) + '\n'
            #else:
            #    line = protein + '\t'+ elm + '\t'+ 'NO_DATA(' + str(len(virus_freqs)) + ',' + str(len(non_virus_freqs)) + ')'
            #lines += line + '\n'
with open(ofile,'w') as f:
    f.write(lines)
        




