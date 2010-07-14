"""Given a set of uniq sequences for human & bird flu, compare
   the residue frequencies for the 2 species."""

import os
from collections import defaultdict

def get_seqs(uniq_file, use_protein):
    """Get uniq seqs on flu"""

    seqs = {}
    with open(uniq_file) as f:
        for line in f:
            protein, elmseq, elm, seq = line.strip().split('\t')
            #if protein == use_protein:
            seqs[seq] = True

    rm_seqs = {}
    for seq in seqs:
        for seq2 in seqs:
            if seq != seq2:
                if seq in seq2:
                    rm_seqs[seq2] = True
    print len(seqs), len(rm_seqs)
    #for seq in rm_seqs:
    #    del seqs[seq]
    print len(seqs)
    return seqs

def get_freqs(file):
    """Get host freqs"""

    freqs = {}
    with open(file) as f:
        for line in f:
            elmseq, freq = line.strip().split('\t')
            freqs[elmseq.split(':')[1]] = float(freq)
    return freqs

def evaluate(name, flu_seqs, host1_freqs, host2_freqs):
    count = 0
    total = 0
    for seq in flu_seqs:
        if seq in host1_freqs and seq in host2_freqs:
            if host1_freqs[seq] > host2_freqs[seq]:
                count += 1
            total += 1
    print name, count, total, float(count)/float(total)
    
protein = 'nonstructural protein 1'
dir = 'working/Jul12/'

mammal_pre = get_seqs(os.path.join(dir, 'mammal_uniq'),
                      protein)
bird_pre = get_seqs(os.path.join(dir, 'bird_uniq'),
                    protein)
mammal_control_pre = get_seqs(os.path.join(dir, 'mammal_control'),
                              protein)
bird_control_pre = get_seqs(os.path.join(dir, 'bird_control'),
                            protein)
mammal_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.init.elm_aa_freq'))
bird_host_freqs = get_freqs(os.path.join(dir, 'Gallus_gallus.init.elm_aa_freq'))

mammal = set(mammal_pre.keys()) - set(bird_pre.keys())
bird = set(bird_pre.keys()) - set(mammal_pre.keys())

mammal_control_pre2 = set(mammal_control_pre.keys()) - mammal
bird_control_pre2 = set(bird_control_pre.keys()) - bird
both_control = bird_control_pre2 & mammal_control_pre2
mammal_control = mammal_control_pre2 - both_control
bird_control = bird_control_pre2 - both_control

print 'intr', len(mammal_control & mammal)
print 'intr', len(bird_control & bird)
print 'intr', len(mammal & bird)
print 'intr', len(mammal_control & bird_control)

evaluate('MAMMAL', mammal, mammal_host_freqs, bird_host_freqs)
evaluate('MAMMAL CONTROL', mammal_control, mammal_host_freqs, bird_host_freqs)

evaluate('BIRD', bird, bird_host_freqs, mammal_host_freqs)
evaluate('BIRD CONTROL', bird_control, bird_host_freqs, mammal_host_freqs)



