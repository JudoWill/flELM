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
    return seqs

def get_freqs(flu_seqs, file):
    """Get host freqs for uniq flu seqs"""

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
            if host1_freqs[seq] > float(1.2)*host2_freqs[seq]:
                count += 1
            total += 1
    print name, count, total, float(count)/float(total)
    
protein = 'nonstructural protein 1'
dir = 'working/Jul12/'

mammal = get_seqs(os.path.join(dir, 'mammal_uniq'),
                  protein)
bird = get_seqs(os.path.join(dir, 'bird_uniq'),
                protein)
mammal_control = get_seqs(os.path.join(dir, 'mammal_control'),
                          protein)
bird_control = get_seqs(os.path.join(dir, 'bird_control'),
                        protein)
mammal_host_freqs = get_freqs(mammal, 
                              os.path.join(dir, 'H_sapiens.init.elm_aa_freq'))
bird_host_freqs = get_freqs(bird, 
                            os.path.join(dir, 'Gallus_gallus.init.elm_aa_freq'))

evaluate('MAMMAL', mammal, mammal_host_freqs, bird_host_freqs)
evaluate('MAMMAL CONTROL', mammal_control, mammal_host_freqs, bird_host_freqs)

evaluate('BIRD', bird, bird_host_freqs, mammal_host_freqs)
evaluate('BIRD CONTROL', bird_control, bird_host_freqs, mammal_host_freqs)



