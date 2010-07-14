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
            seqs[seq] = True

    rm_seqs = {}
    for seq in seqs:
        for seq2 in seqs:
            if seq != seq2:
                if seq in seq2:
                    rm_seqs[seq] = True
    print len(seqs), len(rm_seqs)
    #for seq in rm_seqs:
    #    del seqs[seq]
    print len(seqs)
    return seqs

def get_freqs(file, elmdict_file):
    """Get host freqs"""

    counts = {}
    with open(elmdict_file) as f:
        for line in f:
            elmseq, seq, count, freq = line.strip().split('\t')
            counts[seq] = int(count)

    freqs = {}
    with open(file) as f:
        for line in f:
            elmseq, freq = line.strip().split('\t')
            seq = elmseq.split(':')[1]
            if counts[seq] > 500:
                freqs[seq] = float(freq)
    return freqs

def evaluate(name, flu_seqs, host1_freqs, host2_freqs):
    """Test hypothesis"""

    count = 0
    total = 0
    for seq in flu_seqs:
        f1 = float(0)
        f2 = float(0)
        if seq in host1_freqs:
            f1 = host1_freqs[seq]
        if seq in host2_freqs:
            f2 = host2_freqs[seq]
        if f1 != float(0) and f2 != float(0):
            if f1 > f2:
                count += 1
            total += 1
        # if seq in host1_freqs and seq in host2_freqs:
        #     if host1_freqs[seq] > host2_freqs[seq]:
        #         count += 1
        #     total += 1
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
mammal_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.init.elm_aa_freq'),
                              os.path.join(dir, 'elmdict_H_sapiens.init'))
bird_host_freqs = get_freqs(os.path.join(dir, 'Gallus_gallus.init.elm_aa_freq'),
                            os.path.join(dir, 'elmdict_Gallus_gallus.init'))

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



