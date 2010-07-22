"""Gather ELM sequence hits & construct new
   simple patterns for rescanning"""
import sys, global_settings, os, utils

def acc_hits(hits, afile):
    """Grab all ELM,seq pairs.
       Convert seqs to simple seqs."""

    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            new_seq = utils.mk_sub(seq)
            if new_seq != 'NA':
                hits[elm + ':' + new_seq] = True

in_dir = sys.argv[1]
elms_file = sys.argv[2]

hits = {}
for g in ('Gallus_gallus', 'H_sapiens'):
    afile = os.path.join(in_dir, 'elmdict_'
                         + g + '.init')
    acc_hits(hits, afile)

elms = {}
seqs = {}
for elmseq in hits:
    elm, seq = elmseq.split(':')
    elms[elm] = True
    seqs[seq] = True

with open(elms_file, 'w') as f:
    for elm in elms:
        f.write(elm + '\tpattern\n')

for seq in seqs:
    print seq + '\t' + seq
