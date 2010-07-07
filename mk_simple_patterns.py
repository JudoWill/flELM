"""Gather ELM sequence hits & construct new
   simple patterns for rescanning"""
import sys, global_settings, os, utils

def acc_hits(hits, afile):
    """Grab all ELM,seq pairs.
       Convert seqs to simple seqs."""

    with open(afile) as f:
        for line in f:
            elm, seq, count, freq = line.strip().split('\t')
            if 'X' not in seq:
                hits[elm + ':' + utils.mk_sub(seq)] = True

in_dir = sys.argv[1]
elms_file = sys.argv[2]

hits = {}
for g in global_settings.TEST_GENOMES:
    afile = os.path.join(in_dir, 'elmdict_'
                         + g + '.init')
    acc_hits(hits, afile)

elms = {}
for elmseq in hits:
    print elmseq + '\t' + elmseq.split(':')[1]
    elms[elmseq.split(':')[0]] = True

with open(elms_file, 'w') as f:
    for elm in elms:
        f.write(elm + '\tpattern\n')

