"""Convert raw sequences to simplified versions
   for use in scanning with simple patterns
   (those made with utils.mk_sub)"""
import sys, global_settings, os, utils

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(out_file, 'w') as f:
    for ID, seq in utils.fasta_iter(in_file):
        f.write('>' + ID + '\n')
        f.write(utils.mk_sub(seq) + '\n')
