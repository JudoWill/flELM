"""Substitute residues for properties
   in flu ELMs file"""
import sys, global_settings, utils
from collections import defaultdict

with open(sys.argv[1]) as f:
    for line in f:
        (protein, st, stp, elm, seq, ELM) = line.strip().split('\t')
        if 'X' not in seq and 'B' not in seq and 'Z' not in seq and 'J' not in seq:
            print('\t'.join((protein, st, stp, 
                             elm, utils.mk_sub(seq), ELM)) 
                  + '\n')
