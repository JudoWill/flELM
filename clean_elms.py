"""Make new ELM pattern file that
   1) removes don't care positons from pattern ends
   2) removes patterns less than 4 residues
"""
import sys

with open(sys.argv[1]) as f:
    for line in f:
        elm, pattern = line.strip().split('\t')
        new_pattern = pattern.strip('.')
        if len(new_pattern) > 4:
            print('%s\t%s' %
                  (elm, new_pattern))
