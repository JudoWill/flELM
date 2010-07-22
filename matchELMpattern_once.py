#---------------------matchELMpattern.py-
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#
#----------------------------------------
""" Match regular expressions from ELM given a pattern file and fasta file.
Only matches a pattern in the sequence once for faster matching.
../../../school/Data/Protein_Annotations/ELM/elm2pattern
fasta file """
import string, re, sys, utils
import utils_scripting
from collections import defaultdict

special_p = re.compile("Q.[^FHWY][ILM][^P][^FHILVWYP][DHFM][FMY]..")

# add one to start b/c that part if 0 index based
# end is used for matching, and is not really the end site
def printResult(protein, elm, match, seq, offset):
    print protein + '\t' + str(offset+match.start()+1) \
          + '\t' + str(offset+match.end()) \
          + '\t' + elm + '\t' \
          + seq[int(match.start()):int(match.end())] + '\tELM'

def matchSeq(protein, seq, pattern2regex):
    for elm_pattern in pattern2regex:
        #elm_pattern = elm2pattern[elm]#.replace('(','').replace(')','')
        #p = re.compile(elm_pattern)
        p = pattern2regex[elm_pattern]
        match = p.search(seq)
        
        # there's no need to search for more matches
        # if the match must occur at the amino end
        offset = 0
        if match:
            printResult(protein, elm_pattern, 
                        match, seq, offset)

req_args = ['pattern file',
            'fasta file']
examples = ['../../Data/ELM/elm2pattern',
            '../../Data/FASTA/Human/hprd.intr.fasta']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

input_pattern_file = sys.argv[1]
fasta_file = sys.argv[2]
pattern2regex = {}
with open(input_pattern_file) as f:
    for line in f:
        elm, pattern = line.strip().split('\t')
        pattern2regex[pattern] = re.compile(pattern)

for protein, seq in utils.fasta_iter(fasta_file):
    matchSeq(protein, seq, pattern2regex)


