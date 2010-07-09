#---------------------matchELMpattern.py-
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#
#----------------------------------------
""" Match regular expressions from ELM given a pattern file and fasta file.
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

def matchSeq(protein, seq, pattern2elm, pattern2regex):
    for elm_pattern in pattern2regex:
        #elm_pattern = elm2pattern[elm]#.replace('(','').replace(')','')
        #p = re.compile(elm_pattern)
        p = pattern2regex[elm_pattern]
        match = p.search(seq)
        
        # there's no need to search for more matches
        # if the match must occur at the amino end
        if 'LIG_PCNA' in pattern2elm[elm_pattern]:
            #p = re.compile("Q.[^FHWY][ILM][^P][^FHILVWYP][DHFM][FMY]..")
            p = special_p
            tempSeq = ''
            for s in seq:
                tempSeq = tempSeq + s
            offset = 0
            while match:
                printResult(protein, elm, match, tempSeq, offset)
                tempSeq = tempSeq[int(match.start())+1:]                
                offset += int( match.start() ) + 1
                match = p.search(tempSeq)
        elif elm_pattern[0] == '^' or elm_pattern[0:2]=='(^':
            if match:
                 printResult(protein, elm, match, seq, 0)
        else:
            tempSeq = ''
            for s in seq:
                tempSeq = tempSeq + s
            offset = 0
            while match:
                for elm in pattern2elm[elm_pattern]:
                    printResult(protein, elm, 
                                match, tempSeq, offset)
                tempSeq = tempSeq[int(match.start())+1:]
                offset += int( match.start() ) + 1
                match = p.search(tempSeq)

req_args = ['pattern file',
            'fasta file']
examples = ['../../Data/ELM/elm2pattern',
            '../../Data/FASTA/Human/hprd.intr.fasta']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

input_pattern_file = sys.argv[1]
fasta_file = sys.argv[2]
pattern2regex = {}
pattern2elm = defaultdict(dict)
with open(input_pattern_file) as f:
    for line in f:
        elm, pattern = line.strip().split('\t')
        pattern2elm[pattern][elm] = True
for pattern in pattern2elm:
    pattern2regex[pattern] = re.compile(pattern)

for protein, seq in utils.fasta_iter(fasta_file):
    matchSeq(protein, seq, pattern2elm, pattern2regex)


