#------------------getConserved.py-
#
#  Perry Evans
#
#  2009
#
#----------------------------------
""" Print ELM conservation. 
    normal out gets conserved ELMs
    error out gets all ELMs percentages 
"""
import warnings
warnings.simplefilter('ignore')
import sys, utils_scripting
from collections import defaultdict

def getCounts(annotation_file):
    """Accumulate total proteins and sequences with hits"""

    virus2annotation = {}
    virus2proteinCount = defaultdict(dict)
    pattern2elm = defaultdict(dict)
    with open(annotation_file) as f:
        for line in f:
            protein, st, stp, motif, seq, junk = line.strip().split('\t')
            protein_split = protein.split('.')
            subtype = protein_split[0]
            protein_name = protein_split[-1]
            elm, pattern = motif.split(':')
            virus2proteinCount[protein_name][subtype] = True
            if protein_name not in virus2annotation:
                virus2annotation[protein_name] = {}
            if pattern not in virus2annotation[protein_name]:
                virus2annotation[protein_name][pattern] = {}
            virus2annotation[protein_name][pattern][subtype] = True
            pattern2elm[pattern][elm] = True
            
    return [virus2annotation, virus2proteinCount, pattern2elm]

def main():
    req_args = ['virus annotation file',
                '% MSA cutoff']
    examples = ['../../Data/ProfileScan/hiv.prosite',
                '90']
    utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)
    
    annotation_file = sys.argv[1]
    conserved_cutoff = float(sys.argv[2])
    
    [virus2annotation, virus2proteinCount,
     pattern2elm] = getCounts(annotation_file)

    for vp in virus2annotation.keys():
        for pattern in virus2annotation[vp]:
            percent = (float(100) * float(len(virus2annotation[vp][motif])) / 
                       float(len(virus2proteinCount[vp])))
            for elm in pattern2elm[pattern]:
                motif = elm + ':' + pattern
                if percent >= conserved_cutoff:
                    print vp + '\t0\t0\t' + motif + '\tseq\t' + tool
                sys.stderr.write(vp + '\t' + motif + '\t' + str(percent) + '\n')

if __name__ == '__main__': main()
