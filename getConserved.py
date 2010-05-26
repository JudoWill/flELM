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
import sys, utils_motif, utils_scripting

def getCounts(protein2annotation):
    virus2annotation = {}
    virus2proteinCount = {}
    for protein in protein2annotation.keys():
        protein_split = protein.split('.')
        subtype = protein_split[0]
        protein_name = protein_split[-1]
        if not virus2proteinCount.has_key(protein_name):
            virus2proteinCount[protein_name] = 0
        virus2proteinCount[protein_name] += 1
        for motif in protein2annotation[protein].keys():
            if not virus2annotation.has_key(protein_name):
                virus2annotation[protein_name] = {}
            if not virus2annotation[protein_name].has_key(motif):
                virus2annotation[protein_name][motif] = 0
            virus2annotation[protein_name][motif] += 1
    return [virus2annotation, virus2proteinCount]

def main():
    req_args = ['virus annotation file',
                'annotation tool',
                '% MSA cutoff']
    examples = ['../../Data/ProfileScan/hiv.prosite',
                'ProfileScan',
                '90']
    utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)
    
    annotation_file = sys.argv[1]
    tool = sys.argv[2]
    conserved_cutoff = float(sys.argv[3])
    
    protein2annotation = utils_motif.protein2annotation(annotation_file,
                                                        {tool:True})
    [virus2annotation, virus2proteinCount] = getCounts(protein2annotation)

    for vp in virus2annotation.keys():
        for motif in virus2annotation[vp].keys():
            percent = (float(100) * float(virus2annotation[vp][motif]) / 
                       float(virus2proteinCount[vp]))
            if percent >= conserved_cutoff:
                print vp + '\t0\t0\t' + motif + '\tseq\t' + tool
            sys.stderr.write(vp + '\t' + motif + '\t' + str(percent) + '\n')

if __name__ == '__main__': main()
