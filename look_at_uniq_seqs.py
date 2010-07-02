# it is miss named now
from collections import defaultdict

def get_annotations(afile):
    d = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            protein, elm = line.strip().split('\t')
            d[protein][elm] = True
    return d

def get_uniq(uniq, ls1, ls2):
    for protein in ls1:
        for elm in ls1[protein]:
            if elm not in ls2[protein]:
                uniq[elm] = True

human_elmseqs_file = 'working/Jul1_year/mammal_elmseqs'
bird_elmseqs_file = 'working/Jul1_year/bird_elmseqs'

human_elmseqs = get_annotations(human_elmseqs_file)
bird_elmseqs = get_annotations(bird_elmseqs_file)

uniq = {}
get_uniq(uniq, human_elmseqs, bird_elmseqs)
get_uniq(uniq, bird_elmseqs, human_elmseqs)
print len(uniq)


