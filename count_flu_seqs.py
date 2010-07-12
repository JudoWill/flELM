"""How many flu sequences am I dealing with ? """
import utils, os
from collections import defaultdict

host = 'bird'
dir = 'working/Jul7'
years = range(2000,2011,1)

if host == 'mammal':
    hosts = ('human','swine', 'equine')
    strains = ('H1N1', 'H3N2', 'H5N1', 'H3N8')
elif host == 'bird':
    hosts = ('duck', 'chicken')
    strains = ('H9N2', 'H5N1')


for host in hosts:
    for year in years:
        for strain in strains:
            file = os.path.join(dir, 
                                '.'.join((host, 
                                          strain, str(year))) + '.elms')
            if os.path.exists(file):
                proteins = utils.get_proteins_from_elm_file(file)
                protein_counts = defaultdict(dict)
                for protein in proteins:
                    protein_class = protein.split('.')[-1]
                    protein_counts[protein_class][protein] = True
                for protein in protein_counts:
                    print('%s\t%s\t%s\t%s\t%d' %
                          (host, year, strain, 
                           protein, len(protein_counts[protein])))
