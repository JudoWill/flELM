"""How many flu sequences am I dealing with ? 
   Make a heatplot per flu protein."""
import utils, os, random, utils, global_settings, math
from collections import defaultdict

host = 'mammal'
dir = 'working/Jul12'
years = range(2000,2011,1)

#if host == 'mammal':
hosts = ('human', 'swine', 'equine', 'chicken', 'duck')
strains = ('H1N1', 'H3N2', 'H5N1', 'H3N8', 'H9N2')
#elif host == 'bird':
#    hosts = ('duck', 'chicken')
#    strains = ('H9N2', 'H5N1')

input = 'working/input' + str(random.randint(0,100))
rfile = 'working/rfile' + str(random.randint(0,100))
outfile = 'working/Jul12/flu_seqs.png'
with open(input, 'w') as f:
    f.write('Host\tStrain\tYear\tProtein\tLogCount\n')
    for host in hosts:
        for year in years:
            for strain in strains:
                file = os.path.join(dir, 
                                    '.'.join((host, 
                                              strain, str(year))) + '.fa')
                if os.path.exists(file):
                    protein_counts = defaultdict(dict)
                    for ID, seq in utils.fasta_iter(file):
                        protein_class = ID.split('.')[-1]
                        protein_counts[protein_class][ID] = True
                    for protein in protein_counts:
                        count = len(protein_counts[protein])
                        if count > global_settings.SEQ_LIMIT:
                            val = math.log(count ,10)
                            f.write('%s\t%s\t%s\t%s\t%.10f\n' %
                                    (host, strain, str(year)[2:], 
                                     global_settings.FLU_PROTEINS[protein], 
                                     val))

with open(rfile, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("d<-read.delim('"
            + input + "', header=TRUE, sep='\\t', as.is=TRUE)\n")
    f.write("png('" + outfile + "')\n")
    f.write("ggplot(d,aes(Year,Host)) + geom_tile(aes(fill=LogCount),colour='white') + scale_fill_gradient(low='white',high='red') + facet_grid(Protein~Strain)\n")
    f.write('dev.off()\n')
os.system('R < ' + rfile + ' --no-save')
os.system('rm ' + rfile + ' ' + input)
                        
