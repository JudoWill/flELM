"""Find ELMs that are conserved on 
   mammal flu protiens for all
   years and strains. Also do this
   for ELM sequences."""
import os, global_settings, utils_motif, utils, sys
from collections import defaultdict

out_file_elms = sys.argv[1]
out_file_seqs = sys.argv[2]
out_file_simpleELMseqs = sys.argv[3]
host = sys.argv[4] # bird|mammal

dir = 'working/Jul7'
years = range(2000,2011,1)

if host == 'mammal':
    hosts = ('human',)
    strains = ('H5N1',)
elif host == 'bird':
    hosts = ('chicken',)
    strains = ('H5N1',)

d = {'ELM':True}
suffix = '.elms'
utils.get_cons_elms(dir, hosts, years, strains, '90', d, out_file_elms, suffix)

# d = {'ELMseq':True}
# suffix = '.elmseqs'
# utils.get_cons_elms(dir, hosts, years, strains, '70', d, out_file_seqs, suffix)

# d = {'ELMseq':True}
# suffix = '.simpleelmseqs'
# utils.get_cons_elms(dir, hosts, years, strains, '70', d, out_file_simpleELMseqs, suffix)

                              
