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

dir = 'working/Jul1_year'
years = range(2000,2011,1)

if host == 'mammal':
    hosts = ('human',)
    strains = ('H1N1', 'H3N2', 'H5N1', 'H3N8')
elif host == 'bird':
    hosts = ('duck', 'chicken')
    strains = ('H9N2', 'H5N1')

d = {'ELM':True}
suffix = '.elms'
utils.get_cons_elms(dir, hosts, years, strains, '70', d, out_file_elms, suffix)

d = {'ELMseq':True}
suffix = '.elmseqs'
utils.get_cons_elms(dir, hosts, years, strains, '70', d, out_file_seqs, suffix)

d = {'ELMseq':True}
suffix = '.simpleelmseqs'
utils.get_cons_elms(dir, hosts, years, strains, '70', d, out_file_simpleELMseqs, suffix)

                              
