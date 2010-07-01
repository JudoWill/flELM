"""Find ELMs that are conserved on 
   mammal flu protiens for all
   years and strains."""
import os, global_settings, utils_motif, utils
from collections import defaultdict

dir = 'working/Jul1_year'
mammals = ('human', 'swine', 'horse')
years = range(2000,2011,1)
strains = ('H1N1', 'H3N2', 'H5N1', 'H3N8')
utils.get_cons_elms(dir, mammals, years, strains, '70')

                              
