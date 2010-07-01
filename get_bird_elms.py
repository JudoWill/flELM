"""Find ELMs that are conserved on 
   mammal flu protiens for all
   years and strains."""
import os, global_settings, utils_motif, utils
from collections import defaultdict

dir = 'working/Jul1_year'
birds = ('duck', 'chicken')
years = range(2000,2011,1)
strains = ('H9N2', 'H5N1')

utils.get_cons_elms(dir, birds, years, strains, '70')

                              
