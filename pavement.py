from paver.easy import *
import os, os.path

import sys
sys.path.append('.')
from local_settings import *

@task
def get_elm_patterns():
	""" Grab ELM patterns from the resource """

	sh('python get_elm_patterns.py > elm_expressions.txt')

def get_flu_seq():
	""" Grab flu protein fasta & description file from NCBI """

	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/influenza.faa.gz %s' % DATADIR)
	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/genomeset.dat.gz %s' % DATADIR)
	sh('gunzip -fq %s' % os.path.join(DATADIR, 'influenza.faa.gz'))
	sh('mv ' + os.path.join(DATADIR, 'influenza.faa') + ' '
	   + os.path.join(DATADIR, 'influenza.fa'))
	sh('gunzip -fq %s' % os.path.join(DATADIR, 'genomeset.dat.gz'))

@task
def get_host_seq():
	""" mouse, cow, dog, fish, hourse, chicken, human, rat protein seq from NCBI """
	
	genomes = ('M_musculus', 'Bos_taurus','Canis_familiaris','D_rerio',
		   'Equus_caballus', 'Gallus_gallus', 'H_sapiens', 
		   'Macaca_mulatta', 'R_norvegicus')
	bs = 'ftp.ncbi.nlm.nih.gov::genomes/'
	for genome in GENOMES:
		fname = genome+'.fa.gz'
		sh('rsync -av --size-only %(bs)s%(ome)s/protein/protein.fa.gz %(pth)s' % {'bs':bs, 
											  'ome':genome, 'pth':os.path.join(DATADIR, fname)})
		sh('gunzip -fq %s' % os.path.join(DATADIR, fname))

@task
def process_elm():
	""" Determines (and writes) the ELM dictionary"""

	FORCE_REDO = False
	
	for genome in GENOMES:
		ofile = os.path.join(RESULTSDIR, 'elmdict_'+genome+'.txt')
		ifile = os.path.join(DATADIR, genome+'.fa')
		if not os.path.exists(ofile) or FORCE_REDO:
			#only do if missing or FORCING
			sh('python makeELMdict.py -o %(out)s %(infile)s' % {'out':ofile, 
									    'infile': ifile})

@task
def elm_hist():
	""" Plot host/virus histograms of sequence counts """

	sh('python elm_hists.py '
	   + os.path.join(RESULTSDIR, 'elmdict_M_musculus.txt') + ' '
	   + 'mouse '
	   + os.path.join(RESULTSDIR, 'elmdict_M_musculus.txt') + ' '
	   + 'notmouse '
	   + '.05 '
	   + PLOTDIR)


# @task
# def get_seq():
# 	""" Grab protein fasta & description file from NCBI """

# 	# flu
# 	for afile, file_name in [['influenza.faa','influenza.fa'], 
# 				 ['genomeset.dat','genomeset.dat']]:
# 		dump_file = os.path.join(DATADIR, afile + '.gz')
# 		sh('rsync -av --size-only ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/'
# 		   + afile + '.gz '
# 		   + dump_file)
# 		sh('gunzip -dqc ' + dump_file + ' > '
# 		   + os.path.join(DATADIR,file_name))
