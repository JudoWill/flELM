from paver.easy import *
import os, os.path

import sys
sys.path.append('.')
from local_settings import *
import utils_plot

@task
def get_elm_patterns():
	""" Grab ELM patterns from the resource """

	sh('python get_elm_patterns.py > elm_expressions.txt')

@task
def get_flu_seq():
	""" Grab flu protein fasta & description file from NCBI """

	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/influenza.faa.gz %s' % DATADIR)
	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/genomeset.dat.gz %s' % DATADIR)
	sh('gunzip -fq %s' % os.path.join(DATADIR, 'influenza.faa.gz'))
	sh('gunzip -fq %s'	% os.path.join(DATADIR, 'genomeset.dat.gz'))

@task
def get_host_seq():
	""" mouse, cow, dog, fish, hourse, chicken, human, rat protein seq from NCBI """
	genomes = ('M_musculus', 'Bos_taurus','Canis_familiaris','D_rerio',
				'Equus_caballus', 'Gallus_gallus', 'H_sapiens', 
				'Macaca_mulatta', 'R_norvegicus')
	bs = 'ftp.ncbi.nlm.nih.gov::genomes/'
	for genome in genomes:
		fname = genome+'.fa.gz'
		sh('rsync -av --size-only %(bs)s%(ome)s/protein/protein.fa.gz %(pth)s' % {'bs':bs, 
					'ome':genome, 'pth':os.path.join(DATADIR, fname)})
		sh('gunzip -fq %s' % os.path.join(DATADIR, fname))

@task
def process_elm():
	"""Determines (and writes) the ELM dictionary"""
	
	seq_files = os.listdir(DATADIR)
	
	for f in filter(lambda x: x.endswith('.fa'), seq_files):
		ofile = os.path.join(RESULTSDIR, f.split('.')[0] + '.txt')
		ifile = os.path.join(DATADIR, f)
		sh('python makeELMdict.py -o %(out)s %(infile)s' % {'out':ofile, 
																	'infile': ifile})

@task
def elm_hist():
	""" Plot host/virus histograms of sequence counts """

	sh('python elm_hists.py')
