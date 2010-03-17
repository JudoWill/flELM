from paver.easy import *
import os, os.path

import sys
sys.path.append('.')
from local_settings import *

@task
def get_elm_patterns():
    """ Grab ELM patterns from the resource """

    sh('python get_elm_patterns.py > elm_expressions.txt')

@task
def get_flu_seq():
    """ Grab flu protein fasta & description file from NCBI """

    sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/influenza.faa.gz ./')
    sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/genomeset.dat.gz ./')
    sh('gunzip influenza.faa.gz')
    sh('gunzip genomeset.dat.gz')

@task
def get_host_seq():
    """ mouse, cow, dog, fish, hourse, chicken, human, rat protein seq from NCBI """
    sh('sh wgetgenome.sh')


@task
def process_elm():
	"""Determines (and writes) the ELM dictionary"""
	
	seq_files = os.listdir(DATADIR)
	
	for f in filter(lambda x: x.endswith('.fa'), seq_files):
		ofile = os.path.join(RESULTSDIR, f.split('.')[0] + '.txt')
		ifile = os.path.join(DATADIR, f)
		sh('python makeELMdict.py -o %(out)s %(infile)s' % {'out':ofile, 
														'infile': ifile})