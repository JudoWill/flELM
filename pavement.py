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

    for afile, file_name in [['influenza.faa',
                              'influenza.fa'], 
                             ['genomeset.dat',
                              'genomeset.dat']]:
        dump_file = os.path.join(DATADIR, afile + '.gz')
        sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/'
           + afile + '.gz '
           + dump_file)
        sh('gunzip -dqc ' + dump_file + ' > '
           + os.path.join(DATADIR,file_name))

@task
def get_host_seq():
    """ mouse, cow, dog, fish, hourse, chicken, human, rat, cat protein seq from NCBI """
    species = [#'M_musculus',
               'Bos_taurus',
               'Canis_familiaris',
               'D_rerio',
               'Equus_caballus',
               'Gallus_gallus',
               'H_sapiens',
               'Macaca_mulatta',
               'R_norvegicus']
    for org in species:
        dump_file = os.path.join(DATADIR, org + '.fa.gz')
        sh('wget -nc -O '
           + dump_file
           + ' ftp://ftp.ncbi.nih.gov/genomes/'
           + org + '/protein/protein.fa.gz')
        sh('gunzip -dqc ' + dump_file + ' > '
           + os.path.join(DATADIR, org + '.fa'))

    #sh('sh wgetgenome.sh')

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
