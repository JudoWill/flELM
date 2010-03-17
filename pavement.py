from paver.easy import *

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

