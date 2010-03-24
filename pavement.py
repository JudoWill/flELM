from paver.easy import *
import os, os.path, itertools

import sys
sys.path.append('.')
from local_settings import *
from global_settings import *
import utils

@task
def get_elm_patterns():
	""" Grab ELM patterns from the resource """

	sh('python get_elm_patterns.py > elm_expressions.txt')

@task
def get_flu_seq():
	""" Grab flu protein fasta & description file from NCBI """

	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/influenza.faa.gz %s' % DATADIR)
	sh('rsync -av ftp.ncbi.nlm.nih.gov::genomes/INFLUENZA/influenza_aa.dat.gz %s' % DATADIR)
	sh('gunzip -fq %s' % os.path.join(DATADIR, 'influenza.faa.gz'))
	sh('mv ' + os.path.join(DATADIR, 'influenza.faa') + ' '
	   + os.path.join(DATADIR, 'influenza.fa'))
	sh('gunzip -fq %s' % os.path.join(DATADIR, 'influenza_aa.dat.gz'))

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
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
			('picloud', 'c', 'Use PiCloud')])
def process_elm(options):
	"""Determines (and writes) the ELM dictionary"""
	
	c_arg = ''
	if options.process_elm.get('picloud', False): c_arg = '-c'
	
	for genome in GENOMES:
		ofile = os.path.join(RESULTSDIR, 'elmdict_'+genome+'.txt')
		ifile = os.path.join(DATADIR, genome+'.fa')
		if not os.path.exists(ofile) or options.process_elm.get('forcenew', False):
			#only do if missing or FORCING
			sh('python makeELMdict.py %(c)s -o %(out)s %(infile)s' % {'out':ofile, 
													'c':c_arg, 'infile': ifile})

@task
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
			('picloud', 'c', 'Use PiCloud')])
def process_flu(options):
	"""Determines (and writes) the ELM dictionary for inluenza"""

	c_arg = ''
	if options.process_flu.get('picloud', False): c_arg = '-c'
	

	for org in FLU_NAMES:
		fname = os.path.join(RESULTSDIR, 'flu_elmdict_'+org)
		if os.path.exists(fname) or options.process_flu.get('forcenew', False):
			continue
		#only do if missing or FORCING
		sh('python process_flu.py %(c)s %(name)s' % {'c':c_arg, 'name':org})


@task
@cmdopts([('picloud', 'c', 'Use PiCloud'),
			('forcenew', 'f', 'Force the re-creation of the result files'),])
def subsample_genomes(options):
	"""Determines (and writes) the ELM dictionary for inluenza"""

	arg = '-c' if options.subsample_genomes.get('picloud', False) else ''
	arg += ' -f' if options.subsample_genomes.get('forcenew', False) else ''


	for org in GENOMES:
		#only do if missing or FORCING
		sh('python subsample.py %(arg)s %(name)s' % {'arg':arg, 'name':org})





@task
def elm_hist():
	""" Plot host histograms of sequence frequencies of at least SEQ_FRAC_CUT """
	input_line = ''
	for genome in GENOMES:
		input_line += os.path.join(RESULTSDIR, 'elmdict_'
					   + genome + '.txt') + ' ' + genome + ' '

	sh('python elm_hists.py '
	   + input_line
	   + SEQ_FRAC_CUT + ' '
	   + PLOTDIR)

@task
def elm_hist_2():
	""" Plot host/virus histograms of sequence frequencies of at least .05 """

	input_line = ''
	for genome in GENOMES:
		input_line += os.path.join(RESULTSDIR, 'elmdict_'
					   + genome + '.txt') + ' ' + genome + ' '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_chicken') + ' chicken '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_human') + ' human '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_swine') + ' swine '
	sh('python elm_hists.py '
	   + input_line
	   + SEQ_FRAC_CUT + ' '
	   + os.path.join(PLOTDIR, 'full'))

@task
def barplot():
	""" Plot host/virus barplots for virus ELMs & sequences """

	input_line = ''
	for genome in ('H_sapiens', 'Gallus_gallus', 'Sus_scrofa'):
		input_line += os.path.join(RESULTSDIR, 'elmdict_'
					   + genome + '.txt') + ' ' + genome + ' '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_chicken') + ' chicken '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_human') + ' human '
	input_line += os.path.join(RESULTSDIR, 'flu_elmdict_swine') + ' swine '
	sh('python host_virus_barplot.py '
	   + input_line
	   + SEQ_FRAC_CUT + ' '
	   + os.path.join(PLOTDIR, 'virus_host'))

@task 
def hprdplot():
	sh('python host_virus_barplot.py '
	   + 'results/human.website.elm.elmdict '
	   + 'web '
	   + 'results/hprd_new.regex.elms.elmdict '
	   + '.01 '
	   + '/plots/hprd/')

@task
def redo_elmdict():
	""" subtract counts expected by chance """
	
	for g in GENOMES:
		sh('python get_aa_freq.py '
		   + os.path.join(DATADIR, g + '.fa ')
		   + '> ' + os.path.join(RESULTSDIR, g + '.aa_freq'))
		sh('python prob_of_seq.py '
		   + os.path.join(RESULTSDIR, g + '.aa_freq ')
		   + os.path.join(DATADIR, g + '.fa ')
		   + os.path.join(RESULTSDIR, 'elmdict_' + g + '.txt ')
		   + '> ' + os.path.join(RESULTSDIR, 'elmdict_' + g + '.redo'))

	for org in FLU_NAMES:
		sh('python get_flu_aa_freq.py '
		   + org)
		sh('python prob_of_seq.py '
		   + os.path.join(RESULTSDIR, 'flu.' + org + '.aa_freq ')
		   + org + ' '
		   + os.path.join(RESULTSDIR, 'flu_elmdict_' + org + ' ')
		   + '> ' + os.path.join(RESULTSDIR, 'flu_elmdict_' + org + '.redo'))

@task
def conserved_elms():
	for host in ['human', 'swine', 'equine', 'chicken']:
		sh('python count_flu.py '
		   + host)
		sh('python matchELMpattern.py '
		   + 'elm_expressions.txt '
		   + 'results/' + host + '.H9N2.fa '
		   + '> ' + 'results/' + host + '.H9N2.elms')
		sh('python getConserved.py '
		   + 'results/' + host + '.H9N2.elms '
		   + 'ELM '
		   + '90 '
		   + '1> results/' + host + '.H9N2.elms.90 '
		   + '2> results/' + host + '.H9N2.elms.conservation')
			
@task
def conserved_elms_2():
	for host in ['human', 'swine', 'equine', 'chicken']:
		sh('python mk_freq.py '
		   + 'results/' + host + '.H9N2.elms.90 '
		   + 'results/' + host + '.H9N2.elms '
		   + '> results/' + host + '.H9N2.elms.90.freq')
		#sh('python prob_of_seq.py '
		#   + os.path.join(RESULTSDIR, 'flu.' + host + '.aa_freq ')
		#   + host + ' '
		#   + os.path.join(RESULTSDIR, host + '.elms.90.freq ')
		#   + '> ' + os.path.join(RESULTSDIR, host + '.elms.90.freq.redo'))

@task
def serotypes():
	""" How does ELM conservation change between serotypes (ex H1N1, H2N4) """
	type2protein2gb2seq = utils.get_fluSeqs_by_serotype('swine')
	for t in type2protein2gb2seq:
		for p in type2protein2gb2seq[t]:
			seq_count = len(type2protein2gb2seq[t][p].keys())
			if seq_count > 100:
				print t + '\t' + p + '\t' + str(seq_count)

	

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
