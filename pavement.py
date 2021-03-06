from paver.easy import *
import os, os.path, itertools

import sys
sys.path.append('.')
from local_settings import *
from global_settings import *
import utils

@task
def roundup():
    """Parse human/chicken roundup & get FASTA"""

    sh('python parse_roundup.py '
       + 'working/Jul22/human_chicken.roundup '
       + 'working/Jul22/human_chicken.roundup.parsed')
    for species, name in (("'Homo sapiens'", 'H_sapiens'),
                          ("'Gallus gallus'", 'Gallus_gallus')):
        sh('python get_protein_seq_for_gi.py '
           + 'working/Jul22/human_chicken.roundup.parsed '
           + species + ' '
           + 'working/Jul22/fasta_all/' + name + '.fa')
    sh('python mk_single_member_clusters.py '
       + 'working/Jul22/human_chicken.roundup.parsed '
       + 'working/Jul22/fasta_all/ '
       + 'working/Jul22/fasta_use/')

@task
def check_hypoth():
    """Do individual flu proteins match host/flu hypoth?"""
    for protein in FLU_PROTEINS_LTD:
        sh('python look_at_polPA_all.py '
           + "'" + protein + "'")

@task
def rewrite_simple_elmdict_len():
    """After reducing the aminoacid space, and scanning reduced fasta, I need to redo the elmdicts produced. Split up ELM hits by length"""
    
    for g in TEST_GENOMES:
        sh('python rewrite_simple_elmdict_len.py '
           + os.path.join('working', 'Jul12', 'elmdict_' + g + '.init ')
           + '> ' + os.path.join('working', 'Jul12', 'elmdict_' + g + '.RWlenInit'))

@task
def rewrite_simple_elmdict():
    """After reducing the aminoacid space, and scanning reduced fasta, I need to redo the elmdicts produced"""
    
    for g in TEST_GENOMES:
        sh('python rewrite_simple_elmdict.py '
           + os.path.join('working', 'Jul12', 'elmdict_' + g + '.init ')
           + '> ' + os.path.join('working', 'Jul12', 'elmdict_' + g + '.RWinit'))
        
@task
def mk_simple_elm_patterns():
    """Look at elmdict hits, and make new simple patterns from them, and make new fasta"""

    sh('python mk_simple_patterns.py '
       + 'working/Jul22/ '
       + 'working/Jul22/use_elms '
       + '> working/Jul22/simple_patterns')
    for g in ('H_sapiens', 'Gallus_gallus'):
        sh('python mk_simple_fasta.py '
           + 'working/Jul22/fasta_use/' + g + '.fa '
           + 'working/Jul22/' + g + '.fa')
    # copy from Jul19
    # host_strains = [['human','H1N1'],
    #                 ['human','H3N2'],
    #                 ['human','H5N1'],
                    
    #                 ['swine','H3N2'],
    #                 ['swine','H1N1'],
                    
    #                 ['equine','H3N8'],
                    
    #                 ['chicken','H9N2'],
    #                 ['chicken','H5N1'],
                    
    #                 ['duck','H9N2'],
    #                 ['duck','H5N1']]
    # for host, strain in host_strains:
    #     for year in xrange(2000, 2011):
    #         sh('python mk_simple_fasta.py '
    #            + 'working/Jul22/fasta_use' + host + '.' + strain + '.' + str(year) + '.fa '
    #            + 'working/Jul22/' + host + '.' + strain + '.' + str(year) + '.fa')

@task
def get_mammal_bird_elms():
    """Find ELMs conserved on all mammal
       or bird flu proteins"""

    sh('python get_host_elms.py '
       + 'working/Jul7/mammal_elms '
       + "working/Jul7/mammal_elmseqs "
       + 'working/Jul7/mammal_simpleelmseqs mammal')
    sh('python get_host_elms.py '
       + 'working/Jul7/bird_elms '
       + "working/Jul7/bird_elmseqs "
       + 'working/Jul7/bird_simpleelmseqs bird')

@task 
def simplify_elmdicts():
    """Substitute residues for properties"""

    for host in TEST_GENOMES:
        sh('python simplify_elmdict.py '
           + 'working/Jun29/elmdict_' + host + '.init '
           + '> working/Jun29/elmdict_' + host + '.simple')

@task
def individual_elms():
    """Make host phylogeny for each ELM"""

    elms = {}
    with open('elm_expressions.txt') as f:
        for line in f:
            elm, exp = line.strip().split('\t')
            elms[elm] = True
            
    for elm in elms:
        with open('tmpELM', 'w') as f:
            f.write(elm + '\tstuff\n')
        try:
            sh('python js_elmSeqDist_hosts.py '
               + 'NA '
               + 'working/Jun30/ '
               + elm + '.simpleELMs.png '
               + "0 0 tmpELM '.simple'")
        except: pass
    sh('rm tmpELM')

@task
def find_best_elms():
    """Look for highest ELM frequency that recovers host phylogeny"""

    results_dir = 'working/Jun28/'
    cut = '.00001'
    use_elms_file = results_dir + 'use_elms' + cut

    sh('python threshold_elms.py '
       + cut + ' '
       + results_dir + ' '
       + '> ' + use_elms_file)
    sh('python js_elmDist_host.py '
       + results_dir + ' '
       + results_dir + 'js_host_elmDist_phylogeny' + cut + '.png '
       + 'F '
       + use_elms_file)
    sh('python js_elmSeqDist_hosts.py '
       + 'NA '
       + results_dir + ' '
       + 'js_host_elmSeqDist_phylogeny' + cut + '.png '
       + '2 3 '
       + use_elms_file)

@task
def get_all_roundup_seqs_ncbi():
    """use NCBI eutils to grab protein seqs for roundup orthologs"""

    for species, name in (("'Homo sapiens'", 'H_sapiens'),
                          ("'Mus musculus'", 'M_musculus'),
                          ("'Rattus norvegicus'", 'R_norvegicus'),
                          ("'Pan troglodytes'", 'Pan_troglodytes'),
                          ("'Bos taurus'", 'Bos_taurus'),
                          ("'Taeniopygia guttata'", 'Taeniopygia_guttata'),
                          ("'Gallus gallus'", 'Gallus_gallus'),
                          ("'Canis familiaris'", 'Canis_familiaris')):
        sh('python get_protein_seq_for_gi.py '
           + 'results/Homo_Mus_Pan_Rat_Bos_Can_Gal_Tae_Dan_Mac.roundup.parsed '
           + species + ' '
           + 'data/roundup_all/' + name + '.fa')

@task
def get_paper_roundup_seqs_ncbi():
    """use NCBI eutils to grab protein seqs for roundup orthologs"""

    for species, name in (("'Homo sapiens'", 'H_sapiens'),
                          ("'Mus musculus'", 'M_musculus'),
                          ("'Pan troglodytes'", 'Pan_troglodytes'),
                          ("'Sus scrofa'", 'Sus_scrofa'),
                          ("'Taeniopygia guttata'", 'Taeniopygia_guttata'),
                          ("'Gallus gallus'", 'Gallus_gallus'),
                          ("'Equus caballus'", 'Equus_caballus')):
        sh('python get_protein_seq_for_gi.py '
           + 'results/roundup_paper/roundup.parsed '
           + species + ' '
           + 'data/roundup_paper/' + name + '.fa')

@task
def get_mammal_roundup_seqs_ncbi():
    """use NCBI eutils to grab protein seqs for roundup orthologs"""

    for species, name in (("'Homo sapiens'", 'H_sapiens'),
                          ("'Mus musculus'", 'M_musculus'),
                          ("'Pan troglodytes'", 'Pan_troglodytes'),
                          ("'Sus scrofa'", 'Sus_scrofa'),
                          ("'Equus caballus'", 'Equus_caballus')):
        sh('python get_protein_seq_for_gi.py '
           + 'working/Jun28_mammals/mammal_roundup_clusters '
           + species + ' '
           + 'working/Jun28_mammals/' + name + '.fa')

@task
def other_virus_elms():
	"""format HIV/HCV ELMs"""

	sh('python convert_elm_hits.py '
	   'results/HIV1.clean.elms '
	   + '> results/flu_elmdict_HIV')
	sh('python convert_elm_hits.py '
	   '../../HCVhhp/data/HCV.elms '
	   + '> results/flu_elmdict_HCV')
	sh('python get_cons_elms.py '
	   + 'results/HIV1.clean.elms '
	   + '70 '
	   + '> results/HIV.all.elms.70.controled')
	sh('python get_cons_elms.py '
	   + '../../HCVhhp/data/HCV.elms '
	   + '70 '
	   + '> results/HCV.all.elms.70.controled')

@task
def get_fail_elms():
	host_strains = [['human','H1N1'],
			['human','H3N2'],
			['human','H5N1'],
			
			['swine','H3N2'],
			['swine','H1N1'],

			['equine','H3N8'],
			
			['chicken','H9N2'],
			['chicken','H5N1'],

			['duck','H9N2'],
			['duck','H5N1']]
	for host,strain in host_strains:
		for cint in (70,80,90):
			c = str(cint)
			sh('python findGoodProteins.py '
			   + host + '.' + strain + '.elms.' + c + ' '
			   + '1 '
			   + '> results/' + host + '.' + strain + '.elms.'
			   + c + '.controled')

@task
def host_elms():
    """Find host ELM freqs and redo freqs"""

    for genome in TEST_GENOMES:
        sh('python get_aa_freq.py '
           + 'data/roundup_paper/' + genome + '.fa '
           + '> working/Jun28/' + genome + '.aa_freq')
        sh('python mk_aa_freq.py '
           + 'data/roundup_paper/' + genome + '.fa '
           + 'working/Jun28/elmdict_' + genome + '.init '
           + 'working/Jun28/' + genome + '.init.elm_aa_freq')
        sh('python prob_of_seq.py '
           + os.path.join('working', 'Jun28', genome + '.aa_freq ')
           + os.path.join('data', 'roundup_paper', genome + '.fa ')
           + os.path.join('working', 'Jun28', 'elmdict_' + genome + '.init ')
           + '> ' + os.path.join('working', 'Jun28', 'elmdict_' + genome + '.redo'))
        sh('python mk_aa_freq.py '
           + 'data/roundup_paper/' + genome + '.fa '
           + 'working/Jun28/elmdict_' + genome + '.redo '
           + 'working/Jun28/' + genome + '.redo.elm_aa_freq')
        
@task
def elm_aa_freqs_roundup():
	for genome in TEST_GENOMES:
            sh('python mk_aa_freq.py '
               + 'data/roundup_paper/' + genome + '.fa '
               + 'results/roundup_all/elmdict_' + genome + '.simple '
               + 'results/roundup_all/' + genome + '.init.elm_aa_freq')

@task
def elm_aa_freqs_Jul20():
    for genome in ('Gallus_gallus', 'H_sapiens'):
        sh('python mk_aa_freq.py '
           + 'working/Jul22/' + genome + '.fa '
           + 'working/Jul22/elmdict_' + genome + '.simple '
           + 'working/Jul22/' + genome + '.simple.elm_aa_freq')

#conserved_elms -c 90
@task 
@cmdopts([('cutoff=', 'c', '% cutoff'),])
def ratios():
	"""test hypothesis"""

	cut = options.ratios.get('cutoff')
	for genome in ('H_sapiens', 'Gallus_gallus', 'Sus_scrofa', 'Taeniopygia_guttata'):
		sh('python get_ratios.py '
		   + 'results/elmdict_' + genome + '.redo '
		   + '> ' + genome + '.redo.ratio')
	sh('python triplet.py > new_data.tab')
	sh('python get_mammal_bird_diffs.py '
	   + 'new_data.tab '
	   + cut)

@task
def recount_elm_seqs():
	"""count ELM seqs based on fixed positions"""

	for genome in ('H_sapiens', 'Gallus_gallus', 'Sus_scrofa'):
		sh('python reCountELMseqs.py '
		   + 'results/elmdict_' + genome + '.redo > '
		   + 'results/elmdict_' + genome + '.redo.remakeELMs')

@task
def get_elm_patterns():
	""" Grab ELM patterns from the resource """

        out_file = 'data/elm_expressions.txt'
	sh('python get_elm_patterns.py > '
           + out_file)

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
	
	#genomes = ('Drosophila_melanogaster')
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
	
	for genome in ('H_sapiens', 'Gallus_gallus'):
		ofile = os.path.join('working', 'Jul22', 
                                     'elmdict_'+genome+'.init')
		ifile = os.path.join('data', genome+'.fa')
                st_elm_file = os.path.join('data', 'elm_expressions.txt')
                elms = {}
                with open(st_elm_file) as f:
                    for line in f:
                        elm, pattern = line.strip().split('\t')
                        elms[elm] = pattern
                elm_files = []
                size = 50
                if len(elms) > size:
                    counter = 0
                    for chunk in utils.chunks(elms.keys(), size):
                        new_elm_file = 'working/elm_tmp_file' + str(counter)
                        elm_files.append(new_elm_file)
                        with open(new_elm_file, 'w') as f:
                            for elm in chunk:
                                f.write(elm + '\t' + elms[elm] + '\n')
                        counter += 1
                else:
                    elm_files.append(st_elm_file)
                
                if not os.path.exists(ofile) or options.process_elm.get('forcenew', False):
                    counter = 0
                    for elmfile in elm_files:
                        #only do if missing or FORCING
                        sh('python makeELMdict.py %(c)s -o %(out)s %(infile)s %(elm)s' % {'out':ofile, 
                                                                                          'c':c_arg, 'infile': ifile,
                                                                                          'elm': elmfile})
                        sh('mv ' + ofile + ' ' + ofile + str(counter))
                        counter += 1

@task
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
	  ('picloud', 'c', 'Use PiCloud')])
def process_elm_simple(options):
	"""Determines (and writes) the ELM dictionary"""
	
	c_arg = ''
	if options.process_elm_simple.get('picloud', False): c_arg = '-c'
	
	for genome in ('H_sapiens', 'Gallus_gallus'):
		ofile = os.path.join('working', 'Jul22', 
                                     'elmdict_'+genome+'.simple')
		ifile = os.path.join('working', 'Jul22', genome+'.fa')
                st_elm_file = os.path.join('working', 'Jul22', 
                                           'simple_patterns')
                elms = {}
                with open(st_elm_file) as f:
                    for line in f:
                        elm, pattern = line.strip().split('\t')
                        elms[elm] = pattern
                elm_files = []
                size = 1000
                if len(elms) > size:
                    counter = 0
                    for chunk in utils.chunks(elms.keys(), size):
                        new_elm_file = 'working/elm_tmp_file' + str(counter)
                        elm_files.append(new_elm_file)
                        with open(new_elm_file, 'w') as f:
                            for elm in chunk:
                                f.write(elm + '\t' + elms[elm] + '\n')
                        counter += 1
                else:
                    elm_files.append(st_elm_file)
                
                if not os.path.exists(ofile) or options.process_elm_simple.get('forcenew', False):
                    counter = 0
                    for elmfile in elm_files:
                        #only do if missing or FORCING
                        sh('python makeELMdict.py %(c)s -o %(out)s %(infile)s %(elm)s' % {'out':ofile, 
                                                                                          'c':c_arg, 'infile': ifile,
                                                                                          'elm': elmfile})
                        sh('mv ' + ofile + ' ' + ofile + str(counter))
                        counter += 1

@task
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
	  ('picloud', 'c', 'Use PiCloud')])
def process_elm_roundup(options):
	"""Determines (and writes) the ELM dictionary"""
	
	c_arg = ''
	if options.process_elm_roundup.get('picloud', False): c_arg = '-c'
	
	for genome in MAMMALS2:
		ofile = os.path.join('working', 'Jun28_mammals', 'elmdict_'+genome+'.init')
		ifile = os.path.join('working', 'Jun28_mammals', genome+'.fa')
		if not os.path.exists(ofile) or options.process_elm_roundup.get('forcenew', False):
			# only do if missing or FORCING
			sh('python makeELMdict.py %(c)s -o %(out)s %(infile)s' % {'out':ofile, 
										  'c':c_arg, 'infile': ifile})

@task
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
	  ('picloud', 'c', 'Use PiCloud')])
def process_elm_roundup_sampled(options):
	"""Determines (and writes) the ELM dictionary for a run of sampled sequences from roundup orthologs"""
	
        results_dir = 'working/runs/Jun25_2/'
	c_arg = ''
	if options.process_elm_roundup_sampled.get('picloud', False): c_arg = '-c'
	
        sh('python sample_from_clusters.py '
           + 'results/Homo_Mus_Pan_Rat_Bos_Can_Gal_Tae_Dan_Mac.roundup.parsed '
           + 'data/roundup_all/ '
           + results_dir)
        
	for genome in ('H_sapiens', 'M_musculus', 'Pan_troglodytes', 
                       'R_norvegicus', 'Gallus_gallus', 'Taeniopygia_guttata',
                       'Canis_familiaris', 'Bos_taurus'):
		ofile = os.path.join(results_dir, 'elmdict_'+genome+'.init')
		ifile = os.path.join(results_dir, genome+'.fa')
		if not os.path.exists(ofile) or options.process_elm_roundup_sampled.get('forcenew', False):
			# only do if missing or FORCING
			sh('python makeELMdict.py %(c)s -o %(out)s %(infile)s' % {'out':ofile, 
										  'c':c_arg, 'infile': ifile})

@task
@cmdopts([('forcenew', 'f', 'Force the re-creation of the result files'),
	  ('picloud', 'c', 'Use PiCloud')])
def process_elm_roundup_single(options):
	"""Determines (and writes) the ELM dictionary for single roundup ortholog clusters"""
	
        results_dir = 'working/runs/Jun25/'
	c_arg = ''
	if options.process_elm_roundup_single.get('picloud', False): c_arg = '-c'
	
        sh('python get_single_clusters.py '
           + 'results/Homo_Mus_Pan_Rat_Bos_Can_Gal_Tae_Dan_Mac.roundup.parsed '
           + 'data/roundup_all/ '
           + results_dir)
        
	for genome in ('H_sapiens', 'M_musculus', 'Pan_troglodytes', 
                       'R_norvegicus', 'Gallus_gallus', 'Taeniopygia_guttata',
                       'Canis_familiaris', 'Bos_taurus'):
		ofile = os.path.join(results_dir, 'elmdict_'+genome+'.init')
		ifile = os.path.join(results_dir, genome+'.fa')
		if not os.path.exists(ofile) or options.process_elm_roundup_single.get('forcenew', False):
			# only do if missing or FORCING
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
def redo_elmdict_2():
	""" subtract counts expected by chance, allow negatives """
	
	for g in GENOMES:
		sh('python prob_of_seq_wNeg.py '
		   + os.path.join(RESULTSDIR, g + '.aa_freq ')
		   + os.path.join(DATADIR, g + '.fa ')
		   + os.path.join(RESULTSDIR, 'elmdict_' + g + '.txt ')
		   + '> ' + os.path.join(RESULTSDIR, 'elmdict_' + g + '.redoWNeg'))

@task
def freqs():
     for g in GENOMES:
          ofile = os.path.join(RESULTSDIR, g + '.aa_freq')
	  if not os.path.exists(ofile):
	       sh('python get_aa_freq.py '
		  + os.path.join(DATADIR, g + '.fa ')
		  + '> ' + ofile)
	  ofile =  os.path.join(RESULTSDIR, g + '.diAA_freq')
	  if not os.path.exists(ofile):
               sh('python get_diAA_freq.py '
		  + os.path.join(DATADIR, g + '.fa ')
		  + '> ' + ofile)

@task
def redo_elmdict():
	""" subtract counts expected by chance, no negatives """
	
	for g in GENOMES:
		sh('python get_aa_freq.py '
		   + os.path.join(DATADIR, g + '.fa ')
		   + '> ' + os.path.join('working', 'Jun28_startOver', g + '.init.aa_freq'))
		sh('python prob_of_seq.py '
		   + os.path.join('working', 'Jun28_startOver', g + '.aa_freq ')
		   + os.path.join(DATADIR, g + '.fa ')
		   + os.path.join('working', 'Jun28_startOver', 'elmdict_' + g + '.init ')
		   + '> ' + os.path.join(RESULTSDIR, 'elmdict_' + g + '.redo'))
                sh('python get_aa_freq.py '
		   + os.path.join(DATADIR, g + '.fa ')
		   + '> ' + os.path.join('working', 'Jun28_startOver', g + '.redo.aa_freq'))

	# for org in FLU_NAMES:
	# 	sh('python get_flu_aa_freq.py '
	# 	   + org)
	# 	sh('python prob_of_seq.py '
	# 	   + os.path.join(RESULTSDIR, 'flu.' + org + '.aa_freq ')
	# 	   + org + ' '
	# 	   + os.path.join(RESULTSDIR, 'flu_elmdict_' + org + ' ')
	# 	   + '> ' + os.path.join(RESULTSDIR, 'flu_elmdict_' + org + '.redo'))

@task
def redo_elmdict_realFrac():
	""" subtract counts expected by chance """
	
	for g in GENOMES:
		sh('python real_fraction.py '
		   + os.path.join(RESULTSDIR, g + '.aa_freq ')
		   + os.path.join(DATADIR, g + '.fa ')
		   + os.path.join(RESULTSDIR, 'elmdict_' + g + '.txt ')
		   + os.path.join(RESULTSDIR, g + '.diAA_freq ')
		   + '> ' + os.path.join(RESULTSDIR, 'elmdict_' + g + '.realFraction'))


@task
def simplify_flu_elms():
    """Subs flu residues for properties"""

    host_strains = [['human','H1N1'],
                    ['human','H3N2'],
                    ['human','H5N1'],
                    
                    ['swine','H3N2'],
                    ['swine','H1N1'],
                    
                    ['equine','H3N8'],
                    
                    ['chicken','H9N2'],
                    ['chicken','H5N1'],
                    
                    ['duck','H9N2'],
                    ['duck','H5N1']]
    for host, strain in host_strains:
        sh('python simplify_flu_elms.py '
           + 'results/' + host + '.' + strain + '.elms '
           + '> working/Jul1/' + host  + '.' + strain + '.simpleELMs')

@task
def simplify_flu_elms_2():
    """Subs flu residues for properties"""

    host_strains = [['human','H1N1'],
                    ['human','H3N2'],
                    ['human','H5N1'],
                    
                    ['swine','H3N2'],
                    ['swine','H1N1'],
                    
                    ['equine','H3N8'],
                    
                    ['chicken','H9N2'],
                    ['chicken','H5N1'],
                    
                    ['duck','H9N2'],
                    ['duck','H5N1']]
    for host, strain in host_strains:
        for year in xrange(2000, 2011):
            sh('python simplify_flu_elms.py '
               + 'working/Jul1_year/' + host + '.' + strain + '.' + str(year) + '.elms '
               + '> working/Jul1_year/' + host  + '.' + strain + '.' + str(year) + '.simpleELMs')

@task
@cmdopts([('cutoff=', 'c', '% cutoff'),])
def conserved_elms():
	"""Find ELMs conserved on strains"""

	cut = options.conserved_elms.get('cutoff')
	host_strains = [['human','H1N1'],
			['human','H3N2'],
			['human','H5N1'],

			['swine','H3N2'],
			['swine','H1N1'],

			['equine','H3N8'],
			
			['chicken','H9N2'],
			['chicken','H5N1'],

			['duck','H9N2'],
			['duck','H5N1']]

	for host, strain in host_strains:
		sh('python get_flu_seqs.py '
		   + host + ' '
		   + strain + ' NA results')
		sh('python matchELMpattern.py '
		   + 'elm_expressions.txt '
		   + 'results/' + host + '.' + strain + '.fa '
		   + '> ' + 'results/' + host + '.' + strain + '.elms')
		sh('python getConserved.py '
		   + 'results/' + host + '.' + strain + '.elms '
		   + 'ELM '
		   + cut + ' '
		   + '1> results/' + host + '.' + strain + '.elms.' + cut + ' '
		   + '2> results/' + host + '.' + strain + '.elms.conservation')

@task
def get_flu_seqs():
    """Parse NCBI flu file & make simple fasta"""

    host_strains = [['human','H1N1'],
                    ['human','H3N2'],
                    ['human','H5N1'],
                        
                    ['chicken','H9N2'],
                    ['chicken','H5N1'],
                    ['chicken','H7N2']]
                    
    for host, strain in host_strains:
        for year in xrange(2000, 2010):
            sh('python get_flu_seqs.py '
                + host + ' '
                + strain + ' '
                + str(year) + ' '
                + 'working/Jul22/fasta_use/')
            sh('python mk_simple_fasta.py '
               + 'working/Jul22/fasta_use/' + host + '.' + strain + '.' + str(year) + '.fa '
               + 'working/Jul22/' + host + '.' + strain + '.' + str(year) + '.fa')
    

@task
@cmdopts([('cutoff=', 'c', '% cutoff'),])
def conserved_elms_2():
	"""Find ELMs conserved on strains"""

	cut = options.conserved_elms_2.get('cutoff')
	host_strains = [#['human','H1N1'],
                        #['human','H3N2'],
                        #['human','H5N1'],
                        
                        #['swine','H3N2'],
                        #['swine','H1N1'],
                        
                        #['equine','H3N8'],
			
                        #['chicken','H9N2'],
                        #['chicken','H5N1'],
                        ['chicken','H7N2']]
                        
                        #['duck','H9N2'],
                        #['duck','H5N1']]

	for host, strain in host_strains:
            for year in xrange(2000, 2011):
		#sh('python get_flu_seqs.py '
		#   + host + ' '
		#   + strain + ' '
                #   + str(year) + ' '
                #   + 'working/Jul7/')
		sh('python matchELMpattern_once.py '
		   + 'working/Jul22/simple_patterns '
		   + 'working/Jul22/' + host + '.' + strain + '.' + str(year) + '.fa '
		   + '> ' + 'working/Jul22/' + host + '.' 
                   + strain + '.' + str(year) + '.elms_once')
		sh('python getConserved.py '
		   + 'working/Jul22/' + host + '.' + strain + '.' + str(year) + '.elms_once '
		   + cut + ' '
		   + '1> working/Jul22/' + host + '.' + strain + '.' + str(year) + '.elms.' + cut + ' '
		   + '2> working/Jul22/' + host + '.' + strain + '.' + str(year) + '.elms.conservation')

@task
def serotypes_random():
	""" run serotypes with random ELMs for human """
	
	#type2protein2gb2seq = utils.get_fluSeqs_by_serotype('human')
	for r in xrange(10):
		r_str = str(r)
		#sh('mkdir -p random/' + r_str)
		for t in ['H1N1', 'H5N1']:#, 'H3N2']:
			#sh('python matchELMpattern.py '
			#   + 'elm_exp_random' + r_str +  ' '
			#   + 'human.' + t + '.fa '
			#   + '> random/' + r_str + '/human.' + t + '.elms')
			sh('python getConserved.py '
			   + 'random/' + r_str + '/human.' + t + '.elms '
			   + 'ELM '
			   + '90 '
			   + '1> random/' + r_str + '/human.' + t + '.elms.90 '
			   + '2> random/' + r_str + '/human.' + t + '.elms.conservation')
			sh('python mk_freq.py '
			   + 'random/' + r_str + '/human.' + t + '.elms.90 '
			   + 'random/' + r_str + '/human.' + t + '.elms '
			   + '> random/' + r_str + '/human.' + t + '.elms.90.freq')

@task
def serotypes_random_fasta():
	""" run serotypes with random flu sequences for human """
	
	species = 'swine'
	#type2protein2gb2seq = utils.get_fluSeqs_by_serotype('human')
	
	for r in xrange(10):
		r_str = str(r)
		sh('mkdir -p random_seq/' + r_str)
		for t in ['H3N2','H1N1']:
			#utils.mk_random_fasta('results/' + species + '.' + t + '.fa',
			#		      'random_seq/' + r_str + '/' + species + '.' + t + '.fa')
			#sh('python matchELMpattern.py '
			#   + 'elm_expressions.txt '
			#   + 'random_seq/' + r_str + '/' + species + '.' + t + '.fa '
#			   + '> random_seq/' + r_str + '/' + species + '.' + t + '.elms')
			for cons in (70,80):
				c = str(cons)
				sh('python getConserved.py '
				   + 'random_seq/' + r_str + '/' + species + '.' + t + '.elms '
				   + 'ELM '
				   + str(c) + ' '
				   + '1> random_seq/' + r_str + '/' + species + '.' + t + '.elms.' + c + ' '
				   + '2> random_seq/' + r_str + '/' + species + '.' + t + '.elms.conservation')
				sh('python mk_freq.py '
				   + 'random_seq/' + r_str + '/' + species + '.' + t + '.elms.' + c + ' '
				   + 'random_seq/' + r_str + '/' + species + '.' + t + '.elms '
				   + '> random_seq/' + r_str + '/' + species + '.' + t + '.elms.' + c + '.freq')


