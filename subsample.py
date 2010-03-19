from local_settings import *
from global_settings import *
from utils import DictFromGen, GetFluSeqs
from makeELMdict import SeqGen
from optparse import OptionParser

import os.path, os, random, logging


import cloud


def SampleIterator(gen, fract_return):
	for it in gen:
		if random.random() < fract_return:
			yield it


if __name__ == '__main__':
	logging.basicConfig(level=logging.WARNING)
	parser = OptionParser()
	
	parser.add_option('-c', '--use-cloud', default = False, action = 'store_true',
						dest = 'usecloud', help = 'Use PiCloud computing')
	parser.add_option('-n', '--num-samples', default = 3000, type = 'int',
						help = 'Number of subsamples to create', dest = 'numsamples')
	parser.add_option('-p', '--percentage', default = 0.1, type = 'float',
						help = 'Percentage to sample', dest = 'percentage')
	parser.add_option('-f', '--forcenew', default = False, action = 'store_true',
						help = 'Force overwritting of previously generated samples', 
						dest = 'forcenew')
						
	(options, args) = parser.parse_args()
	
	cloud.setkey(CLOUD_KEY, CLOUD_SECRET)
	if not options.usecloud: cloud.start_simulator()
	
	revAliases = {}
	for key, val in ALIASES.items():
		revAliases[val.lower()] = key
	
	print ALIASES	
	
	for genome in args:
		
		short_name = ALIASES[genome]
		outdirect = os.path.join(RESULTSDIR, 'sampling', short_name)
		if not os.path.exists(outdirect):
			os.mkdir(outdirect)
		ifile = os.path.join(DATADIR, genome + '.fa')
		for i in xrange(options.numsamples):
			outfile = os.path.join(outdirect, short_name + '%05d.txt' % (i+1,))
			if os.path.exists(outfile) and not options.forcenew:
				continue
			
			logging.warning('Processing subsample %d' % (i+1,))
			
			fgen = SeqGen(ifile)
			sub_gen = SampleIterator(fgen, options.percentage)
			d = DictFromGen(sub_gen, label = short_name, chunk_size = 200)
			
			
			logging .warning('Writting file: %s' % outfile)
			with open(outfile, 'w') as handle:
				l = d.keys()
				l.sort()
				for key in l:
					count, frac = d[key]
					elm, spec = key
					handle.write('\t'.join(map(str, [elm, spec, count, '%.10f' % frac]))+'\n')
				
				
	
	
	
	
	
	
	
	