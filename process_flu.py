from local_settings import *
from global_settings import *
from itertools import *
from utils import GetFluSeqs, DictFromGen
from optparse import OptionParser
from collections import defaultdict

import os, os.path, re, logging

import cloud


def take(n, iterable):
	"Return first n items of the iterable as a list"
	return list(islice(iterable, n))

	
if __name__ == '__main__':
	
	parser = OptionParser()
	
	parser.add_option('-c', '--use-cloud', default = False, action = 'store_true',
						dest = 'usecloud', help = 'Use PiCloud computing')
	
						
	(options, args) = parser.parse_args()
	
	cloud.setkey(CLOUD_KEY, CLOUD_SECRET)
	if not options.usecloud: cloud.start_simulator()
	
	for short_name in args:
		logging.warning('Processing Genome: %s' % short_name)
		flu_gen = GetFluSeqs(organism = FLU_NAMES[short_name])
		outdata = DictFromGen(flu_gen, label = short_name)
	
		logging.warning('writting data for %s' % short_name)
		with open(os.path.join(RESULTSDIR, 'flu_elmdict_'+short_name), 'w') as handle:
			l = outdata.keys()
			l.sort()
			for key in l:
				count, frac = outdata[key]
				elm, spec = key
				handle.write('\t'.join(map(str, [elm, spec, count, '%.10f' % frac]))+'\n')
	
	
	
	
	