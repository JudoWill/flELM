from local_settings import *
from global_settings import *
from itertools import *
from utils import GetFluSeqs
from optparse import OptionParser
from collections import defaultdict

import os, os.path, re

import cloud




def take(n, iterable):
	"Return first n items of the iterable as a list"
	return list(islice(iterable, n))


def ReadELMs_nocompile(filename):
	"""Reads ELM file and returns dictionary of (name:regxp)"""

	outdict = dict()
	with open(filename) as handle:
		for line in handle:
			parts = line.strip().split()
			outdict[parts[0]] = parts[1]

	return outdict
	

def ProcessFlu(genome_tup):
	
	
	
	def ChunkGen(genome_tup, per_block):
		gen = GetFluSeqs(organism = genome_tup)
		block = take(per_block, gen)
		c = per_block
		while block:
			#print 'yeilding block %d, len: %d' % (c, len(block))
			#raise KeyError
			yield block
			block = take(per_block, gen)
			c += per_block
			#if c > 100: break
	
	def ProcessSeq(inp):
		seq_chunk, d = inp
		count_dict = defaultdict(int)
		elm_dict = {}
		for key, val in d.items():
			elm_dict[key] = re.compile(val)
		for i, seq in enumerate(seq_chunk):		
			#print i
			for elm, reg in elm_dict.items():
				m = reg.search(seq)
				while m:
					count_dict[(elm, m.group())] += 1
					m = reg.search(seq, m.start()+1)
		return count_dict
		
	label = genome_tup[0]
	elm_dict = ReadELMs_nocompile('elm_expressions.txt')

	count_dict = defaultdict(int)
	elm_count = defaultdict(int)
	jids = []
	for block in ChunkGen(genome_tup, 500):
		try:
			jids.append(cloud.call(ProcessSeq, (block, elm_dict), _label=label))
		except cloud.cloud.CloudException:
			try:
				print 'tooooo big'
				jids.append(cloud.call(ProcessSeq, (block[0::2], elm_dict), _label=filename))
				jids.append(cloud.call(ProcessSeq, (block[1::2], elm_dict), _label=filename))
			except cloud.cloud.CloudException:
				print 'really tooo big'
				jids.append(cloud.call(ProcessSeq, (block[0::3], elm_dict), _label=filename))
				jids.append(cloud.call(ProcessSeq, (block[1::3], elm_dict), _label=filename))
				jids.append(cloud.call(ProcessSeq, (block[2::3], elm_dict), _label=filename))

	#print str(jids)
	print 'waiting!'
	for i, res in enumerate(cloud.iresult(jids)):
		if i % 4 == 0: print 'processing result %d' % i
		for key, item in res.iteritems():
			elm, spec = key
			count_dict[key] += item
			elm_count[elm] += item

	cloud.delete(jids)

	outdict = {}
	for key, count in count_dict.iteritems():
		elm, spec = key
		outdict[key] = (count, float(count) / float(elm_count[elm]))

	return outdict	
	
if __name__ == '__main__':
	
	parser = OptionParser()
	
	parser.add_option('-c', '--use-cloud', default = False, action = 'store_true',
						dest = 'usecloud', help = 'Use PiCloud computing')
						
	(options, args) = parser.parse_args()
	
	cloud.setkey(CLOUD_KEY, CLOUD_SECRET)
	if not options.usecloud: cloud.start_simulator()
	
	for genome_tup in FLU_NAMES:
		outdata = ProcessFlu(genome_tup)
		print 'writting data'
		counter = 0
		with open(os.path.join(RESULTSDIR, 'flu_elmdict_'+genome_tup[0]), 'w') as handle:
			l = outdata.keys()
			l.sort()
			for key in l:
				count, frac = outdata[key]
				elm, spec = key

				#print key, val
				counter += 1
				#print key[0]
				#if counter > 100: break
				handle.write('\t'.join(map(str, [elm, spec, count, '%.10f' % frac]))+'\n')
	
	
	
	
	