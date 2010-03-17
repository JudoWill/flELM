import re
from collections import defaultdict
from Bio import SeqIO
from optparse import OptionParser
from itertools import islice, repeat, izip

import cloud

from local_settings import CLOUD_KEY, CLOUD_SECRET

def take(n, iterable):
	"Return first n items of the iterable as a list"
	return list(islice(iterable, n))


def SeqGen(filename):
	"""A generator which returns the raw sequence from a fasta file"""
	
	with open(filename) as handle:
		for seq_rec in SeqIO.parse(handle, 'fasta'):
			yield str(seq_rec.seq)

def ReadELMs(filename):
	"""Reads ELM file and returns dictionary of (name:regxp)"""
	
	outdict = dict()
	with open(filename) as handle:
		for line in handle:
			parts = line.strip().split()
			outdict[parts[0]] = re.compile(parts[1])
			
	return outdict

def ReadELMs_nocompile(filename):
	"""Reads ELM file and returns dictionary of (name:regxp)"""

	outdict = dict()
	with open(filename) as handle:
		for line in handle:
			parts = line.strip().split()
			outdict[parts[0]] = parts[1]

	return outdict


			
def ProcessFile(filename):
	"""Finds the frequencey of ELM instances in a FASTA file.
	
	Returns a defaultdict key-ed by ('ELM_NAME', 'ELM_SEQ') and has a value of
	counts/AA
	"""
	
	elm_dict = ReadELMs('elm_expressions.txt')
	
	count_dict = defaultdict(int)
	
	len_counter = 0
	for i, seq in enumerate(SeqGen(filename)):
		if i % 1000 == 0: print 'processing seq: %d' % i
		#if i > 1000: break
		len_counter += len(seq)
		for elm, reg in elm_dict.items():
			m = reg.search(seq)
			while m:
				count_dict[(elm, m.group())] += 1
				m = reg.search(seq, m.start()+1)
	
	#genome_total = float(len_counter)
	#for key, val in count_dict.iteritems():
	#	count_dict[key] = float(val) // genome_total
	
	return count_dict
	
def ProcessFile_withCloud(filename):
	"""Finds the frequencey of ELM instances in a FASTA file.

	Returns a defaultdict key-ed by ('ELM_NAME', 'ELM_SEQ') and has a value of
	counts/AA
	"""

	def ProcessSeq(inp):
		seq_chunk, d = inp
		count_dict = defaultdict(int)
		elm_dict = {}
		for key, val in d.items():
			elm_dict[key] = re.compile(val)
		for i, seq in enumerate(seq_chunk):		
			print i
			for elm, reg in elm_dict.items():
				m = reg.search(seq)
				while m:
					count_dict[(elm, m.group())] += 1
					m = reg.search(seq, m.start()+1)
		return count_dict

	def ChunkGen(filename, per_block):
		gen = SeqGen(filename)
		block = take(per_block, gen)
		c = per_block
		while block:
			#print 'yeilding block %d' % c
			yield block
			block = take(per_block, gen)
			c += per_block
			#if c > 5000: break

	elm_dict = ReadELMs_nocompile('elm_expressions.txt')

	count_dict = defaultdict(int)

	jids = cloud.map(ProcessSeq, izip(ChunkGen(filename, 100), repeat(elm_dict)))
	#print str(jids)
	print 'waiting!'
	for i, res in enumerate(cloud.iresult(jids)):
		print 'processing result %d' % i
		for key, item in res.iteritems():
			count_dict[key] += item
	
	
	return count_dict



if __name__ == '__main__':
	
	parser = OptionParser()
	
	parser.add_option('-o', '--outfile', dest = 'outfile',
						help = 'The destination file too write data')
	parser.add_option('-c', '--use-cloud', default = False, action = 'store_true',
						dest = 'usecloud', help = 'Use PiCloud computing')
						
	(options, args) = parser.parse_args()
	
	cloud.setkey(CLOUD_KEY, CLOUD_SECRET)
	if not options.usecloud: cloud.start_simulator()
	
	outdata = ProcessFile_withCloud(args[0])
	print 'writting data'
	counter = 0
	with open(options.outfile, 'w') as handle:
		l = outdata.keys()
		l.sort()
		for key in l:
			val = outdata[key]
			counter += 1
			#print key[0]
			#if counter > 100: break
			handle.write('%s\t%s\t%d\n' % (key[0], key[1], val))
	
	
	
	
	
	
	