from local_settings import *
from Bio import SeqIO
from copy import deepcopy
import re, os.path, os, logging
from functools import partial



def GetFluSeqs(*args, **kwargs):
	"""A generator which returns flu sequences based on provided requests
	
	yields: sequence
	
	Valid KWARGS:
	
	year (int)	1980
	strain (string)	H1N1
	source (string) Avian
	country (string) Hong Kong
	organism (string) Black Duck
	genbank (string/regexp) CY\d* NOT COMPILED!!!
	
	Any other kwargs raise a KeyError!
	
	"""
	
		
	def TruthTest(inp, line):
		for i in inp:
			if i in line:
				return True
		return False

	vdict = {}
	always_true = lambda x: True
	easy = ('year','strain','source','country','organism')
	for key in easy:
		if key in kwargs:
			vdict[key] = partial(TruthTest, deepcopy(kwargs[deepcopy(key)]))
			kwargs.pop(key)
		else:
			vdict[key] = always_true
	if 'genbank' in kwargs:
		vdict['genbank'] = lambda x: re.match(kwargs['genbank'],x)
		kwargs.pop('genbank')
	else:
		vdict['genbank'] = always_true
	
	if len(kwargs) > 0:
		raise KeyError, '%s are not valid kwargs' % ','.join(kwargs.keys())
	
	
	valid_seqs = set()
	
	logging.info('Getting Information')
	with open(os.path.join(DATADIR, 'influenza_aa.dat')) as handle:
		for line in handle:
			l_line = line.lower()
			bad = False
			for fun in vdict.values():
				if not fun(l_line):
					bad = True
					break
			if not bad:
				valid_seqs.add(line.split('\t')[0])
	
	logging.warning('valid gbs: %d' % len(valid_seqs))
	logging.info('yielding sequences')
	with open(os.path.join(DATADIR, 'influenza.fa')) as handle:
		for seqreq in SeqIO.parse(handle, 'fasta'):
			gb = seqreq.id.split('|')[3]
			if gb in valid_seqs:
				yield str(seqreq.seq)
				valid_seqs.remove(gb)
		
def ReadELMs_nocompile(filename):
	"""Reads ELM file and returns dictionary of (name:regxp)"""

	outdict = dict()
	with open(filename) as handle:
		for line in handle:
			parts = line.strip().split()
			outdict[parts[0]] = parts[1]

	return outdict


		
def DictFromGen(GEN, label = None, chunk_size = 500, stop_count = None):		
	"""Creates a dictionary of ELM frequencies from a generator"""
	
	
	
	def ChunkGen(gen, per_block, stop_count):
		block = take(per_block, gen)
		c = per_block
		while block:
			logging.debug('yeilding block %d, len: %d' % (c, len(block)))
			yield block
			block = take(per_block, gen)
			c += per_block
			if stop_count != None and stop_count < c: break
		
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

	elm_dict = ReadELMs_nocompile('elm_expressions.txt')
	jids = []
	
	for chunk in ChunkGen(GEN, chunk_size, stop_count):
		submitted = False
		s = 1
		tids = []
		while not submitted:
			c = 1
			while c<=s:
				try:
					#try to submit the current slice of the block
					logging.debug('Submitting a chunk to the cloud')
					id = cloud.call(ProcessSeq, chunk[c-1::s], _label = label)
					tids.append(id)
					c+=1
					submitted = True
				except cloud.cloud.CloudException:
					#if there is an exception because there is too much info
					#then kill, delete the cloud.call and then increase the
					#slicing
					submitted = False
					cloud.kill(tids)
					cloud.delete(tids)
					logging.warning('Chunk was too big at slice: %d' %c)
					s += 1
					break
		jids += tids
	
	count_dict = defaultdict(int)
	elm_count = defaultdict(int)
	logging.warning('Waiting for the cloud')
	for i, res in enumerate(cloud.iresult(jids)):
		logging.debug('Processing result: %s' % i)
		for key, item in res.iteritems():
			elm, spec = key
			count_dict[key] += item
			elm_count[elm] += item
		
	logging.info('Deleting jobs')
	cloud.delete(jids)

	logging.info('Creating output dictionary')
	outdict = {}
	for key, count in count_dict.iteritems():
		elm, spec = key
		outdict[key] = (count, float(count) / float(elm_count[elm]))

	return outdict		
		