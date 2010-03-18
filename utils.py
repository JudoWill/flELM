from local_settings import *
from Bio import SeqIO
from copy import deepcopy
import re, os.path, os
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
	print 'valid gbs:', len(valid_seqs)
	#print valid_seqs
	with open(os.path.join(DATADIR, 'influenza.fa')) as handle:
		for seqreq in SeqIO.parse(handle, 'fasta'):
			gb = seqreq.id.split('|')[3]
			#print gb
			if gb in valid_seqs:
				yield str(seqreq.seq)
				valid_seqs.remove(gb)
				#print 'yielding:', gb
		
		
		
		
		
		
		