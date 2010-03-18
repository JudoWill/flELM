from local_settings import *
from Bio import SeqIO
from copy import deepcopy
import re, os.path




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
	
	def ProcessLine(line):
		"""Processes the genomeset.dat line and returns a dict or properties"""
		out = {}
		parts = line.split('\t')
		out['genbank'] = parts[0]
		out['source'] = parts[1]
		out['strain'] = parts[3]
		out['country'] = parts[4]
		if parts[1] == 'Human': 
			out['organism'] = 'Human' 
		else:
			out['organism'] = parts[7].split('/')[1]
		return out
		
		
	
	vdict = {}
	always_true = lambda x: True
	easy = ('year','strain','source','country','organism')
	for key in easy:
		if key in kwargs:
			vdict[key] = lambda x: x == deepcopy(kwargs[key])
			kwargs.pop(key)
		else:
			vdcit[key] = always_true
	if 'genbank' in kwargs:
		vdict['genbank'] = lambda x: re.match(kwargs['genbank'],x)
		kwargs.pop('genbank')
	else:
		vdict['genbank'] = always_true
	
	if len(kwargs) > 0:
		raise KeyError, '%s are not valid kwargs' % ','.join(kwargs.keys())
	
	
	valid_seqs = set()
	
	with open(os.path.join(DATADIR, 'genomeset.dat')) as handle:
		for line in handle:
			odict = ProcessLine(line.strip())
			all_true = True
			for key, fun in vdict.items():
				if not fun(odict[key]):
					all_true = False
					break
			if all_true:
				valid_seqs.add(odict['genbank'])
	
	with open(os.path.join(DATADIR, 'influenza.faa')) as handle:
		for seqreq in SeqIO.parse(handle, 'fasta'):
			gb = seqreq.id.split('|')[3]
			if gb in valid_seqs:
				yield str(seqreq.seq)
				valid_seqs.remove(gb)
		
		
		
		
		
		
		