import re
from collections import defaultdict
from Bio import SeqIO
from optparse import OptionParser


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
	
	
if __name__ == '__main__':
	
	parser = OptionParser()
	
	parser.add_option('-o', '--outfile', dest = 'outfile',
						help = 'The destination file too write data')
						
	(options, args) = parser.parse_args()
	
	
	outdata = ProcessFile(args[0])
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
	
	
	
	
	
	
	