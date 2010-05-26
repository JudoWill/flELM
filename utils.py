from local_settings import *
from global_settings import *
from Bio import SeqIO
from copy import deepcopy
import re, os.path, os, logging, math
from functools import partial
from itertools import *
from collections import defaultdict
import markov_chain, utils_graph, math
from scipy.spatial import distance
from numpy import *

import cloud

def seqDistance(virus_d, host_d):
    seqs = utils_graph.intersectLists([virus_d,
                                       host_d])
    sum = float(0)
    for seq in seqs:
        sum += virus_d[seq]
    return sum

norm = lambda x : sqrt(square(x).sum())

def chiDistance(virus_d, host_d, means):
    """chi distance function
       http://www.econ.upf.edu/~michael/stanford/maeb4.pdf"""

    dis = float(0)
    for seq in means:
        x = float(0)
        y = float(0)
        if seq in host_d:
            x = host_d[seq]
        if seq in virus_d:
            y = virus_d[seq]
        dis += (x-y)**2/means[seq]
        print (x-y)**2, means[seq]
    return math.sqrt(dis)    

def renorm(seqs, d):
    new_d = {}
    total = float(0)
    for seq in seqs:
        if seq in d:
            total += d[seq]
    for seq in seqs:
        if seq in d:
            new_d[seq] = d[seq]/total
    return new_d

def klDistance(virus_d, host_d, use_seqs):
    d1 = float(0)
    seqs = utils_graph.intersectLists([virus_d,
                                       use_seqs])

    #print seqs
    virus_d_rn = renorm(seqs, virus_d)
    host_d_rn = renorm(seqs, host_d)
    for seq in virus_d_rn:
        p_x = virus_d_rn[seq]
        if seq in host_d_rn:
            q_x = host_d_rn[seq]
            d1 += p_x * log(p_x/q_x)
    if len(seqs.keys()) != 0:
        return d1
    else:
        return NaN

def getDistance(virus_d, host_d):
    seqs = utils_graph.unionLists([virus_d,
                                   host_d])
    #virus_d_rn = renorm(seqs, virus_d)
    #host_d_rn = renorm(seqs, host_d)

    host_v = []
    virus_v = []
    for seq in seqs:
        for v,d in ( (host_v, host_d),
                     (virus_v, virus_d) ):
            if seq in d:
                v.append(d[seq])
            else:
                v.append(float(0))

    # host_norm = norm(host_v)
    # virus_norm = norm(virus_v)
    # if host_norm:
    #     host_u = host_v/host_norm
    # else:
    #     host_u = host_v
    # if virus_norm:
    #     virus_u = virus_v/virus_norm
    # else:
    #     virus_u = virus_v
    dis = distance.cosine(host_v, virus_v)
    #print dis, virus_v, host_v
    return dis

def loadFASTA(fasta_file):
    """ Make a {} from a FASTA file.

    @param fasta_file: >name
                       seq
                       seq ...
    @return: fasta[gene] = seq (oneline)
    """

    fasta = {}
    fasta_f = open(fasta_file)
    name = ''
    seq = ''
    line = fasta_f.readline()
    while line != '':
        if line[0] == '>':
            if name != '':
                fasta[ name ] = seq
            name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = fasta_f.readline()
    fasta_f.close()

    if name != '':    
        fasta[name] = seq

    return fasta

def mk_random_fasta(fasta_file, out_fasta_file):
    """ given a flu FASTA file of ID.protein,
        use the markov chain to produce
        a random FASTA file.
        Keep the same # of protein types
    """

    fasta = loadFASTA(fasta_file)
    proteins = defaultdict(list)
    for k in fasta:
        protein = k.split('.')[1]
        proteins[protein].append(fasta[k])
    with open(out_fasta_file, 'w') as f:
        for protein in proteins:
            chain = markov_chain.MarkovChain()
            for seq in proteins[protein]:
                chain.add(seq)
            for i in xrange(len(proteins[protein])):
                f.write('>' + str(i) + '.' + protein + '\n')
                f.write("".join(chain.random_output()) + '\n')

def take(n, iterable):
	"Return first n items of the iterable as a list"
	return list(islice(iterable, n))

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
	handle = os.path.join(DATADIR, 'influenza.fa')
	fasta = loadFASTA(handle)
	for seqreq in fasta:
		gb = seqreq.split('|')[3]
		protein = seqreq.split('|')[4].split('[')[0].strip()
		if protein in PROTEIN_ALIAS:
			protein = PROTEIN_ALIAS[protein]
		if gb in valid_seqs:
			yield [str(fasta[seqreq]), protein, gb]
			valid_seqs.remove(gb)
		
def ReadELMs_nocompile(filename):
	"""Reads ELM file and returns dictionary of (name:regxp)"""

	outdict = dict()
	with open(filename) as handle:
		for line in handle:
			parts = line.strip().split()
			outdict[parts[0]] = parts[1]

	return outdict


		
def DictFromGen(GEN, label = None, chunk_size = 10, stop_count = None):		
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
		
	def ProcessSeq(seq_chunk, d):
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
					id = cloud.call(ProcessSeq, chunk[c-1::s], elm_dict, _label = label)
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
		
def get_seq2count_dict(elm_file, cutoff):
    """ grab ELM seq frequency data """

    elm2seq2count = defaultdict(dict)
    with open(elm_file) as f:
        for line in f:
            [elm, seq, count, frac_st] = line.strip().split('\t')
            frac = float(frac_st)
            if frac >= cutoff:
                elm2seq2count[elm][seq] = frac
    return elm2seq2count

def init_zero(): return 0

def get_seq2count_dict_for_seqs(elm_file, cutoff, virus2elms2seq):
    """ grab ELM seq frequency data 
        ONLY for these ELM sequences
	renormalize for these ELM sequences
    """

    elm2seq2count = defaultdict(dict)
    elm_totals = defaultdict(init_zero)
    with open(elm_file) as f:
	    for line in f:
		    [elm, seq, count_st, frac_st] = line.strip().split('\t')
		    for virus in virus2elms2seq:
			    if elm in virus2elms2seq[virus]:
				    if seq in virus2elms2seq[virus][elm]:
					    count = float(count_st)
					    elm_totals[elm] += count
					    elm2seq2count[elm][seq] = count
				    
    # renormalize
    to_remove = {}
    for elm in elm2seq2count:
	    for seq in elm2seq2count[elm]:
		    frac = float(elm2seq2count[elm][seq])/float(elm_totals[elm])
		    if frac >= cutoff:
			    elm2seq2count[elm][seq] = frac
		    else: # ignore this fraction
			    to_remove[elm+':'+seq] = True
    for k in to_remove:
	    [elm, seq] = k.split(':')
	    del elm2seq2count[elm][seq]
    return elm2seq2count

def check_ones(species2elms, elm):
    """ Is there only one sequence for this ELM
        in all species? 
    """

    not_one = False

    for species in species2elms:
        if len(species2elms[species][elm].keys()) > 1:
            not_one = True
            break
    return not_one

def calc_entropy(prob_list):
    """ given a list of probabilities,
	find the entropy
    """
    
    entropy = float(0)
    for prob in prob_list:
        entropy -= prob * math.log(prob, 2)
    return entropy

def get_species_entropy(elm2seq2prob):
    """ find entropy for all ELMs
        given {} of ELM 2 seq
    """

    entropy = {}
    for elm in elm2seq2prob:
        prob_ls = []
	for seq in elm2seq2prob[elm]:
            prob_ls.append(elm2seq2prob[elm][seq])
	entropy[elm] = calc_entropy(prob_ls)
    return entropy

def calc_elm_frequency(elmdict_file):
	""" find the fraction of ELM occurance
            out of total ELMs """
	counts = defaultdict(init_zero)
	total = 0
	with open(elmdict_file) as f:
		for line in f:
			[elm, seq, count, frac] = line.strip().split('\t')
			counts[elm] += int(count)
			total += int(count)
	for elm in counts:
		counts[elm] = float(counts[elm])/float(total)
	return counts

def get_fluSeqs_by_serotype(flu_shortname):
    """ return {} of serotype 2 proteinName 2 seq """

    gb2type = {}
    with open(os.path.join(DATADIR, 'influenza_aa.dat')) as f:
        for line in f:
            for shortname in FLU_NAMES[flu_shortname]:
                if shortname in line.lower() and 'Influenza A virus' in line:
                    if shortname == 'human' and line.split('\t')[1] == 'Avian':
                        pass
                    else:
                        gb = line.split('\t')[0]
                        serotype = line.split('\t')[3]
                        gb2type[gb] = serotype
    fasta = loadFASTA(os.path.join(DATADIR, 'influenza.fa'))
    type2protein2gb2seq = {}
    for key in fasta:
        gb = key.split('|')[3]
        protein = key.split('|')[4].split('[')[0].strip()
        if protein in PROTEIN_ALIAS:
            protein = PROTEIN_ALIAS[protein]
        if gb in gb2type:
            serotype = gb2type[gb]
            if not serotype in type2protein2gb2seq:
                type2protein2gb2seq[serotype] = {}
            if not protein in type2protein2gb2seq[serotype]:
                type2protein2gb2seq[serotype][protein] = {}
            type2protein2gb2seq[serotype][protein][gb] = fasta[key]
    return type2protein2gb2seq
                
