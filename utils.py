from local_settings import *
from global_settings import *
from Bio import SeqIO
import Bio.Cluster, Levenshtein
from copy import deepcopy
import re, os.path, os, logging, math, MySQLdb, local_settings
from functools import partial
from itertools import *
from collections import defaultdict
import markov_chain, utils_graph, math
from scipy.spatial import distance
import cloud, random, numpy, utils_motif

def init_zero(): return 0

def get_proteins_from_elm_file(afile):
     """Grab the proteins with ELM hits"""

     proteins = {}
     with open(afile) as f:
          for line in f:
               proteins[line.split('\t')[0]] = True
     return proteins

def count_cons(use_files, protein_counts_pass, f, d, new_f):
     protein_counts = defaultdict(dict)
     proteinName2motif = get_proteins_from_elm_file(f)
     for proteinName in proteinName2motif:
         protein = proteinName.split('.')[1]
         if protein in FLU_PROTEINS_LTD:
             protein_counts[protein][proteinName] = True
     for protein in protein_counts:
         if len(protein_counts[protein]) > 50:
             use_files[new_f] = True
             protein_counts_pass[protein][new_f] = True

def get_cons_elms(dir, hosts, years, strains, per, d, out_file, suffix):
    """Find ELMs that are consered at some per
       for all host/strain/year combinations w/
       at least 50 sequences"""

    d1 = {'ELM':True}
    d2 = d
    use_files = {}
    protein_counts_pass = defaultdict(dict)
    for host in hosts:
        for year in years:
            for strain in strains:
                f = os.path.join(dir, '.'.join((host, strain, str(year))) + '.elms')
                new_f = os.path.join(dir, '.'.join((host, strain, str(year))) + suffix + '.' + per)
                try:
                    count_cons(use_files, protein_counts_pass, f, d1, new_f)
                             #print host, year, strain
                except: pass
    for f in use_files:
        use_files[f] = utils_motif.protein2annotation(f, d2)
#    pass_elms = 

    with open(out_file, 'w') as afile:
         for protein in protein_counts_pass:
              #print protein + '\t' + str(len(protein_counts_pass[protein])) + '\t' + str([x.split('/')[2].split('.')[0:3] for x in protein_counts_pass[protein].keys()])
              elm_counts_local = defaultdict(init_zero)
              for f in protein_counts_pass[protein]:
                   if protein in use_files[f]:
                        for elm in use_files[f][protein]:
                             elm_counts_local[elm] += 1
              for elm in elm_counts_local:
                   if len(protein_counts_pass[protein]) == elm_counts_local[elm]:
                        afile.write(protein + '\t' + elm + '\n')
              # #      #else:
              # #      #     afile.write(protein + '\t' + elm + '\tFAIL\t' + str(elm_counts_local[elm]) + '\t' + elm + '\n')

def mk_sub(seq):
    """Make substitutions based on
       residue properties"""
    
    new_seq = ''.join([AA_SUB_2[c] 
                       for c in seq])
    return new_seq

def init_mysql(database):
    conn = MySQLdb.connect(host='localhost',
                           user=local_settings.MYSQL_USR,
                           passwd=local_settings.MYSQL_PASS,
                           db=database)
    cur = conn.cursor()
    return (conn, cur)

# def fasta_inter_mysql(database):
#     """iterator for fasta seqs in database"""

#     (conn, cur) = init_mysql('fasta')
#     line = 'SELECT * from ' + database
#     cur.execute(line)

def fasta_iter(afile, getID=None):
    """generator for FASTA files
       You can define your own header parser getID.
       
       Ex.
       D = dict(fasta_itr('empty', getID=lambda line: line.split('|')[1]))
       print ', '.join(D.keys())"""

    if not getID:
        getID = lambda line: line[1:].strip()
    isheader = lambda line: line[0] == '>'
    with open(afile) as f:
        for header, group in groupby(f, isheader):
            if header:
                line = group.next()
                ID = getID(line)
            else:
                seq = ''.join(line.strip() for line in group)
                yield ID, seq

def chunks(ls, n):
    """Yield n-sized chunks for ls"""
    for i in xrange(0, len(ls), n):
        yield ls[i:i+n]

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

def getDistFromCount(counts):
    s = float(sum(counts))
    if s == float(0):
        return [float(0) for x in counts]
    else:
        return [float(x)/s for x in counts]

def klDistance(dist1, dist2):
    """KL divergence between 2 distributions. No 0 counts."""

    dis = float(0)
    for p1, p2 in izip(dist1, dist2):
        if p1:
            dis += p1 * numpy.log(p1/p2)
    return dis

def jensen_shannon(counts1, counts2):
    """Find jensen shannon divergence between
       lists of counts.
       Make sure the lists count the same
       features"""

    # make counts into distributions
    d1 = getDistFromCount(counts1)
    d2 = getDistFromCount(counts2)
    davg = [numpy.average(x) for x in izip(d1,d2)]

    # find KL divergence from average
    return numpy.average([klDistance(d1, davg),
                          klDistance(d2, davg)])

def jensen_shannon_dists(d1, d2):
    """Find jensen shannon divergence between
       distributions.
       Make sure the lists count the same
       features.
       Use this version to avoid calculating
       the same distribution many times."""

    davg = [numpy.average(x) for x in izip(d1,d2)]

    # find KL divergence from average
    return numpy.average([klDistance(d1, davg),
                          klDistance(d2, davg)])

def test_jensen_shannon():
    c1=[10,20,30]
    c2=[15,10,25]
    print jensen_shannon(c1,c2)

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


		
def DictFromGen(GEN, elmfile, label = None, chunk_size = 10, stop_count = None):		
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

	elm_dict = ReadELMs_nocompile(elmfile)
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



def dict_init_zero(): return defaultdict(init_zero)

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
                
############### cluster ELM sequences #############
def fix_overlap(elm, mapping, overlap, new_clusters, dis_cutoff):
    """Assign sequences in overlap
       based on edit distance. Assume
       there are no ties to break."""

    for elmSeq in overlap:
        dis_cluster = []
        for cluster in new_clusters:
            dis = numpy.average([Levenshtein.distance(a_elmSeq.split(':')[1],
                                                      elmSeq.split(':')[1])
                                 for a_elmSeq in new_clusters[cluster]])
            dis_cluster.append([dis, cluster])
        dis_cluster.sort()
        best_cluster = dis_cluster[0]
        #print dis_cluster[0], dis_cluster[1]
        if best_cluster[0] < dis_cutoff:
            mapping[elm][elmSeq] = elm + ':' + str(best_cluster[1])
        else: # make new cluster
            mapping[elm][elmSeq] = elm + ':' + str(len(new_clusters)+1)
            new_clusters[len(new_clusters)+1][elmSeq] = True

def mk_mapping(elm, clusters, overlap, mapping, dis_cutoff):
    """Update a map of sequences to cluster.
       Ignore sequences in overlap that cannot
       be assigned."""

    new_clusters = defaultdict(dict)
    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                mapping[elm][seq] = elm + ':' + str(cluster)
                new_clusters[cluster][seq] = True
    fix_overlap(elm, mapping, overlap, new_clusters, dis_cutoff)

def get_initial_clusters(distance_file, dis_cutoff):
    """Make a cluster for each flu sequence.
       Place in the cluster any host sequence
       below some threshold (dis_cutoff)."""

    # map each flu elmSeq to a host elmSeq
    flu_host_elmSeq_mapping = {}

    with open(distance_file) as f:
        for line in f:
            (elm, flu_seq, host_seq, distance) = line.strip().split('\t')
            dis = int(distance)
            host_elmSeq = elm + ':' + host_seq
            flu_elmSeq = elm + ':' + flu_seq

            if dis < dis_cutoff:
                if elm not in flu_host_elmSeq_mapping:
                    flu_host_elmSeq_mapping[elm] = {}

                if flu_elmSeq not in flu_host_elmSeq_mapping[elm]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq] = {}

                if host_elmSeq not in  flu_host_elmSeq_mapping[elm][flu_elmSeq]:
                    flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

                flu_host_elmSeq_mapping[elm][flu_elmSeq][host_elmSeq] = True

    return flu_host_elmSeq_mapping

def get_meta_clusters(flu_host_elmSeq_mapping, dis_cutoff):
    """Cluster flu sequence clusters based on 
       overlapping host sequences."""

    mapping = defaultdict(dict)
    for elm in flu_host_elmSeq_mapping:
        dis_mat = []
        for seq1,seq2 in combinations(flu_host_elmSeq_mapping[elm], 2):
            match = float(len(set(seq1) & set(seq2)))
            dis = numpy.average([match/float(x) for x in 
                                 [len(seq1), len(seq2)]])
            dis_mat.append(dis)
        mat = numpy.array(dis_mat)
        num_clusters = 5
        percent_overlap = float(1)
        while num_clusters > 1 and percent_overlap > float(.1):
            if len(flu_host_elmSeq_mapping[elm]) > num_clusters:
                ans, error, nfound = Bio.Cluster.kmedoids(mat, 
                                                          nclusters=num_clusters, 
                                                          npass=20)
                clusters = defaultdict(dict)
                total = {}
                for flu_seq, cluster_id in zip(flu_host_elmSeq_mapping[elm],
                                               ans):
                    clusters[cluster_id][flu_seq] = True
                    total[flu_seq] = True
                    for host_seq in flu_host_elmSeq_mapping[elm][flu_seq]:
                        clusters[cluster_id][host_seq] = True  
                        total[host_seq] = True
                overlap = {}
                for cluster1, cluster2 in combinations(clusters, 2):
                    for gene in set(clusters[cluster1]) & set(clusters[cluster2]):
                        overlap[gene] = True
                percent_overlap = float(len(overlap))/float(len(total))
                if percent_overlap < float(.1):
                    #print_results(elm, clusters, overlap)
                    mk_mapping(elm, clusters, overlap, mapping, dis_cutoff)
                    break
                else:
                    num_clusters -= 1
            else:
                break
    return mapping

def get_clusters(distance_file, dis_cutoff_init, dis_cutoff_meta):
    """Return map of sequences to cluster name.
       dis_cutoff_init decides when to add to the
       initial clustering.
       dis_cutoff_meta decides the second clustering"""

    flu_host_elmSeq_mapping = get_initial_clusters(distance_file,
                                                   dis_cutoff_init)
    mapping = get_meta_clusters(flu_host_elmSeq_mapping, dis_cutoff_meta)
    return mapping

############# Make ELM sequence vectors for JS divergence

def count_flu(protein2counts, mapping, all_elmSeqs, do_clusters):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                if do_clusters:
                    elm = elmSeq.split(':')[0]
                    if elm in mapping:
                        if elmSeq in mapping[elm]:
                            key = mapping[elm][elmSeq]
                            counts[key] += protein2counts[protein][seq][elmSeq]
                            all_elmSeqs[key] = True
                else:
                    counts[elmSeq] += protein2counts[protein][seq][elmSeq]
                    all_elmSeqs[elmSeq] = True
    return counts

def get_flu_counts(afile, proteins):
    """Make protein_name -> seq_name -> elm_seq_counts"""

    counts = {}
    with open(afile) as f:
        for line in f:
            (protein, st, stp,
             elm, seq, junk) = line.strip().split('\t')
            name = protein.split('.')[-1]
            elmSeq = elm + ':' + seq
            if name in proteins:
                if name not in counts:
                    counts[name] = {}
                if protein not in counts[name]:
                    counts[name][protein] = {}
                if elmSeq not in counts:
                    counts[name][protein][elmSeq] = 0
                counts[name][protein][elmSeq] += 1
    return counts

def count_flu_sampled(flu, elm_file, flu_counts, seen_seqs, mapping, do_clusters):
    """Sample flu sequences for each flu protein
       so that all proteins are equally represented.
       Then do the ELM:seq counts.
       Accumulate results in flu_counts & seen_seqs"""

    pre = get_flu_counts(elm_file, 
                         FLU_PROTEINS)
    flu_counts_sampled = {}
    flu_proteins_sampled = {}
    protein_counts = []
    for protein in pre:
        protein_counts.append(len(pre[protein]))
    m = min(protein_counts)
    for protein in pre:
        if m == len(pre[protein]):
            flu_proteins_sampled[protein] = pre[protein].keys()
        else:
            flu_proteins_sampled[protein] = random.sample(pre[protein], m)
    for protein in flu_proteins_sampled:
        flu_counts_sampled[protein] = {}
        for sampled_protein in flu_proteins_sampled[protein]:
            flu_counts_sampled[protein][sampled_protein] = pre[protein][sampled_protein]
    seen_seqs[flu] = {}
    flu_counts[flu] = count_flu(flu_counts_sampled, 
                                mapping, seen_seqs[flu],
                                do_clusters)

def mk_vec(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for this host's counts"""
    
    vec = []
    for elmseq in all_elmSeqs:
        if elmseq in counts:
            vec.append(counts[elmseq])
        else:
            vec.append(float(0))
    return vec

def mk_count_vecs(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for all hosts"""

    vecs = {}
    for host in counts:
        vecs[host] = mk_vec(counts[host],
                            all_elmSeqs)
    return vecs

def mk_count_dists(vecs):
    """change count vectors into distributions"""

    dists = {}
    for host in vecs:
        dists[host] = getDistFromCount(vecs[host])
    return dists

def count_0s(ls):
    """How many 0s are in this vector?"""

    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

def print_it(name, vec):
    print name, float(count_0s(vec))/float(len(vec))

def print_results(elm, clusters, overlap):
    """Print out clustering of ELM sequences.
       Ignore sequences in overlap that cannot
       be assigned."""

    for cluster in clusters:
        for seq in clusters[cluster]:
            if seq not in overlap:
                print('%s\t%d\t%s' %
                      (elm, cluster, seq))

def count_host_elmSeqs(hosts, do_clustering, mapping, results_dir, use_elms, suffix):
    """Count elm:seq occurence.
       Choose to use mapping or not
       Only look at ELMs in use_elms
       
       Return host->ELMseq>count."""

    counts = {}
    for host in hosts:
        counts[host] = defaultdict(init_zero)
        elm_file = os.path.join(results_dir, 
                                'elmdict_' + host + suffix)
        with open(elm_file) as f:
            for line in f:
                (elm, seq, count, fq) = line.strip().split('\t')
                if elm in use_elms:
                    elmSeq = elm + ':' + seq
                    if do_clustering:
                        if elmSeq in mapping[elm]:
                            key = mapping[elm][elmSeq]
                            counts[host][key] += int(count)
                    else:
                        counts[host][elmSeq] += int(count)
    return counts
