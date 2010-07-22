"""Split ELM/seqs into those on all strains,
   and those specific to 90% mammal but not
   in 90% bird."""
import os, utils, global_settings, sys, utils_stats
from collections import defaultdict

outfile = sys.argv[1]

def write_latex(g1, g2, outfile):
    """Make table for paper"""

    with open(outfile, 'w') as f:
        f.write('\\par\n')
        f.write('\\mbox{\n')
        f.write('\\begin{tabular}{|c|c|c||c|c|}\n')
        f.write('\\hline & \\multicolumn{2}{|c||}{Mammal} & \\multicolumn{2}{|c|}{Bird} \\\\ \\hline\n')
        f.write('Hypoth & True  &  False & True & False \\\\ \\hline\n')
        f.write('Unique & %d  & %d & %d & %d \\\\ \\hline\n' %
                (g1[0], g1[1], g2[0], g2[1]))
        f.write(' Control & %d & %d & %d & %d \\\\ \\hline\n' % 
                (g1[2], g1[3], g2[2], g2[3]))
        f.write('\\hline P-value & \\multicolumn{2}{|c||}{%.4f} & \\multicolumn{2}{|c|}{%.4f} \\\\ \\hline\n' %
                (g1[-1], g2[-1]))
        f.write('\\end{tabular}}\n')

def get_seqs(uniq_file, use_protein):
    """Get uniq seqs on flu"""

    seqs = {}
    with open(uniq_file) as f:
        for line in f:
            protein, elmseq, elm, seq = line.strip().split('\t')
            seqs[seq] = True

    rm_seqs = {}
    for seq in seqs:
        for seq2 in seqs:
            if seq != seq2:
                if seq in seq2:
                    rm_seqs[seq] = True
    print len(seqs), len(rm_seqs)
    #for seq in rm_seqs:
    #    del seqs[seq]
    print len(seqs)
    return seqs

def get_freqs(file, elmdict_file):
    """Get host freqs"""

    counts = {}
    with open(elmdict_file) as f:
        for line in f:
            elmseq, seq, count, freq = line.strip().split('\t')
            counts[seq] = int(count)

    freqs = {}
    with open(file) as f:
        for line in f:
            seq, freq = line.strip().split('\t')
            #seq = elmseq.split(':')[1]
            if counts[seq] > 500:
                freqs[seq] = float(freq)
    return freqs

def count_evaluate(name, flu_seqs, host1_freqs, host2_freqs):
    """Test hypothesis"""

    host1_present = 0
    host2_present = 0
    for seq in flu_seqs:
        if seq in host1_freqs:
            host1_present += 1
        if seq in host2_freqs:
            host2_present += 1
    print name, host1_present, host2_present, len(flu_seqs)

def evaluate(name, flu_seqs, host1_freqs, host2_freqs):
    """Test hypothesis"""

    count = 0
    total = 0
    k = {}
    for f in flu_seqs:
        k[f] = True
    print 'TESTING', len(flu_seqs), len(k)
    for seq in flu_seqs:
        f1 = float(0)
        f2 = float(0)
        if seq in host1_freqs:
            f1 = host1_freqs[seq]
        if seq in host2_freqs:
            f2 = host2_freqs[seq]
        if f1 != float(0) and f2 != float(0):
            if f1 > f2:
                count += 1
            total += 1
        # if seq in host1_freqs and seq in host2_freqs:
        #     if host1_freqs[seq] > host2_freqs[seq]:
        #         count += 1
        #     total += 1
    print name, count, total-count, total, float(count)/float(total)
    return (count, total-count)

def get_protein_counts(dir, hosts, strains, years):
    """Count seqs you have per protein"""

    counts = {}
    for host in hosts:
        for strain in strains:
            for year in years:
                key = '.'.join((host, strain, str(year)))
                f = os.path.join(dir, key + '.fa')
                if os.path.exists(f):
                    counts[key] = defaultdict(dict)
                    for ID, seq in utils.fasta_iter(f):
                        protein = ID.split('.')[-1]
                        counts[key][protein][ID] = True
    return counts

def dump_it(dir, protein_elm, name):
    """Print annotation results"""

    with open(os.path.join(dir, name), 'w') as f:
        for protein in protein_elm:
            for seq in protein_elm[protein]:
                #elm, seq = elmseq.split(':')
                f.write(protein + '\t' + seq + '\t'
                        + seq + '\t' + seq + '\n')

def dump_seqs(dir, protein_elm, name):
    """Print annotation results. Just seqs"""

    seqs = {}
    for protein in protein_elm:
        for seq in protein_elm[protein]:
            #elm, seq = elmseq.split(':')
            seqs[seq] = True
    with open(os.path.join(dir, name), 'w') as f:
        for seq in seqs:
            f.write(seq + '\n')

def get_cons_elms(dir, protein_counts, limit):
    """Get all ELM conservation for proteins w/ over LIMIT seqs"""

    cons = {}
    for key in protein_counts:
        for protein in protein_counts[key]:
            if len(protein_counts[key][protein]) > limit:
                cons[key] = defaultdict(dict)
                f = os.path.join(dir, key + '.elms.conservation')
                with open(f) as elms:
                    for line in elms:
                        protein, elm, cons_st = line.strip().split('\t')
                        cons[key][protein][elm] = float(cons_st)
    return cons

def get_all_cons(cons, protein_counts, seq_limit, cons_limit):
    """Get protein/elm pairs that are annotated or all flu strains"""

    protein_elm_count = defaultdict(utils.init_zero)
    protein_count = defaultdict(utils.init_zero)
    for key in protein_counts:
        for protein in protein_counts[key]:
            if len(protein_counts[key][protein]) > seq_limit:
                protein_count[protein] += 1
                if protein in cons[key]:
                    for elm in cons[key][protein]:
                        if cons[key][protein][elm] > cons_limit:
                            protein_elm_count[protein + '=' + elm] += 1
    protein_elm_pass = defaultdict(dict)
    for protein_elm in protein_elm_count:
        protein, elm = protein_elm.split('=')
        if protein_elm_count[protein_elm] == protein_count[protein]:
            protein_elm_pass[protein][elm] = True

    return protein_elm_pass

def get_uniq(all_cons, cmp_cons):
    """What protein/ELM pairs are in all_cons, but are not in cmp_cons?
       Ex. mammal is all_cons, and bird is cmp_cons

       This leaves out proteins that are not in both sets. """

    uniq = defaultdict(dict)
    for protein in all_cons:
        if protein in cmp_cons:
            for elm in all_cons[protein]:
                if elm not in cmp_cons[protein]:
                    uniq[protein][elm] = True
        #else:
        #    print 'ignore', protein
    return uniq

def count_it(uniq):
    count = 0
    for protein in uniq:
        for elm in uniq[protein]:
            count += 1
    return count

def mk_control(species_all_cons, species_uniq):
    """Make a control using ELMs that are on 90% of this flu,
       but not unique and not conserved on all bird flus."""

    control = defaultdict(dict)
    for protein in species_all_cons:
        if protein in species_uniq:
            for elm in species_all_cons[protein]:
                if elm not in species_uniq[protein]:
                    control[protein][elm] = True
    return control    

dir = 'working/Jul20/'
years = range(2000,2011,1)
mammal_hosts = ('human',)#'swine','horse'
mammal_strains = ('H5N1','H1N1','H3N2')#,'H3N8','H1N1'
bird_hosts = ('chicken',)
bird_strains = ('H5N1','H9N2')#'H9N2'

limit = global_settings.SEQ_LIMIT
cons_cut = float(90)
low_cons_cut = float(60)
mammal_protein_counts = get_protein_counts(dir, mammal_hosts,
                                           mammal_strains, years)
bird_protein_counts = get_protein_counts(dir, bird_hosts,
                                         bird_strains, years)
mammal_cons = get_cons_elms(dir, mammal_protein_counts, limit)
bird_cons = get_cons_elms(dir, bird_protein_counts, limit)
mammal_all_cons = get_all_cons(mammal_cons, mammal_protein_counts, 
                               limit, cons_cut)
bird_all_cons = get_all_cons(bird_cons, bird_protein_counts, limit, cons_cut)
mammal_uniq = get_uniq(mammal_all_cons, bird_all_cons)
bird_uniq = get_uniq(bird_all_cons, mammal_all_cons)

# find what is conserved on everything
mammal_protein_counts.update(bird_protein_counts)
mammal_cons.update(bird_cons)
all_cons = get_all_cons(mammal_cons, mammal_protein_counts, limit, cons_cut)

print 'MAMMAL', count_it(mammal_uniq)
print 'BIRD', count_it(bird_uniq)
print 'ALL', count_it(all_cons)
print 'MAMMAL ALL', count_it(mammal_all_cons)
print 'BIRD ALL', count_it(bird_all_cons)

dump_seqs(dir, mammal_uniq, 'mammal_uniq_venn')
dump_seqs(dir, bird_uniq, 'bird_uniq_venn')

# what is not considered uniq?
mammal_control = mk_control(mammal_all_cons, mammal_uniq)
bird_control = mk_control(bird_all_cons, bird_uniq)

dump_seqs(dir, mammal_control, 'mammal_control_venn')
dump_seqs(dir, bird_control, 'bird_control_venn')

print 'MAMMAL CONTROL', count_it(mammal_control)
print 'BIRD CONTROL', count_it(bird_control)

dump_it(dir, mammal_uniq, 'mammal_uniq')
dump_it(dir, bird_uniq, 'bird_uniq')
dump_it(dir, mammal_all_cons, 'control')
dump_it(dir, mammal_control, 'mammal_control')
dump_it(dir, bird_control, 'bird_control')

protein = ''

mammal_pre = get_seqs(os.path.join(dir, 'mammal_uniq'),
                      protein)
bird_pre = get_seqs(os.path.join(dir, 'bird_uniq'),
                    protein)
mammal_control_pre = get_seqs(os.path.join(dir, 'mammal_control'),
                              protein)
bird_control_pre = get_seqs(os.path.join(dir, 'bird_control'),
                            protein)
mammal_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.simple.elm_aa_freq'),
                              os.path.join(dir, 'elmdict_H_sapiens.simple'))
bird_host_freqs = get_freqs(os.path.join(dir, 
                                         'Gallus_gallus.simple.elm_aa_freq'),
                            os.path.join(dir, 'elmdict_Gallus_gallus.simple'))

mammal = set(mammal_pre.keys()) - set(bird_pre.keys())
bird = set(bird_pre.keys()) - set(mammal_pre.keys())

mammal_control_pre2 = set(mammal_control_pre.keys()) - set(mammal_pre.keys())
bird_control_pre2 = set(bird_control_pre.keys()) - set(bird_pre.keys())
both_control = bird_control_pre2 & mammal_control_pre2
mammal_control = mammal_control_pre2 - both_control
bird_control = bird_control_pre2 - both_control

print 'intr', len(mammal_control & mammal)
print 'intr', len(bird_control & bird)
print 'intr', len(mammal & bird)
print 'intr', len(mammal_control & bird_control)

m1, m2 = evaluate('MAMMAL', mammal, mammal_host_freqs, bird_host_freqs)
c1, c2 = evaluate('MAMMAL CONTROL', mammal_control, mammal_host_freqs, 
                  bird_host_freqs)
mammal_pval = utils_stats.fisher_positive_pval((m1, m2),
                                               (c1, c2))
print mammal_pval

b1, b2 = evaluate('BIRD', bird, 
                  bird_host_freqs, mammal_host_freqs)
bc1, bc2 = evaluate('BIRD CONTROL', bird_control, 
                    bird_host_freqs, mammal_host_freqs)
bird_pval = utils_stats.fisher_positive_pval((b1, b2),
                                             (bc1, bc2))
print bird_pval

write_latex((m1, m2, c1, c2, mammal_pval),
            (b1, b2, bc1, bc2, bird_pval),
            outfile)

f = os.path.join(dir, 'test_host_seqs')
with open(f, 'w') as f:
    for d in (mammal, bird):
        for seq in d:
            f.write(seq + '\n')

count_evaluate('MAMMAL', mammal, mammal_host_freqs, bird_host_freqs)
count_evaluate('MAMMAL CONTROL', mammal_control, mammal_host_freqs, 
               bird_host_freqs)
count_evaluate('BIRD', bird, 
               bird_host_freqs, mammal_host_freqs)
count_evaluate('BIRD CONTROL', bird_control, 
               bird_host_freqs, mammal_host_freqs)

