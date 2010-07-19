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
            elmseq, freq = line.strip().split('\t')
            seq = elmseq.split(':')[1]
            if counts[seq] > 500:
                freqs[seq] = float(freq)
    return freqs

def evaluate(name, flu_seqs, host1_freqs, host2_freqs):
    """Test hypothesis"""

    count = 0
    total = 0
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
    print name, count, total, float(count)/float(total)
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
            for elmseq in protein_elm[protein]:
                elm, seq = elmseq.split(':')
                f.write(protein + '\t' + elmseq + '\t'
                        + elm + '\t' + seq + '\n')

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

    other_seqs = {}
    for protein in cmp_cons:
        for elm in cmp_cons[protein]:
            seq = elm.split(':')[-1]
            other_seqs[seq] = True
        
    uniq = defaultdict(dict)
    for protein in all_cons:
        for elm in all_cons[protein]:
            seq = elm.split(':')[-1]
            if seq not in other_seqs:
                uniq[seq] = True
    return uniq

def count_it(uniq):
    return len(uniq)

def mk_control(species_all_cons, species_uniq):
    """Make a control using ELMs that are on 90% of this flu,
       but not unique and not conserved on all bird flus."""

    control = defaultdict(dict)
    for protein in species_all_cons:
        for elm in species_all_cons[protein]:
            seq = elm.split(':')[-1]
            if seq not in species_uniq:
                control[seq] = True
    return control    

dir = 'working/Jul19/'
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

print 'MAMMAL', len(mammal_uniq)
print 'BIRD', len(bird_uniq)
# print 'ALL', count_it(all_cons)
# print 'MAMMAL ALL', count_it(mammal_all_cons)
# print 'BIRD ALL', count_it(bird_all_cons)
# sys.exit(0)
# what is not considered uniq?
mammal_control = mk_control(mammal_all_cons, mammal_uniq)
bird_control = mk_control(bird_all_cons, bird_uniq)

print 'MAMMAL CONTROL', len(mammal_control)
print 'BIRD CONTROL', len(bird_control)

# dump_it(dir, mammal_uniq, 'mammal_uniq')
# dump_it(dir, bird_uniq, 'bird_uniq')
# dump_it(dir, mammal_all_cons, 'control')
# dump_it(dir, mammal_control, 'mammal_control')
# dump_it(dir, bird_control, 'bird_control')

protein = ''

mammal_pre = mammal_uniq
bird_pre = bird_uniq
mammal_control_pre = mammal_control
bird_control_pre = bird_control

mammal_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.init.elm_aa_freq'),
                              os.path.join(dir, 'elmdict_H_sapiens.init'))
bird_host_freqs = get_freqs(os.path.join(dir, 'Gallus_gallus.init.elm_aa_freq'),
                            os.path.join(dir, 'elmdict_Gallus_gallus.init'))

mammal = set(mammal_pre.keys()) - set(bird_pre.keys())
bird = set(bird_pre.keys()) - set(mammal_pre.keys())

mammal_control_pre2 = set(mammal_control_pre.keys()) - set(mammal_pre.keys())
bird_control_pre2 = set(bird_control_pre.keys()) - set(bird_pre.keys())
both_control = bird_control_pre2 & mammal_control_pre2
mammal_control = mammal_control_pre2# - both_control
bird_control = bird_control_pre2# - both_control

print 'intr', len(mammal_control & mammal)
print 'intr', len(bird_control & bird)
print 'intr', len(mammal & bird)
print 'intr', len(mammal_control & bird_control)

m1, m2 = evaluate('MAMMAL', mammal, mammal_host_freqs, bird_host_freqs)
c1, c2 = evaluate('MAMMAL CONTROL', mammal_control, mammal_host_freqs, bird_host_freqs)
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
