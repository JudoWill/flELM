"""Split ELM/seqs into those on all strains,
   and those specific to 90% mammal but not
   in 90% bird."""
import os, utils
from collections import defaultdict

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
            for elm in protein_elm[protein]:
                f.write(protein + '\t' + elm + '\n')

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
       Ex. mammal is all_cons, and bird is cmp_cons"""

    uniq = defaultdict(dict)
    for protein in all_cons:
        if protein in cmp_cons:
            for elm in all_cons[protein]:
                if elm not in cmp_cons[protein]:
                    uniq[protein][elm] = True
    return uniq

def count_it(uniq):
    count = 0
    for protein in uniq:
        for elm in uniq[protein]:
            count += 1
    return count

dir = 'working/Jul12/'
years = range(2000,2011,1)
mammal_hosts = ('human','swine','horse')
mammal_strains = ('H5N1','H3N2','H3N8','H1N1')
bird_hosts = ('chicken','duck')
bird_strains = ('H5N1','H9N2')

limit = 50
cons_cut = float(90)
low_cons_cut = float(60)
mammal_protein_counts = get_protein_counts(dir, mammal_hosts,
                                           mammal_strains, years)
bird_protein_counts = get_protein_counts(dir, bird_hosts,
                                         bird_strains, years)
mammal_cons = get_cons_elms(dir, mammal_protein_counts, limit)
bird_cons = get_cons_elms(dir, bird_protein_counts, limit)
mammal_all_cons = get_all_cons(mammal_cons, mammal_protein_counts, limit, cons_cut)
bird_all_cons = get_all_cons(bird_cons, bird_protein_counts, limit, cons_cut)
mammal_uniq = get_uniq(mammal_all_cons, bird_all_cons)
bird_uniq = get_uniq(bird_all_cons, mammal_all_cons)

mammal_protein_counts.update(bird_protein_counts)
mammal_cons.update(bird_cons)
all_cons = get_all_cons(mammal_cons, mammal_protein_counts, limit, float(95))

print 'MAMMAL', count_it(mammal_uniq)
print 'BIRD', count_it(bird_uniq)
print 'ALL', count_it(all_cons)
print 'MAMMAL ALL', count_it(mammal_all_cons)
print 'BIRD ALL', count_it(bird_all_cons)

dump_it(dir, mammal_uniq, 'mammal_uniq')
dump_it(dir, bird_uniq, 'bird_uniq')
dump_it(dir, mammal_all_cons, 'control')
