"""What ELM sequences are present on all mammal/bird flus,
   but missing from all of the opposite class?"""
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

def get_uniq(all_cons, cons, low_cons_cut):
    """What protein/ELM pairs are in all_cons, but are below 
       low_cons_cut for the other flu strains in cons?

       Ex. mammal is all_cons, and bird is cons"""

    uniq = defaultdict(dict)
    for protein in all_cons:
        for elm in all_cons[protein]:
            found = False
            for key in cons:
                if protein in cons[key]:
                    if elm in cons[key][protein]:
                        if cons[key][protein][elm] > low_cons_cut:
                            found = True
                            break
            if not found:
                uniq[protein][elm] = True
    return uniq

def count_it(uniq):
    count = 0
    for protein in uniq:
        for elm in uniq[protein]:
            count += 1
    return count

dir = 'working/Jul12'
years = range(2000,2011,1)
mammal_hosts = ('human',)
mammal_strains = ('H5N1',)
bird_hosts = ('chicken',)
bird_strains = ('H5N1',)

limit = 25
cons_cut = float(80)
low_cons_cut = float(75)
mammal_protein_counts = get_protein_counts(dir, mammal_hosts,
                                           mammal_strains, years)
bird_protein_counts = get_protein_counts(dir, bird_hosts,
                                           bird_strains, years)
mammal_cons = get_cons_elms(dir, mammal_protein_counts, limit)
bird_cons = get_cons_elms(dir, bird_protein_counts, limit)
mammal_all_cons = get_all_cons(mammal_cons, mammal_protein_counts, limit, cons_cut)
bird_all_cons = get_all_cons(bird_cons, bird_protein_counts, limit, cons_cut)
mammal_uniq = get_uniq(mammal_all_cons, bird_cons, low_cons_cut)
bird_uniq = get_uniq(bird_all_cons, mammal_cons, low_cons_cut)

print 'MAMMAL', count_it(mammal_uniq)
print 'BIRD', count_it(bird_uniq)


