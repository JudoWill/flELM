from collections import defaultdict
import utils_motif, sys, utils_graph
import utils_stats

def get_aa_freqs(afile):
    d = {}
    elms = {}
    with open(afile) as f:
        for line in f:
            [elm, v] = line.strip().split('\t')
            d[elm] = float(v)
            elms[elm] = True
    return (d, elms)

def check_gtr(elm, elm2fracs):
    return (elm2fracs[elm]['human'] > elm2fracs[elm]['chicken'] and elm2fracs[elm]['human'] > elm2fracs[elm]['finch'] and elm2fracs[elm]['swine'] > elm2fracs[elm]['chicken'] and elm2fracs[elm]['swine'] > elm2fracs[elm]['finch'])

def check_less(elm, elm2fracs):
    return (elm2fracs[elm]['human'] < elm2fracs[elm]['chicken'] and elm2fracs[elm]['human'] < elm2fracs[elm]['finch'] and elm2fracs[elm]['swine'] < elm2fracs[elm]['chicken'] and elm2fracs[elm]['swine'] < elm2fracs[elm]['finch'])

elm2fracs = {}
with open(sys.argv[1]) as f:
    f.readline()
    for line in f:
        [elm, human, swine, chicken, finch] = line.strip().split('\t')
        elm2fracs[elm] = {'human':float(human),
                          'swine':float(swine),
                          'chicken':float(chicken),
                          'finch':float(finch)}
# del elm2fracs['LIG_PDZ_3']
# #del elm2fracs['MOD_GSK3_1']
# #del elm2fracs['MOD_CK1_1']
# #del elm2fracs['MOD_CK2_1']

aa_freqs = {'human':get_aa_freqs('results/H_sapiens.elm_aa_freq'),
            'chicken':get_aa_freqs('results/Gallus_gallus.elm_aa_freq'),
            'finch':get_aa_freqs('results/Taeniopygia_guttata.elm_aa_freq'),
            'swine':get_aa_freqs('results/Sus_scrofa.elm_aa_freq')}
freq_elms = {}
for k in aa_freqs:
    d, e = aa_freqs[k]
    for elm in e:
        freq_elms[elm] = True
# del freq_elms['LIG_PDZ_3']
# del freq_elms['MOD_CK1_1']
# del freq_elms['MOD_CK2_1']
# del freq_elms['MOD_GSK3_1']
elm2freq = {}
for elm in freq_elms:
    elm2freq[elm] = {}
    for s in aa_freqs:
        if elm in aa_freqs[s][0]:
            elm2freq[elm][s] = aa_freqs[s][0][elm]
        else:
            elm2freq[elm][s] = float(0)

cut = sys.argv[2]
d = {'ELM':True}
swine_H1N1_elms = utils_motif.protein2annotation('results/swine.H1N1.elms.' + cut, d)
swine_H3N2_elms = utils_motif.protein2annotation('results/swine.H3N2.elms.' + cut, d)
swine = [swine_H1N1_elms, swine_H3N2_elms]

human_H1N1_elms = utils_motif.protein2annotation('results/human.H1N1.elms.' + cut, d)
human_H3N2_elms = utils_motif.protein2annotation('results/human.H3N2.elms.' + cut, d)
human_H5N1_elms = utils_motif.protein2annotation('results/human.H5N1.elms.' + cut, d)
human = [human_H1N1_elms, human_H3N2_elms, human_H5N1_elms]

chicken_H5N1_elms = utils_motif.protein2annotation('results/chicken.H5N1.elms.' + cut, d)
chicken_H9N2_elms = utils_motif.protein2annotation('results/chicken.H9N2.elms.' + cut, d)
chicken = [chicken_H5N1_elms, chicken_H9N2_elms]

duck_H5N1_elms = utils_motif.protein2annotation('results/duck.H5N1.elms.' + cut, d)
duck_H9N2_elms = utils_motif.protein2annotation('results/duck.H9N2.elms.' + cut, d)
duck = [duck_H5N1_elms, duck_H9N2_elms]

# these {}s are not completely right
# b/c they may miss proteins not present
# in the starting {}
common_all = defaultdict(dict)
common_all_elms = {}
for protein in human[0]:
    for elm in human[0][protein]:
        not_found = False
        for h in human[1:]:
            if protein in h:
                if not elm in h[protein]:
                    not_found = True
        for s in swine:
            if protein in s:
                if not elm in s[protein]:
                    not_found = True
        for c in chicken:
            if protein in c:
                if not elm in c[protein]:
                    not_found = True
        for d in duck:
            if protein in d:
                if not elm in d[protein]:
                    not_found = True
        if not not_found:
            common_all[protein][elm] = True
            common_all_elms[elm] = True
print common_all.keys()
def test_it(elm, d, out):
    same = 0
    if elm in d:
        if check_gtr(elm, d):
            out.write(elm + '\tGTG\n')
        elif check_less(elm, d):
            out.write(elm + '\tLESS\n')
        else:
            same += 1
            out.write(elm + '\tSAME\n')
    return same

test_elms = common_all_elms
virus_elms_same = 0
non_virus_elms_same = 0
with open('mammal_bird.same.' + str(cut) + '.test', 'w') as f:
    for elm in test_elms:
        virus_elms_same += test_it(elm, elm2freq, f)

non_virus_elms = 0
with open('mammal_bird.same.' + str(cut) + '.notest', 'w') as f:       
    for elm in elm2freq:
        if not elm in test_elms:
            non_virus_elms += 1
            non_virus_elms_same += test_it(elm, elm2freq, f)

virus_elm_count = len(test_elms.keys())
with open(str(cut) + '.same.results', 'w') as f:
    p = utils_stats.fisher_positive_pval([virus_elm_count-
                                          virus_elms_same,
                                          virus_elms_same],
                                         [non_virus_elms-
                                          non_virus_elms_same,
                                          non_virus_elms_same])
    f.write('pvalue\t' + str(p) + '\n')
    f.write('virus\t' + str(virus_elm_count-virus_elms_same)
            + '\t' + str(virus_elms_same) + '\n')
    f.write('nvirus\t' + str(non_virus_elms-non_virus_elms_same)
            + '\t' + str(non_virus_elms_same) + '\n')
