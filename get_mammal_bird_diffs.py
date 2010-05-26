from collections import defaultdict
import utils_motif, sys, utils_graph

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

common_mammal = defaultdict(dict)
mammal_elms = {}
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
        if not not_found:
            common_mammal[protein][elm] = True
            mammal_elms[elm] = True

common_bird = defaultdict(dict)
bird_elms = {}
for protein in chicken[0]:
    for elm in chicken[0][protein]:
        not_found = False
        for c in chicken[1:]:
            if protein in c:
                if not elm in c[protein]:
                    not_found = True
        for d in duck:
            if protein in d:
                if not elm in d[protein]:
                    not_found = True
        if not not_found:
            common_bird[protein][elm] = True
            bird_elms[elm] = True

use_elms = {}
control_elms = {}
for protein in common_mammal:
    if protein in common_bird:
        for elm in common_mammal[protein]:
            if not elm in common_bird[protein]:
                use_elms[elm] = True#print protein + '\t' + elm
            else:
                control_elms[elm] = True
for protein in common_bird:
    if protein in common_mammal:
        for elm in common_bird[protein]:
            if not elm in common_mammal[protein]:
                use_elms[elm] = True
            else:
                control_elms[elm] = True

# test_elms = {}
# with open('mammal_bird.test', 'w') as f:
#     for elm in use_elms:
#         if elm in elm2fracs:
#             if elm in mammal_elms and elm in bird_elms:
#                      control_elms[elm] = True
#             else:
#                 test_elms[elm] = True
#                 if check_gtr(elm, elm2fracs):
#                     f.write(elm + '\tGTR\n')
#                 elif check_less(elm, elm2fracs):
#                     f.write(elm + '\tLESS\n')
#                 else:
#                     f.write(elm + '\tSAME\n')
#         else:
#             if (elm in mammal_elms and not elm in bird_elms) or (elm in bird_elms and not elm in mammal_elms):
#                 print elm

# with open('mammal_bird.notTest', 'w') as f:
#     for elm in control_elms:
#         if not elm in test_elms:
#             if elm in elm2fracs:
#                 if check_gtr(elm, elm2fracs):
#                     f.write(elm + '\tGTR\n')
#                 elif check_less(elm, elm2fracs):
#                     f.write(elm + '\tLESS\n')
#                 else:
#                     f.write(elm + '\tSAME\n')

test_elms = {}
with open('mammal_bird.' + cut + '.test', 'w') as f:
    for elm in utils_graph.intersectLists([use_elms,freq_elms]):
         if elm in mammal_elms and elm in bird_elms:
                     control_elms[elm] = True
         elif elm in elm2fracs:
             test_elms[elm] = True
             if check_gtr(elm, elm2fracs):
                 f.write(elm + '\tGTR\n')
             elif check_less(elm, elm2fracs):
                 f.write(elm + '\tLESS\n')
             else:
                 f.write(elm + '\tSAME\n')
         else:
             test_elms[elm] = True
             if check_gtr(elm, elm2freq):
                 f.write(elm + '\tGTR\n')
             elif check_less(elm, elm2freq):
                 f.write(elm + '\tLESS\n')
             else:
                 f.write(elm + '\tSAME\n')

with open('mammal_bird.' + cut + '.notTest', 'w') as f:
     for elm in freq_elms:
         if not elm in test_elms:
             if elm in elm2fracs:
                 if check_gtr(elm, elm2fracs):
                     f.write(elm + '\tGTR\n')
                 elif check_less(elm, elm2fracs):
                     f.write(elm + '\tLESS\n')
                 else:
                     f.write(elm + '\tSAME\n')
             else:
                 if check_gtr(elm, elm2freq):
                  f.write(elm + '\tGTR\n')
                 elif check_less(elm, elm2freq):
                     f.write(elm + '\tLESS\n')
                 else:
                     f.write(elm + '\tSAME\n')

# with open('mammal_bird.notTest', 'w') as f:
#     for elm in utils_graph.intersectLists([control_elms,freq_elms]):
#         if not elm in test_elms:
#             if elm in elm2fracs:
#                 if check_gtr(elm, elm2fracs):
#                     f.write(elm + '\tGTR\n')
#                 elif check_less(elm, elm2fracs):
#                     f.write(elm + '\tLESS\n')
#                 else:
#                     f.write(elm + '\tSAME\n')
#             else:
#                 if check_gtr(elm, elm2freq):
#                  f.write(elm + '\tGTR\n')
#                 elif check_less(elm, elm2freq):
#                     f.write(elm + '\tLESS\n')
#                 else:
#                     f.write(elm + '\tSAME\n')

utils_graph.dumpNodes('common_bird', utils_graph.intersectLists([elm2fracs,bird_elms]))
utils_graph.dumpNodes('common_mammal', utils_graph.intersectLists([elm2fracs,mammal_elms]))
        

