from collections import defaultdict
import utils_motif, sys, utils_graph
import utils_stats

ignore_elms = {}#'LIG_PDZ_3':True,
               #'MOD_CK1_1':True,
               #'MOD_CK2_1':True,
               #'MOD_GSK3_1':True}

def get_unique_elms(common_mammal, common_bird):
    use_mammal = {}
    for protein in common_mammal:
        if protein in common_bird:
            for elm in common_mammal[protein]:
                if not elm in common_bird[protein]:
                    use_mammal[elm] = True
    return use_mammal


def get_common_elms(d):
    """enter [] of species 2 [] of strains2protein2elms"""

    protein2elm = defaultdict(dict)
    elms = {}
    start_strain = d[0][0] # needs to have all proteins
    for protein in start_strain:
        for elm in start_strain[protein]:
            not_found = False
            for a_strain in d[0][1:]:
                if protein in a_strain:
                    if not elm in a_strain[protein]:
                        not_found = True
            for species in d[1:]:
                for a_strain in species:
                    if protein in a_strain:
                        if not elm in a_strain[protein]:
                            not_found = True
            if not not_found:
                protein2elm[protein][elm] = True
                elms[elm] = True
    return (protein2elm, elms)

def test_it(elm, d, out):
    same = 0
    count = 0
    if elm in d:
        count = 1
        if check_gtr(elm, d):
            out.write(elm + '\tGTR\n')
        elif check_less(elm, d):
            out.write(elm + '\tLESS\n')
        else:
            same = 1
            out.write(elm + '\tSAME\n')
    return (count, same)

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
    return (elm2fracs[elm]['human'] > elm2fracs[elm]['chicken'] and elm2fracs[elm]['human'] > elm2fracs[elm]['finch'] and elm2fracs[elm]['swine'] > elm2fracs[elm]['chicken'] and elm2fracs[elm]['swine'] > elm2fracs[elm]['finch'] and elm2fracs[elm]['horse'] > elm2fracs[elm]['chicken'] and elm2fracs[elm]['horse'] > elm2fracs[elm]['finch'])

def check_gtr_2fold(elm, elm2fracs):
    return (elm2fracs[elm]['human']/elm2fracs[elm]['chicken'] > float(1.05) and elm2fracs[elm]['human']/elm2fracs[elm]['finch'] > float(1.05) and elm2fracs[elm]['swine']/elm2fracs[elm]['chicken'] > float(1.05) and elm2fracs[elm]['swine']/elm2fracs[elm]['finch'] > float(1.05))

def check_less(elm, elm2fracs):
    return (elm2fracs[elm]['human'] < elm2fracs[elm]['chicken'] and elm2fracs[elm]['human'] < elm2fracs[elm]['finch'] and elm2fracs[elm]['swine'] < elm2fracs[elm]['chicken'] and elm2fracs[elm]['swine'] < elm2fracs[elm]['finch'] and elm2fracs[elm]['horse'] < elm2fracs[elm]['chicken'] and elm2fracs[elm]['horse'] < elm2fracs[elm]['finch'])

def check_less_2fold(elm, elm2fracs):
    return (elm2fracs[elm]['human']/elm2fracs[elm]['chicken'] > float(1.05) and elm2fracs[elm]['human']/elm2fracs[elm]['finch'] > float(1.05) and elm2fracs[elm]['swine']/elm2fracs[elm]['chicken'] > float(1.05) and elm2fracs[elm]['swine']/elm2fracs[elm]['finch'] > float(1.05))

elm2fracs = {}
with open(sys.argv[1]) as f:
    f.readline()
    for line in f:
        [elm, human, swine, chicken, finch] = line.strip().split('\t')
        elm2fracs[elm] = {'human':float(human),
                          'swine':float(swine),
                          'chicken':float(chicken),
                          'finch':float(finch)}

aa_freqs = {'human':get_aa_freqs('results/H_sapiens.elm_aa_freq'),
            'chicken':get_aa_freqs('results/Gallus_gallus.elm_aa_freq'),
            'finch':get_aa_freqs('results/Taeniopygia_guttata.elm_aa_freq'),
            'swine':get_aa_freqs('results/Sus_scrofa.elm_aa_freq'),
            'horse':get_aa_freqs('results/Equus_caballus.elm_aa_freq')}

freq_elms = {}
for k in aa_freqs:
    d, e = aa_freqs[k]
    for elm in e:
        freq_elms[elm] = True

elm2freq = {}
for elm in freq_elms:
    elm2freq[elm] = {}
    for s in aa_freqs:
        if elm in aa_freqs[s][0]:
            elm2freq[elm][s] = aa_freqs[s][0][elm]
        else:
            elm2freq[elm][s] = float(0.0000000000000000000001)

cut = sys.argv[2]

d = {'ELM':True}
swine_H1N1_elms = utils_motif.protein2annotation('results/swine.H1N1.elms.' + cut, d)
swine_H3N2_elms = utils_motif.protein2annotation('results/swine.H3N2.elms.' + cut, d)
swine = [swine_H1N1_elms, swine_H3N2_elms]

human_H1N1_elms = utils_motif.protein2annotation('results/human.H1N1.elms.' + cut, d)
human_H3N2_elms = utils_motif.protein2annotation('results/human.H3N2.elms.' + cut, d)
human_H5N1_elms = utils_motif.protein2annotation('results/human.H5N1.elms.' + cut, d)
human = [human_H1N1_elms, human_H3N2_elms, human_H5N1_elms]

horse_H3N8_elms = utils_motif.protein2annotation('results/equine.H3N8.elms.' + cut, d)
horse = [horse_H3N8_elms]

chicken_H5N1_elms = utils_motif.protein2annotation('results/chicken.H5N1.elms.' + cut, d)
chicken_H9N2_elms = utils_motif.protein2annotation('results/chicken.H9N2.elms.' + cut, d)
chicken = [chicken_H5N1_elms, chicken_H9N2_elms]

duck_H5N1_elms = utils_motif.protein2annotation('results/duck.H5N1.elms.' + cut, d)
duck_H9N2_elms = utils_motif.protein2annotation('results/duck.H9N2.elms.' + cut, d)
duck = [duck_H5N1_elms, duck_H9N2_elms]

swine_H1N1_elms_controled = utils_motif.protein2annotation('results/swine.H1N1.elms.' + cut + '.controled', d)
swine_H3N2_elms_controled = utils_motif.protein2annotation('results/swine.H3N2.elms.' + cut + '.controled', d)
swine_controled = [swine_H1N1_elms_controled, swine_H3N2_elms_controled]

human_H1N1_elms_controled = utils_motif.protein2annotation('results/human.H1N1.elms.' + cut + '.controled', d)
human_H3N2_elms_controled = utils_motif.protein2annotation('results/human.H3N2.elms.' + cut + '.controled', d)
human_H5N1_elms_controled = utils_motif.protein2annotation('results/human.H5N1.elms.' + cut + '.controled', d)
human_controled = [human_H1N1_elms_controled, human_H3N2_elms_controled, human_H5N1_elms_controled]

horse_H3N8_elms_controled = utils_motif.protein2annotation('results/equine.H3N8.elms.' + cut + '.controled', d)
horse_controled = [horse_H3N8_elms_controled]

chicken_H5N1_elms_controled = utils_motif.protein2annotation('results/chicken.H5N1.elms.' + cut + '.controled', d)
chicken_H9N2_elms_controled = utils_motif.protein2annotation('results/chicken.H9N2.elms.' + cut + '.controled', d)
chicken_controled = [chicken_H5N1_elms_controled, chicken_H9N2_elms_controled]

duck_H5N1_elms_controled = utils_motif.protein2annotation('results/duck.H5N1.elms.' + cut + '.controled', d)
duck_H9N2_elms_controled = utils_motif.protein2annotation('results/duck.H9N2.elms.' + cut + '.controled', d)
duck_controled = [duck_H5N1_elms_controled, duck_H9N2_elms_controled]

# these {}s are not completely right
# b/c they may miss proteins not present
# in the starting {}
ls = [human, swine, horse, chicken, duck]
common_all, common_all_elms = get_common_elms(ls)

ls = [human, swine, horse]
common_mammal, mammal_elms = get_common_elms(ls)

ls = [chicken, duck]
common_bird, bird_elms = get_common_elms(ls)

use_mammal = get_unique_elms(common_mammal, common_bird)
use_bird = get_unique_elms(common_bird, common_mammal)
utils_graph.dumpNodes('mammal' + str(cut), use_mammal)
utils_graph.dumpNodes('bird' + str(cut), use_bird)
utils_graph.dumpNodes('common' + str(cut), common_all_elms)

ls = [human_controled, swine_controled, 
      horse_controled, chicken_controled, duck_controled]
common_all_controled, common_all_elms_controled = get_common_elms(ls)

ls = [human_controled, swine_controled, horse_controled]
common_mammal_controled, mammal_elms_controled = get_common_elms(ls)

ls = [chicken_controled, duck_controled]
common_bird_controled, bird_elms_controled = get_common_elms(ls)

use_mammal_controled = get_unique_elms(common_mammal_controled, 
                                       common_bird_controled)
use_bird_controled = get_unique_elms(common_bird_controled, 
                                     common_mammal_controled)
utils_graph.dumpNodes('mammal_controled' + str(cut), use_mammal_controled)
utils_graph.dumpNodes('bird_controled' + str(cut), use_bird_controled)
utils_graph.dumpNodes('common_controled' + str(cut), common_all_elms_controled)

test_elms = utils_graph.unionLists([use_mammal, use_bird]) 
virus_elms_same = 0
virus_elm_count = 0
non_virus_elms_same = 0
with open('mammal_bird.different.' + str(cut) + '.test', 'w') as f:
    for elm in test_elms:
        if elm not in ignore_elms:
            count,same = test_it(elm, elm2freq, f)
            virus_elms_same += same
            virus_elm_count += count

non_virus_elms = 0
non_virus_elms_all = 0
control_elms = {}
with open('mammal_bird.different.' + str(cut) + '.notest', 'w') as f:       
    for elm in elm2freq:
        if not elm in test_elms:
            control_elms[elm] = True
            non_virus_elms_all += 1
            if not elm in ignore_elms:
                count,same = test_it(elm, elm2freq, f)
                non_virus_elms += count
                non_virus_elms_same += same            
                
diff_diff = virus_elm_count-virus_elms_same
diff_same = virus_elms_same + len(test_elms.keys())-virus_elm_count-len(utils_graph.intersectLists([test_elms,ignore_elms]))
diff_bg_diff = non_virus_elms-non_virus_elms_same
diff_bg_same = non_virus_elms_same + non_virus_elms_all-non_virus_elms-len(utils_graph.intersectLists([control_elms,ignore_elms]))
with open(str(cut) + '.different.results', 'w') as f:
    p = utils_stats.fisher_positive_pval([diff_diff,diff_same],
                                         [diff_bg_diff,diff_bg_same])
    f.write('pvalue\t' + str(p) + '\n')
    f.write('virus\t' + str(diff_diff)
            + '\t' + str(diff_same) + '\n')
    f.write('nvirus\t' + str(diff_bg_diff)
            + '\t' + str(diff_bg_same) + '\n')

test_elms_2 = {}
for c in common_all_elms:
    if not c in test_elms:
        test_elms_2[c] = True

virus_elms_same = 0
virus_elm_count = 0
non_virus_elms_same = 0
with open('mammal_bird.same.' + str(cut) + '.test', 'w') as f:
    for elm in test_elms_2:
        if not elm in ignore_elms:
            count,same = test_it(elm, elm2freq, f)
            virus_elms_same += same
            virus_elm_count += count

non_virus_elms = 0
non_virus_elms_all = 0
control_elms = {}
with open('mammal_bird.same.' + str(cut) + '.notest', 'w') as f:       
    for elm in elm2freq:
        if not elm in test_elms_2:
            control_elms[elm] = True
            non_virus_elms_all += 1
            if not elm in ignore_elms:
                count,same = test_it(elm, elm2freq, f)
                non_virus_elms += count
                non_virus_elms_same += same

same_diff = virus_elm_count-virus_elms_same
same_same = virus_elms_same + len(test_elms_2.keys())-virus_elm_count-len(utils_graph.intersectLists([test_elms_2,ignore_elms]))
same_bg_diff = non_virus_elms-non_virus_elms_same
same_bg_same = non_virus_elms_same + non_virus_elms_all-non_virus_elms-len(utils_graph.intersectLists([control_elms,ignore_elms]))
with open(str(cut) + '.same.results', 'w') as f:
    p = utils_stats.fisher_positive_pval([same_diff,same_same],
                                         [same_bg_diff,same_bg_same])
    f.write('pvalue\t' + str(p) + '\n')
    f.write('virus\t' + str(same_diff)
            + '\t' + str(same_same) + '\n')
    f.write('nvirus\t' + str(same_bg_diff)
            + '\t' + str(same_bg_same) + '\n')

with open(str(cut) + '.results', 'w') as f:
    p = utils_stats.fisher_positive_pval([diff_diff,diff_same],
                                         [same_diff,same_same])
    f.write('pvalue\t' + str(p) + '\n')
    f.write('virus\t' + str(diff_diff)
            + '\t' + str(diff_same) + '\n')
    f.write('nvirus\t' + str(same_diff)
            + '\t' + str(same_same) + '\n')
    
