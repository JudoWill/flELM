""" Plot comparitive histograms
    of sequence distributions.

    Focus on virus ELMs (flu in filename)
    and renomalize host ELMs
    based on virus ELMs

    This will output out seq count
    bar graphs for each ELM.
"""
import utils_plot, sys, os, utils, utils_distance
from collections import defaultdict

def main(args):
    file_species_pairs = []
    i = 1
    while i < len(args)-2:
        file_species_pairs.append([args[i], args[i+1]])
        i += 2

    cutoff = float(sys.argv[-2])
    plot_dir = sys.argv[-1]

    species2elms = {}
    virus2elms = {}
    # first grab virus ELMs
    for file, species in file_species_pairs:
        if file.find('flu') != -1:
            virus2elms[species] = utils.get_seq2count_dict(file, cutoff)
        else:
            species2elms[species] = True
    elms = {}
    for species in virus2elms:
        for elm in virus2elms[species]:
            elms[elm] = True
    for species in species2elms:
        species2elms[species] = utils.get_seq2count_dict_for_seqs(file, 
                                                                  cutoff,
                                                                  virus2elms)
    for virus in virus2elms:
        species2elms[virus] = virus2elms[virus]
    for elm in elms:
        if utils.check_ones(species2elms, elm):
            if utils_distance.distance_elms(species2elms['Sus_scrofa'][elm],
                                            species2elms['H_sapiens'][elm]) > float(-1) or utils_distance.distance_elms(species2elms['Sus_scrofa'][elm],
                                                                                                                        species2elms['Gallus_gallus'][elm]) > float(0):
                utils_plot.elm_host_barplot(species2elms, elm,
                                            os.path.join(plot_dir,
                                                         elm + '.virus_hosts.png'))

if __name__ == '__main__': main(sys.argv)

