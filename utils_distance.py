""" Functions for weights and distances
    between ELM freq vectors
"""
import itertools, math
from collections import defaultdict

def get_elements(dict1, dict2):
    """ this should be replaced with a list union """

    elements = {}
    for a_dict in [dict1, dict2]:
        for element in a_dict:
            elements[element] = True
    return elements

def distance_elms(seq2freq_dict_1, seq2freq_dict_2):
    """ compute the distance btwn 2 species
        for an ELM as
        sqrt( sum (seq_freq1-seq_freq2)**2 ) """

    distance = float(0)
    seqs = get_elements(seq2freq_dict_1,
                        seq2freq_dict_2)
    
    for seq in seqs:
        if seq in seq2freq_dict_1 and seq in seq2freq_dict_2:
            distance += (seq2freq_dict_1[seq]
                         - seq2freq_dict_1[seq])**2
        elif seq in seq2freq_dict_1:
            distance += seq2freq_dict_1[seq]**2
        else:
            distance += seq2freq_dict_2[seq]**2
    return math.sqrt(distance)

def distance_elms_weighted(seq2freq_dict_1, seq2freq_dict_2):
    """ compute the distance btwn 2 species
        for an ELM as
        sqrt( sum (seq_freq1-seq_freq2)**2 ) 
        weight the result by the sum of 
        the frequencies 
    """

    distance = float(0)
    seqs = get_elements(seq2freq_dict_1,
                        seq2freq_dict_2)
    
    for seq in seqs:
        if seq in seq2freq_dict_1 and seq in seq2freq_dict_2:
            distance += (seq2freq_dict_1[seq]
                         + seq2freq_dict_1[seq])*(seq2freq_dict_1[seq]
                                                  - seq2freq_dict_1[seq])**2
        elif seq in seq2freq_dict_1:
            distance += seq2freq_dict_1[seq]**3
        else:
            distance += seq2freq_dict_2[seq]**3
    return math.sqrt(distance)

def distance_species(species2elms_1, species2elms_2):
    """ compute the distance between species as
        a sum of the ELM distances """

    distance = float(0)
    elms = get_elements(species2elms_1,
                        species2elms_2)
    for elm in elms:
        if elm in species2elms_1 and elm in species2elms_2:
            distance += distance_elms(species2elms_1[elm],
                                      species2elms_2[elm])
        elif elm in species2elms_1:
            for seq in species2elms_1[elm]:
                distance += species2elms_1[elm][seq]
        else:
            for seq in species2elms_2[elm]:
                distance += species2elms_2[elm][seq]
    return distance

def distance_species_weighted(species2elms_1, species2elms_2):
    """ compute the distance between species as
        a sum of the ELM distances.
        use weighted_elms """

    distance = float(0)
    elms = get_elements(species2elms_1,
                        species2elms_2)
    for elm in elms:
        if elm in species2elms_1 and elm in species2elms_2:
            distance += distance_elms_weighted(species2elms_1[elm],
                                               species2elms_2[elm])
        elif elm in species2elms_1:
            for seq in species2elms_1[elm]:
                distance += species2elms_1[elm][seq]
        else:
            for seq in species2elms_2[elm]:
                distance += species2elms_2[elm][seq]
    return distance

def elm_distance_matrix(species2elms):
    """ return {} of elm to distances between genome pairs
        ex. d[elm][dog][human] = .4
    """
    
    mat = defaultdict(dict)
    for species1, species2 in itertools.combinations(species2elms.keys(), 2):
        for elm in get_elements(species2elms[species1],
                                species2elms[species2]):
            if elm in species2elms[species1] and elm in species2elms[species2]:
                distance = distance_elms(species2elms[species1][elm],
                                          species2elms[species2][elm])
            elif elm in species2elms[species1]:
                for seq in species2elms[species1][elm]:
                    distance = species2elms[species1][elm][seq]
            else:
                for seq in species2elms[species2][elm]:
                    distance = species2elms[species2][elm][seq]
            mat[elm][species1+':'+species2] = distance
    return mat

def distance_matrix(species2elms):
    """ return {} of species distances """

    mat = defaultdict(dict)
    for species1, species2 in itertools.combinations(species2elms.keys(), 2):
        distance = distance_species_weighted(species2elms[species1],
                                             species2elms[species2])
        mat[species1][species2] = distance
        mat[species2][species1] = distance
    return mat
            
