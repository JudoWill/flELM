"""Try to recover host phylogeny by
   finding an ELM frequency cutoff.
"""
import sys, global_settings, os
from collections import defaultdict

elm_freq_cutoff = float(sys.argv[1])
results_dir = sys.argv[2]

def test_distances(distances):
    """Do these distances match the phylogeny?
       This function is not used yet."""

    distance_is_correct = True
    
    # check human/chimp

    # check chicken/finch

    return distance_is_correct

def do_elm_cutoff(results_dir, elm_freq_cutoff):
    """Return ELMs that are below the freq cutoff
       for all hosts."""

    suffix = '.init.elm_aa_freq'
    elm_passes = defaultdict(dict)
    for genome in global_settings.TEST_GENOMES:
        file = os.path.join(results_dir, genome + suffix)
        with open(file) as f:
            for line in f:
                elm, freq = line.strip().split('\t')
                if float(freq) < elm_freq_cutoff:
                    elm_passes[elm][genome] = True
    for elm in elm_passes:
        if len(elm_passes[elm]) == len(global_settings.TEST_GENOMES):
            print elm + '\tstuff'

do_elm_cutoff(results_dir, elm_freq_cutoff)

    
    
