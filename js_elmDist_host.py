"""Count ELMs (not ELM sequences) on hosts
   and make phylogeny use Jensen-Shannon
   divergence."""
import itertools, utils, global_settings, sys, utils_plot, os
from collections import defaultdict

results_dir = sys.argv[1] # 'working/runs/Jun24/'
out_file = sys.argv[2] # 'working/runs/Jun24/js_elm_host_phylogeny_new.png'
use_freqs = sys.argv[3] # T|F file is init.elm_aa_freq

def get_host_freqs(ls_of_hosts):
    """Grab ELM frequencies"""
    print 'here'
    host_elmFreqs = defaultdict(dict)
    seen_elms = defaultdict(dict)
    for host in ls_of_hosts:
        fname = os.path.join(results_dir,
                             host + '.init.elm_aa_freq')
        with open(fname) as f:
            for line in f:
                (elm, freq) = line.strip().split('\t')
                host_elmFreqs[host][elm] = float(freq)
                seen_elms[elm][host] = True
    use_elms = {}
    for elm in seen_elms:
#        if len(seen_elms[elm]) == len(ls_of_hosts):
        use_elms[elm] = True
    return (host_elmFreqs, use_elms)

def get_host_counts(ls_of_hosts):
    """Count total ELM hits for each
       host in ls_of_hosts.
       Return counts & ELMs found."""

    host_elmCounts = {}
    seen_elms = defaultdict(dict)
    for host in ls_of_hosts:
        host_elmCounts[host] = defaultdict(utils.init_zero)
        with open(results_dir + 'elmdict_' + host + '.init') as f:
            for line in f:
                (elm, seq, count, fq) = line.strip().split('\t')
                host_elmCounts[host][elm] += int(count)
                seen_elms[elm][host] = True
    use_elms = {}
    for elm in seen_elms:
#        if len(seen_elms[elm]) == len(ls_of_hosts):
        use_elms[elm] = True
    print len(use_elms)
    return (host_elmCounts, use_elms)

hosts = global_settings.TEST_GENOMES#('H_sapiens', 'Macaca_mulatta', 'Pan_troglodytes', 'R_norvegicus', 'M_musculus')
if use_freqs == 'T':
    host_elmCounts, elms = get_host_freqs(hosts)
else:
    host_elmCounts, elms = get_host_counts(hosts)
host_vecs = utils.mk_count_vecs(host_elmCounts, elms)
host_dists = utils.mk_count_dists(host_vecs)
utils_plot.phylogeny_js(out_file, host_dists)
