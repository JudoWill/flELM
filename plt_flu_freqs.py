"""Plot sequence frequencies for flu"""
import global_settings, utils, os, random, math
from collections import defaultdict

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
            seq, freq = line.strip().split('\t')
            #seq = elmseq.split(':')[1]
            if float(freq) > float(0.0001):
                freqs[seq] = float(freq)
            #if float(freq) == float(.0001):
            #    print 'count', counts[seq]#freqs[seq] = float(freq)
            # elif counts[seq] == 1000:
            #     print 'freq', freq
                
    return freqs

def mk_plot(dir, protein, human_freqs, chicken_freqs, limit):
    """For this flu protein, plot MSA % coverage if above limit"""

    interesting_human = {}
    interesting_chicken = {}
    control = {}

    tmp_r = 'rfile' + str(random.randint(0,100))
    tmp_input1 = '1rinput' + str(random.randint(0,100))
    tmp_input2 = '2rinput' + str(random.randint(0,100))
    diff_input = 'diff_input' + str(random.randint(0,100))
    with open(tmp_input1, 'w') as f1:
        with open(tmp_input2, 'w') as f2:
            with open(diff_input, 'w') as df:
                for seq in set(human_freqs.keys()) | set(chicken_freqs.keys()):
                    if human_freqs[seq] > limit or chicken_freqs[seq] > limit:
                        diff = human_freqs[seq] - chicken_freqs[seq]
                        df.write(str(diff) + '\n')
                        if abs(diff) > float(10):
                            f1.write('%.10f\t%.10f\n' %
                                     (human_freqs[seq], chicken_freqs[seq]))
                            if diff > float(0):
                                interesting_human[seq] = True
                            else:
                                interesting_chicken[seq] = True
                        elif abs(diff) < float(1):
                            control[seq] = True
                        f2.write('%.10f\t%.10f\n' %
                                 (human_freqs[seq], chicken_freqs[seq]))
    outfile = os.path.join(dir, '_'.join(protein.split()) + '.png')
    with open(tmp_r, 'w') as f:
        f.write("d1<-read.delim('"
                + tmp_input1 + "', header=FALSE, sep='\\t')\n")
        f.write("d2<-read.delim('"
                + tmp_input2 + "', header=FALSE, sep='\\t')\n")
        f.write("png('" + outfile + "')\n")
        f.write("plot(d2,col='black',xlab='Human',ylab='Chicken',main='"
                + protein + "')\n")
        f.write("points(d1,col='red')\ndev.off()\nq()\n")
    os.system('R < ' + tmp_r + ' --no-save')
    os.system('rm ' + tmp_r + ' ' + tmp_input1 + ' ' + tmp_input2)
    
    outfile = os.path.join(dir, 'diff_' + '_'.join(protein.split()) + '.png')
    with open(tmp_r, 'w') as f:
        f.write("d<-read.delim('"
                + diff_input + "', header=FALSE, sep='\\t')\n")
        f.write("png('" + outfile + "')\n")
        f.write("plot(density(d[,1]),main='"
                + protein + "')\ndev.off()\nq()\n")
    os.system('R < ' + tmp_r + ' --no-save')
    #os.system('rm ' + tmp_r + ' ' + diff_input)

    return (interesting_human, interesting_chicken)

def mk_host_plot(dir, protein, human_freqs, chicken_freqs, limit,
                 human_host_freqs, chicken_host_freqs, interesting):
    """Plot host differences against flu differences."""

    interesting_human, interesting_chicken = interesting

    seqs = set(human_freqs.keys()) | set(chicken_freqs.keys())
    diffs_human = []
    diffs_control = []
    for seq in seqs:
        if human_freqs[seq] > limit or chicken_freqs[seq] > limit:
            if seq in interesting_human:
                if seq in human_host_freqs and seq in chicken_host_freqs:
                    virus_diff = human_freqs[seq]-chicken_freqs[seq]
                    host_diff =  float(-1) * (math.log(human_host_freqs[seq]) - 
                                              math.log(chicken_host_freqs[seq]))
                    diffs_human.append((virus_diff, host_diff))
            if seq in human_host_freqs and seq in chicken_host_freqs:
                virus_diff = human_freqs[seq]-chicken_freqs[seq]
                host_diff =  float(-1) * (math.log(human_host_freqs[seq]) - 
                                          math.log(chicken_host_freqs[seq]))
                if seq in interesting_human:
                    diffs_control.append((virus_diff, host_diff))
                elif virus_diff > float(0) and virus_diff < float(10):
                    diffs_control.append((virus_diff, host_diff))

    tmp_r = 'rfile' + str(random.randint(0,100))
    tmp_input1 = '1rinput' + str(random.randint(0,100))
    tmp_input2 = '2rinput' + str(random.randint(0,100))
    with open(tmp_input1, 'w') as f:
        for virus_diff, host_diff in diffs_human:
            f.write('%.10f\t%.10f\n' %
                    (virus_diff,
                     host_diff))
    with open(tmp_input2, 'w') as f:
        for virus_diff, host_diff in diffs_control:
            f.write('%.10f\t%.10f\n' %
                    (virus_diff,
                     host_diff))
    outfile = os.path.join(dir, 'HOST_' + '_'.join(protein.split()) + '.png')
    with open(tmp_r, 'w') as f:
        f.write("d1<-read.delim('"
                + tmp_input1 + "', header=FALSE, sep='\\t')\n")
        f.write("control<-read.delim('"
                + tmp_input2 + "', header=FALSE, sep='\\t')\n")
        f.write("png('" + outfile + "')\n")
        f.write("plot(control,xlab='Virus diff',ylab='Host diff', main='"
                + protein + "')\n")
        f.write("points(d1,col='red')\ndev.off()\nq()\n")
    os.system('R < ' + tmp_r + ' --no-save')
    os.system('rm ' + tmp_r + ' ' + tmp_input1 + ' ' + tmp_input2)

def get_avg_cons(use_strains, use_files):
    """Make {} of protein to avg SEQ coverage"""

    cons = {}
    for file in use_files:
        with open(file) as f:
            for line in f:
                protein, elm, cons_st = line.strip().split('\t')
                if file in use_strains[protein]:
                    if file not in cons:
                        cons[file] = defaultdict(dict)
                    cons[file][protein][elm] = float(cons_st)
    freqs = {}
    for protein in use_strains:
        freqs[protein] = defaultdict(utils.init_zero)
        for file in use_strains[protein]:
            for elm in cons[file][protein]:
                freqs[protein][elm] += cons[file][protein][elm]
        for elm in freqs[protein]:
            freqs[protein][elm] /= float(len(use_strains[protein]))
    return freqs

def get_use_strains(dir, hosts, strains, years):
    """Find proteins on strains with enough sequences."""

    use_strains = defaultdict(dict)
    use_files = {}
    for host in hosts:
        for strain in strains:
            for year in years:
                key = '.'.join((host, strain, str(year)))
                f = os.path.join(dir, key + '.fa')
                if os.path.exists(f):
                    counts = defaultdict(dict)
                    for ID, seq in utils.fasta_iter(f):
                        protein = ID.split('.')[-1]
                        counts[protein][ID] = True
                    for protein in counts:
                        if len(counts[protein]) > global_settings.SEQ_LIMIT:
                            cons_file = f.replace('.fa',
                                                  '.elms.conservation')
                            use_strains[protein][cons_file] = True
                            use_files[cons_file] = True
    return (use_strains, use_files)

dir = 'working/Jul22/'
years = range(2000,2011,1)
human_hosts = ('human',)
human_strains = ('H5N1','H1N1','H3N2')
chicken_hosts = ('chicken',)
chicken_strains = ('H5N1','H9N2''H7N2')

human_protein2strain, human_files = get_use_strains(dir, human_hosts,
                                                    human_strains, years)
chicken_protein2strain, chicken_files = get_use_strains(dir, chicken_hosts,
                                                        chicken_strains, years)
human_cons = get_avg_cons(human_protein2strain, human_files)
chicken_cons = get_avg_cons(chicken_protein2strain, chicken_files)
limit = float(0)
interesting = {}
for protein in human_cons:
    interesting[protein] = mk_plot(dir, protein, human_cons[protein],
                                   chicken_cons[protein], limit)

human_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.simple.elm_aa_freq'),
                             os.path.join(dir, 'elmdict_H_sapiens.simple'))
chicken_host_freqs = get_freqs(os.path.join(dir, 
                                            'Gallus_gallus.simple.elm_aa_freq'),
                               os.path.join(dir, 'elmdict_Gallus_gallus.simple'))
# for protein in human_cons:
#      mk_host_plot(dir, 'human ' + protein, human_cons[protein], chicken_cons[protein], 
#                   limit,
#                   human_host_freqs, chicken_host_freqs,
#                   interesting[protein])
#      mk_host_plot(dir, 'chicken ' + protein, chicken_cons[protein], human_cons[protein], 
#                   limit,
#                   chicken_host_freqs, human_host_freqs,
#                   (interesting[protein][1], interesting[protein][0]))





