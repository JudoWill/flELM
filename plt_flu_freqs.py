"""Plot sequence frequencies for flu"""
import global_settings, utils, os, random, math, utils_graph, itertools, utils_stats
from collections import defaultdict

def get_x(ls, ten_per):
    """Grab top ten_per percent from this list"""

    counter = 0
    fill = {}
    vals = []
    while counter < ten_per:
        fill[ ls[counter][1] ] = True
        val = ls[counter][0]
        vals.append(val)
        counter += 1
    if val == ls[counter][0]:
        while True:
            fill[ ls[counter][1] ] = True
            vals.append( ls[counter][0] )
            counter += 1
            if val != ls[counter][0]:
                break
    return (fill, vals)

def write_vals(vals, file):
    """Write values in this list"""

    with open(file, 'w') as f:
        for v in vals:
            f.write(str(v) + '\n')

def plot_extremes(foreground_vals, control_vals, outfile):
    """Plot extemes of flu diff distribution"""

    fore = 'tmpfore' + str(random.randint(0,100))
    control = 'tmpcontrol' + str(random.randint(0,100))
    write_vals(foreground_vals, fore)
    write_vals(control_vals, control)
    tmpr = 'tmpr' + str(random.randint(0,100))
    with open(tmpr, 'w') as f:
        f.write("fore<-read.delim('"
                + fore + "', header=FALSE, sep='\\t')\n")
        f.write("control<-read.delim('"
                + control + "', header=FALSE, sep='\\t')\n")
        f.write("png('" + outfile + "')\n")
        f.write("plot.multi.dens(list(fore[,1],control[,1]))\ndev.off()\nq()\n")
    os.system('R < ' + tmpr + ' --no-save')
    os.system('rm ' + tmpr + ' ' + fore + ' ' + control)

def get_fore_control(dir, protein, this_freqs, that_freqs, outfile):
    """For this flu protein, grab top/bottom 10% of MSA % diffs"""

    ls = []
    for seq in this_freqs:
        diff = this_freqs[seq] - that_freqs[seq]
        ls.append((diff,seq))
    ls.sort()
    ten_per = int(float(.05)*float(len(ls)))
    
    # get bottom
    control, control_vals = get_x(ls, ten_per)

    # get top
    ls.reverse()
    foreground, foreground_vals = get_x(ls, ten_per)
    
    plot_extremes(foreground_vals, control_vals, outfile)

    return (foreground, control)

def dump_input(seqs, f1, f2, file):
    """Write host diffs one per line for seq in seqs"""

    vals = []
    with open(file, 'w') as f:
        for seq in seqs:
            if seq in f1 and seq in f2:
                val = float(-1)*(math.log(f1[seq],100)-math.log(f2[seq],100))
                f.write('%.20f\n' %
                        (val))
                vals.append(val)
    return vals

def plot_fore_control(outfile, protein, fore, control, 
                      this_host_freqs, that_host_freqs):
    """Plot density for flu foreground vs control for host differences"""

    fore_input = 'fore' + str(random.randint(0,100))
    control_input = 'control' + str(random.randint(0,100))
    fore_vals = dump_input(fore, this_host_freqs, that_host_freqs, fore_input)
    control_vals = dump_input(control, this_host_freqs, that_host_freqs, control_input)

    rfile = 'tmpr' + str(random.randint(0,100))
    with open(rfile, 'w') as f:
        f.write("fore<-read.delim('"
                + fore_input + "', header=FALSE, sep='\\t')\n")
        f.write("control<-read.delim('"
                + control_input + "', header=FALSE, sep='\\t')\n")
        f.write("png('" + outfile + "')\n")
        f.write("plot.multi.dens(list(fore[,1],control[,1]))\ndev.off()\nq()\n")
    os.system('R < ' + rfile + ' --no-save')
    os.system('rm ' + rfile + ' ' + fore_input + ' ' + control_input)
    return (fore_vals, control_vals)

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
            #if float(freq) > float(0.0001):
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
human_host_freqs = get_freqs(os.path.join(dir, 'H_sapiens.simple.elm_aa_freq'),
                             os.path.join(dir, 'elmdict_H_sapiens.simple'))
chicken_host_freqs = get_freqs(os.path.join(dir, 
                                            'Gallus_gallus.simple.elm_aa_freq'),
                               os.path.join(dir, 'elmdict_Gallus_gallus.simple'))

masters_human = []
masters_chicken = []
for protein in human_cons:
    # interesting[protein] = mk_plot(dir, protein, human_cons[protein],
    #                                chicken_cons[protein], limit)
    human_fore, human_control = get_fore_control(dir, protein, 
                                                 human_cons[protein],
                                                 chicken_cons[protein],
                                                 os.path.join('human_' + protein + '.extreme.png'))
    chicken_fore, chicken_control = get_fore_control(dir, protein, 
                                                     chicken_cons[protein],
                                                     human_cons[protein],
                                                     os.path.join(dir,
                                                                  'chicken_' + protein + '.extreme.png'))
    outfile = os.path.join(dir, 'human_' + protein + '.png')
    human_fore_vals, human_control_vals = plot_fore_control(outfile, protein, 
                                                            human_fore, human_control, 
                                                            human_host_freqs, 
                                                            chicken_host_freqs)
    outfile = os.path.join(dir, 'chicken_' + protein + '.png')
    chicken_fore_vals, chicken_control_vals  = plot_fore_control(outfile, protein, 
                                                                 chicken_fore, 
                                                                 chicken_control, 
                                                                 chicken_host_freqs, 
                                                                 human_host_freqs)
    if protein in ('polymerase PB1', 'neuraminidase', 'polymerase PB2',
                   'nonstructural protein 2', 'polymerase PA'):
        human_pval = 'fore gtr ' + str(len(human_fore_vals)) + ' ' + str(len(human_control_vals) )+ ' ' + str(utils_stats.wilcox_gtr(human_fore_vals, human_control_vals))
        chicken_pval = 'fore gtr ' + str(len(chicken_fore_vals)) + ' ' + str(len(chicken_control_vals)) + ' ' + str(utils_stats.wilcox_gtr(chicken_fore_vals, chicken_control_vals))
    else:
        human_pval = 'ctrl gtr ' + str(len(human_control_vals)) + ' ' + str(len(human_fore_vals)) + ' ' + str(utils_stats.wilcox_gtr(human_control_vals, human_fore_vals))
        chicken_pval = 'ctrl gtr ' + str(len(chicken_control_vals)) + ' ' + str(len(chicken_fore_vals)) + ' ' + str(utils_stats.wilcox_gtr(chicken_control_vals, chicken_fore_vals))
    print 'human pval', protein, human_pval
    print 'chicken pval', protein, chicken_pval
    # print(protein, len(human_fore), len(human_control), 
    #       len(chicken_fore), len(chicken_control),
    #       len(set(human_fore.keys()) & set(chicken_fore.keys())))
    masters_human.append((human_fore, human_control))
    masters_chicken.append((chicken_fore, chicken_control))
# for p1,p2 in itertools.combinations(masters_human, 2):
#     print len(utils_graph.intersectLists([p1[0],p2[1]])), len(utils_graph.intersectLists([p2[0],p1[1]]))
    
     



# for protein in human_cons:
#      mk_host_plot(dir, 'human ' + protein, human_cons[protein], chicken_cons[protein], 
#                   limit,
#                   human_host_freqs, chicken_host_freqs,
#                   interesting[protein])
#      mk_host_plot(dir, 'chicken ' + protein, chicken_cons[protein], human_cons[protein], 
#                   limit,
#                   chicken_host_freqs, human_host_freqs,
#                   (interesting[protein][1], interesting[protein][0]))





