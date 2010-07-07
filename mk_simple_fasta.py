"""Convert raw sequences to simplified versions
   for use in scanning with simple patterns
   (those made with utils.mk_sub)"""
import sys, global_settings, os, utils

in_dir = sys.argv[1]
out_dir = sys.argv[2]

for g in global_settings.TEST_GENOMES:
    afile = os.path.join(in_dir, g + '.fa')
    ofile = os.path.join(out_dir, g + '.fa')
    with open(ofile, 'w') as f:
        for ID, seq in utils.fasta_iter(afile):
            f.write('>' + ID + '\n')
            f.write(utils.mk_sub(seq) + '\n')
