import utils
fasta = utils.loadFASTA('human.H3N2.fa')
for key in fasta:
    if key.split('.')[1] == 'hemagglutinin':
        print fasta[key][3], fasta[key][121:125]
