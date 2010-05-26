import random

residues = ['A','C','E','D','G','F','I',
            'H','K','M','L','N','Q','P',
            'S','R','T','W','V','Y']

def draw():
    return residues[random.randint(0,19)]

files = []
for x in xrange(10):
    files.append(open('elm_exp_random' 
                      + str(x), 'w'))
with open('elm_expressions.txt') as f:
    for line in f:
        [elm, exp] = line.strip().split('\t')
        for afile in files:
            new_exp = ''
            for chr in exp:
                if chr in residues:
                    new_exp += draw()
                else:
                    new_exp += chr
            afile.write(elm + '\t' + new_exp + '\n')
for f in files:
    f.close()
                
        
