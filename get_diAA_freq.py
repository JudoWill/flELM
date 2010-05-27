import sys, utils

def check(c):
    return c != 'X' and c != 'U' and seq.find('B') == -1 and seq.find('J') == -1 and seq.find('Z') == -1

def countTransitions(line, counts):
    for i in xrange(len(line)-1):
        c1 = line[i]
        c2 = line[i+1]
        if check(c1) and check(c2):
            if c1 not in counts:
                counts[c1] = {}
            if c2 not in counts[c1]:
                counts[c1][c2] = 0
            counts[c1][c2] += 1

counts = {}
seq = ''
with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == '>':
            if seq:
                countTransitions(seq, counts)
            seq = ''
        else:
            seq += line.strip()
            
for c1 in counts:
    total = 0
    for c2 in counts[c1]:
        total += counts[c1][c2]
    for c2 in counts[c1]:
        print c1 + c2 + '\t' + str(float(counts[c1][c2])/float(total))
