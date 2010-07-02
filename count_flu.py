"""Tool to count number of flu sequences"""

d = {}
with open('data/influenza.fa') as f:
    for line in f:
        if '>' in line:
            if 'chicken' in line:
                if 'Influenza A' in line:
                    protein = '_'.join(line.split('|')[-1].split('[')[0].strip().lower().split())
                    type = line.split('/')[-1].split('(')[-1].split(')')[0]
                    date = line.split('/')[-1].split('(')[0]
                    if protein not in d:
                        d[protein] = {}
                    if type not in d[protein]:
                        d[protein][type] = 0
                    #if date not in d[protein][type]:
                    #    d[protein][type][date] = 0
                    d[protein][type] += 1
for protein in d:
    for type in d[protein]:
        print protein + '\t' + type + '\t' + str(d[protein][type])
