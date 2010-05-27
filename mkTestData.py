features = {'A':(1,1,1),
            'B':(1,1,1.5),
            'C':(0,1,2),
            'D':(0,1.5,2),
            'E':(2,1,.5)}
with open('test.data', 'w') as f:
#    f.write('Name\tVal1\tVal2\tVal3\n')
    for k in features:
        f.write('\t'.join(str(x) for x in features[k]) + '\n')
with open('labels', 'w') as f:
    f.write('A\tB\tC\tD\tE\n')
