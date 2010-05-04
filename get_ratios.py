from collections import defaultdict
import sys, utils

classes = defaultdict(dict)
elm2class = {}
with open('elm_classes') as f:
    for line in f:
        elm = line.strip()
        sp = elm.split('_')
        a_class = sp[0] + '_' + sp[1]
        classes[a_class][elm] = True
        elm2class[elm] = a_class

member_counts = {}
for a_class in classes:
    member_counts[a_class] = defaultdict(utils.init_zero)
class_counts = defaultdict(utils.init_zero)
with open(sys.argv[1]) as f:
    for line in f:
        elm, seq, count_st, freq = line.strip().split('\t')
        
        sp = elm.split('_')
        a_class = sp[0] + '_' + sp[1]
        if a_class in classes:
            count = int(count_st)
            class_counts[a_class] += count
            member_counts[a_class][elm] += count
for a_class in classes:
    for elm in classes[a_class]:
        if elm in member_counts[a_class]:
            val = float(member_counts[a_class][elm])/float(class_counts[a_class])
        else:
            val = float(0)
        print(a_class + '\t' + elm + '\t'
              + str(val))
                
