import os
xml = ''
with open('monkey_biomart.xml') as f:
    for line in f:
        xml += line.strip()
line = "wget -O results.txt 'http://www.biomart.org/biomart/martservice?query=" + xml + "'"
print line
os.system(line)
