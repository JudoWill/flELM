""" Given a species in the parsed roundup file,
    get protein sequences for each GenInfo ID
    using NCBI's eutils.
"""
import sys, os, utils, random

parsed_roundup_file = sys.argv[1]
species = sys.argv[2]
output_file = sys.argv[3]

query='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id='

genes = {}
with open(parsed_roundup_file) as f:
    for line in f:
        if species in line:
            genes[line.strip().split('\t')[-1]] = True
chunks = utils.chunks(genes.keys(), 500)
files = []
file_name_prefix = 'chunk' + str(random.randint(0,100))
for count, a_chunk in enumerate(chunks):
    file_name = file_name_prefix + str(count)
    os.system("wget --tries=10 '"
              + query
              + ','.join(a_chunk) 
              + "' --output-document="
              + file_name)
    files.append(file_name)

file_line = ' '.join(files)
# make one large file from tmp files
# remove empty lines w/ grep
os.system('cat ' + file_line
          + " | grep -v '^$' > "
          + output_file)
# clean up
os.system('rm ' + file_line)
