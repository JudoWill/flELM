"""Use biomart martservice to get protein sequence data.
   IDs must be in ensmbl protein ID form & you need to know
   your species database, which you can find here
   http://www.biomart.org/biomart/martview/"""
import sys, os, utils, random

parsed_roundup_file = sys.argv[1]
species = sys.argv[2]
dataset = sys.argv[3] # mmulatta_gene_ensembl or drerio_gene_ensembl
output_file = sys.argv[4]

query_prefix= "'http://www.biomart.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"FASTA\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"" + dataset + "\" interface = \"default\" ><Filter name = \"ensembl_peptide_id\" value = \""

query_suffix = "\"/><Attribute name = \"peptide\" /><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"ensembl_peptide_id\" /></Dataset></Query>'"

genes = {}
with open(parsed_roundup_file) as f:
    for line in f:
        if species in line:
            genes[line.strip().split('\t')[-1]] = True
chunks = utils.chunks(genes.keys(), 200)
files = []
file_name_prefix = 'chunk' + str(random.randint(0,100))
for count, a_chunk in enumerate(chunks):
    file_name = file_name_prefix + str(count)
    os.system("wget --tries=10 "
              + query_prefix
              + ','.join(a_chunk) 
              + query_suffix 
              + " --output-document="
              + file_name)
    files.append(file_name)
    if len(files) == 2:
        break

file_line = ' '.join(files)
# make one large file from tmp files
# remove empty lines w/ grep
os.system('cat ' + file_line
          + " | grep -v '^$' > "
          + output_file)
# clean up
#os.system('rm ' + file_line)
