""" Take the parsed roundup file, grab
    the monkey IDs, and use biomaRt to
    get protein sequences for them.
"""
import sys, os, utils, random

parsed_roundup_file = sys.argv[1]
output_file = sys.argv[2]

genes = {}
with open(parsed_roundup_file) as f:
    for line in f:
        sp = line.strip().split('\t')
        if sp[1] == 'Macaca mulatta':
            genes[sp[-1]] = True

# biomaRt appends to the output
# fasta file,
# so it must be removed
if os.path.exists(output_file):
    os.system('rm '
              + output_file)

chunks = utils.chunks(genes.keys(), 200)
rfiles = []
outfiles = []
file_name_prefix = 'r_input' + str(random.randint(0,100))
for count, a_chunk in enumerate(chunks):
    file_name = file_name_prefix + str(count)
    file_name_out = file_name + '.out'
    with open(file_name, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write("ls<-c('"
                + "','".join([str(x) for x in a_chunk])
                + "')\n")
        f.write("get_monkey(ls,"
                + "'" + file_name_out + "')\n")
        f.write("q()\n")
    os.system('R < ' + file_name + ' --no-save')
    rfiles.append(file_name)
    outfiles.append(file_name_out)

file_line = ' '.join(outfiles)
# remove empty lines
os.system('cat ' + file_line
          + " | grep -v '^$' > "
          + output_file)

# clean up
os.system('rm ' + ' '.join(files) + ' '
          + file_line)
