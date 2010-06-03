""" Create a table of FASTA sequences.
    ID -> sequence
"""
import utils, sys, random, os

fasta_file = sys.argv[1]
table_name = sys.argv[2]

max_len = 0
tmp = 'tmp' + str(random.randint(0,100))
with open(tmp, 'w') as f:
    for protein, seq in utils.fasta_iter(fasta_file, lambda line: line.split('|')[1]):
        f.write('%s\t%s\n' %
                (protein, seq))
        max_len = max([len(seq), max_len])
    
(conn, cur) = utils.init_mysql('fasta')
line = "CREATE TABLE " + table_name + " ( seq_id CHAR(100), seq TEXT(" + str(max_len+100) + ") )"
cur.execute(line)
line = "LOAD DATA LOCAL INFILE '" + tmp + "' INTO TABLE " + table_name + " FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n'"
cur.execute(line)
conn.commit()
cur.close()
conn.close()

os.system('rm ' + tmp)


