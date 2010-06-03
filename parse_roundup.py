""" For clusters with all species,
    print clusterID, species, geneIdentifier.
    A species can have more than one ID in a
    cluster

    Create Homo_Mus_Mac_Pan_Rat_5 table in 
    database roundup
"""
import sys, utils
from collections import defaultdict

roundup_file = sys.argv[1]
species = int(sys.argv[2])
outfile = sys.argv[3]

cluster = defaultdict(dict)
current_cluster = ''
with open(outfile, 'w') as outf:
    with open(roundup_file) as f:
        for line in f:
            if 'Gene Cluster' in line:
                if len(cluster.keys()) == species:
                    for a_species in cluster:
                        for ID in cluster[a_species]:
                            outf.write('%s\t%s\t%s\n' %
                                       (current_cluster, a_species, ID))
                current_cluster = line.split()[2][1:]
                cluster = defaultdict(dict)
            elif 'Genome' not in line and line.strip() != '' and current_cluster:
                sp = line.split('\t')
                if sp[0].strip() != '-':
                    cluster[sp[1]][sp[0]] = True
    # catch the last one
    if len(cluster.keys()) == species:
        for a_species in cluster:
            for ID in cluster[a_species]:
                outf.write('%s\t%s\t%s\n' %
                           (current_cluster, a_species, ID))

# put in database
(conn, cur) = utils.init_mysql('roundup')
table_name = 'Homo_Mus_Mac_Pan_Rat_5'
line = "CREATE TABLE " + table_name + "(cluster_id INT, species CHAR(30), seq_id CHAR(50))"
cur.execute(line)
line = "LOAD DATA LOCAL INFILE '" + outfile + "' INTO TABLE " + table_name + " FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n'"
cur.execute(line)
conn.commit()
cur.close()
conn.close()
