""" For clusters with all species,
    print clusterID, species, geneIdentifier.
    A species can have more than one ID in a
    cluster

    The save TEXT option failed me on roundup,
    so I copied the cluster view table, and that
    is what is given to this script as input.

    Create roundup10 table in 
    database roundup
"""
import sys, utils
from collections import defaultdict

roundup_file = sys.argv[1] # Homo_Mus_Pan_Rat_Bos_Can_Gal_Tae_Dan_Mac.roundup
species_count = int(sys.argv[2]) # 10
outfile = sys.argv[3]

with open(outfile, 'w') as outf:
    with open(roundup_file) as f:
        species_ls = f.readline().strip().split('\t')[5:]
        for line in f:
            found = {}
            for species, ID in zip(species_ls, [x.strip() for x in line.split('\t')[5:]]):
                if ID:
                    if ',' in ID:
                        found[species] = [x.strip() for x in ID.split(",")]
                    else:
                        found[species] = [ID]
            if len(found.keys()) == species_count:
                cluster = line.split('\t')[0]
                for species in found:
                    for ID in found[species]:
                        outf.write('%s\t%s\t%s\n' %
                                   (cluster, species, ID))

# put in database
# (conn, cur) = utils.init_mysql('roundup')
# table_name = 'roundup10'
# line = "CREATE TABLE " + table_name + "(cluster_id INT, species CHAR(30), seq_id CHAR(50))"
# cur.execute(line)
# line = "LOAD DATA LOCAL INFILE '" + outfile + "' INTO TABLE " + table_name + " FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n'"
# cur.execute(line)
# conn.commit()
# cur.close()
# conn.close()
