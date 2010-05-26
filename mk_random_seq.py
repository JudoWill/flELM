import markov_chain, utils_fasta

chain = markov_chain.MarkovChain()
f = utils_fasta.loadFASTA('/home/perry/bioperry/Projects/Thesis/Data/FASTA/Human/hprd_new.intr.fasta')
for k in f:
    chain.add(f[k])
for i in range(10):
    print "".join(chain.random_output())
