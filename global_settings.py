TEST_GENOMES = ('H_sapiens', 'M_musculus', 'Pan_troglodytes', 
                'R_norvegicus', 'Gallus_gallus', 'Taeniopygia_guttata',
                'Canis_familiaris', 'Bos_taurus')
GENOMES = ('M_musculus', 'Bos_taurus','Canis_familiaris','D_rerio',
           'Equus_caballus', 'Gallus_gallus', 'H_sapiens',
           'Macaca_mulatta', 'R_norvegicus', 'Sus_scrofa',
           'Taeniopygia_guttata', 'Pan_troglodytes')
ALIASES = {'M_musculus':'Mouse', 'Bos_taurus':'Cow','Pan_troglodytes':'Chimp',
           'Canis_familiaris':'Cat','D_rerio':'Fish',
           'Equus_caballus':'Horse', 'Gallus_gallus':'Chicken', 
           'H_sapiens':'Human', 'Sus_scrofa':'Swine',
           'Macaca_mulatta':'Monkey', 'R_norvegicus':'Rat',
           'Taeniopygia_guttata':'finch',
           'swine':'swineFlu', 'human':'humanFlu', 
           'chicken':'chickenFlu',
           'web':'web', 'regex':'regex'}

FLU_NAMES = {'duck':('duck', 'mallard'), 'chicken':('chicken',), 
             'goose':('goose', 'swan'), 'swine':('swine', 'pig'), 
             'equine':('equine', 'horse'), 'turkey':('turkey',), 
             'human':('human',)}

PROTEIN_ALIAS = {'nuclear export protein':'nonstructural protein 2',
                 'nucleoprotein':'nucleocapsid protein',
                 'haemagglutinin':'hemagglutinin',
                 'M2 protein':'matrix protein 2',
                 'matrix protein M1':'matrix protein 1',
                 'M1 protein':'matrix protein 1'}
VIRUS_PROTEINS = ['hemagglutinin', 'nonstructural protein 1', 'nonstructural protein 2',
                  'matrix protein 1', 'neuraminidase', 'matrix protein 2',
                  'nucleocapsid protein', 'polymerase PB1', 'polymerase PA', 
                  'polymerase PB2', 'PB1-F2 protein']

SEQ_FRAC_CUT = '.05'
