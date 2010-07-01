FLU_PROTEINS = ('hemagglutinin', 'neuraminidase', 'nucleocapsid protein',
                'matrix protein 1', 'nonstructural protein 1', 
                'matrix protein 2',
                'nonstructural protein 2', 'polymerase PA', 'polymerase PB2',
                'polymerase PB1', 'PB1-F2 protein')

TEST_GENOMES = ('H_sapiens', 'M_musculus', 'Pan_troglodytes', 
                'Gallus_gallus', 'Taeniopygia_guttata',
                'Sus_scrofa', 'Equus_caballus')

TEST_GENOMES2 = ('H_sapiens', 'M_musculus', 'Pan_troglodytes', 
                'Gallus_gallus', 'Taeniopygia_guttata',
                'Bos_taurus', 'Canis_familiaris')

MAMMALS = ('Homo sapiens', 'Mus musculus', 'Pan troglodytes', 
           'Sus scrofa', 'Equus caballus')
MAMMALS2 = ('H_sapiens', 'M_musculus', 'Macaca_mulatta', 
           'Sus_scrofa', 'Equus_caballus')

MAMMALS3 = ('H_sapiens', 'M_musculus', 'Macaca_mulatta', 
           'Sus_scrofa', 'Equus_caballus', 'Gallus_gallus')

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

AA_SUB_8 = {'L':'l', 'V':'l', 'I':'l', 'M':'l', 'C':'l',
          'A':'a', 'G':'a',
          'S':'s', 'T':'s',
          'P':'P',
          'F':'f', 'Y':'f', 'W':'f',
          'E':'e', 'D':'e', 'N':'e', 'Q':'e',
          'K':'k', 'R':'k',
          'H':'H'}

AA_SUB_4 = {'L':'l', 'V':'l', 'I':'l', 'M':'l', 'C':'l',
          'A':'a', 'G':'a', 'S':'a', 'T':'a', 'P':'a',
          'F':'f', 'Y':'f', 'W':'f',
          'E':'e', 'D':'e', 'N':'e', 'Q':'e', 'K':'e', 'R':'e', 'H':'e'}

AA_SUB_2 = {'L':'l', 'V':'l', 'I':'l', 'M':'l', 'C':'l', 'A':'l', 'G':'l', 'S':'l', 'T':'l', 'P':'l', 'F':'l', 'Y':'l', 'W':'l',
            'E':'e', 'D':'e', 'N':'e', 'Q':'e', 'K':'e', 'R':'e', 'H':'e'}

SEQ_FRAC_CUT = '.05'
