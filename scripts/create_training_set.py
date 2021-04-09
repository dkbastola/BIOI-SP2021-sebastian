#!/usr/bin/python3

import gffutils
import pyfaidx

import os


species_dir = ['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria',
                'epsilonproteobacteria', 'zetaproteobacteria', 'oligoflexia',
                'acidithiobacillia', 'hydrogenophilalia']

path = '/home/sespinoza/genome_datasets/'
t_file = open('training_set.csv', 'w')
codon_dict = {}

def create_codon_dict(cd):

    # Read file with codons
    with open('codonDic.txt', 'r') as f:
        codons = f.read().strip('\n').split(',')
    
    # Add codons to a dict
    for codon in codons:
        cd.update({codon: 0})

    return cd

# Find a file with a certain extension
def find_file(p, pattern):
    for filename in os.listdir(p):
        if filename.endswith(pattern):
            return p + '/' + filename

def record_name(r):
    with open(r, 'r') as f:
        return f.read().strip('\n').split(' ')[1]



def count_codons(f, db, cd):
    for cds in db.features_of_type('CDS', order_by='start'):
        try:
            temp_seq = cds.sequence(f)
        except:
            temp_seq = ''
        
        invalid_nts = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']
        valid_nts = all(nt not in temp_seq for nt in invalid_nts)
        if not len(temp_seq) % 3 and valid_nts:
            # iterate through every codon

            for n in range(0, len(temp_seq), 3):
                codon = temp_seq[n: n + 3]
                # count the codon
                cd.update({codon: cd[codon] + 1})
        else:
            print(cds.source)
            print(cds.sequence(f))
            print(len(temp_seq))

        print(cd)
    return cd


# create codon dict
codon_dict = create_codon_dict(codon_dict)

# Output header
header = 'species,' + ','.join(list(codon_dict.keys())) + ',label' + '\n'
t_file.write(header)

for species in species_dir:
    for filename in os.listdir(species + '/ncbi_dataset/data'):
        
        gdir = path + species + '/ncbi_dataset/data/' + filename
        print(gdir)

        if os.path.isdir(gdir):

            # Find annotation and fasta files
            temp_gff = find_file(gdir, 'gff')
            temp_fna = find_file(gdir, 'fna')
            rname = record_name(temp_fna)

            # Build annotation database
            temp_db = gffutils.create_db(temp_gff, gdir + '/a.db', merge_strategy='merge', force=True)
            temp_fna = pyfaidx.Fasta(temp_fna)

            # count codons and write results
            codon_dict = count_codons(temp_fna, temp_db, codon_dict)
            
            res_counts = rname + ',' + ','.join(list(map(str, codon_dict.values()))) + ',' + species + '\n'
            t_file.write(res_counts)

            # reset all dictionary counts to 0
            codon_dict = dict.fromkeys(codon_dict, 0)


t_file.close()


