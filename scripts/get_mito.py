#!/usr/bin/python3
from Bio import Entrez

import os
import time

APIkey = '21f9bebb6c243e6f7cdb3e38e088012b0709'
email = 'se.espinoza132@gmail.com'


# Mandatory parameters
Entrez.api_key = APIkey
Entrez.email = email

species_dir = 'testing/'

def get_mito_genomes():
    with open('species_list.txt', 'r') as f:
        for line in f:

            # Remove newline character from species name
            species = line.strip('\n')

            handle = Entrez.esearch(db='Taxonomy', term=species)
            record = Entrez.read(handle)

            species_id = record['IdList'][0]
            fname = species_dir + species.replace(' ', '_') + '.zip'

            # Download reference Mt genome
            cmd = './datasets download genome taxon ' + species_id + ' --refseq --reference --chromosomes chrMT ' \
                '--exclude-gff3 --exclude-protein --exclude-rna' \
                ' --filename ' + fname
            
            os.system(cmd)

# Unzip all species genomes
def unzip_all(directory):

    cmd = 'unzip '

    for filename in os.listdir(directory):

        # Unzip all .zip files
        if filename.endswith('.zip'):

            temp_dir = directory + filename.split('.')[0]

            os.system(cmd + directory + filename + ' -d ' + temp_dir)

def download_gb(ddir, mt_fna):
    mt_path = ddir + '/' + mt_fna

    with open(mt_path, 'r') as f:
        header = f.readline()
        
        # Get the mitochondrial accession number
        acc_no = header.split(' ')[0].strip('>')

    # Download mito genome annotation
    handle = Entrez.efetch(db='nuccore', id=acc_no, rettype='gbwithparts', retmode='text')

    # write genome annotation to file
    with open(ddir + '/' + acc_no, 'w') as mf:

        print(acc_no)
        for line in handle:
            mf.write(line)


 

#get_mito_genomes()

#unzip_all(species_dir)
#for refseq in refList:
#    handle = Entrez.efetch(db='nuccore', id=refseq, rettype='gbwithparts', retmode='text')
#
#    with open(refseq, 'w') as f:
#        for line in handle:
#            f.write(line)

for species in os.listdir(species_dir):

    current_dir = species_dir + species + '/ncbi_dataset/data/' 
    
    # Check if current dir is a directory
    if os.path.isdir(current_dir):


        for f in os.listdir(current_dir):
            if f.startswith('GCF'):
                print(f)
                
                for Mt in os.listdir(current_dir + f):
                    if Mt.startswith('chrMT'):
                        download_gb(current_dir + f, Mt)


