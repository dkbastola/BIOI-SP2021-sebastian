#!/usr/bin/python3
from Bio import Entrez

import os
import time

APIkey = '21f9bebb6c243e6f7cdb3e38e088012b0709'
email = 'se.espinoza132@gmail.com'


# Mandatory parameters
Entrez.api_key = APIkey
Entrez.email = email

# Search ncbi for genome
search_terms = ['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria',
                'epsilonproteobacteria', 'zetaproteobacteria', 'oligoflexia', 
                'acidithiobacillia', 'hydrogenophilalia']

gamma_terms = ['Acidiferrobacterales', 'Aeromonadales', 'Alteromonadales',
               'Arenicellales', 'Cardiobacteriales', 'Cellvibrionales',
               'Chromatiales', 'Enterobacterales', 'Immundisolibacterales',
               'Legionellales', 'Methylococcales', 'Nevskiales',
               'Oceanospirillales', 'Orbales', 'Pasteurellales',
               'Pseudomonadales', 'Salinisphaerales', 'Thiotrichales',
               'Vibrionales', 'Xanthomonadales']

# Search for all proteobacteria and download genomes
def get_terms(term_list):
    for term in term_list:

        # Get the taxonomy ids for all proteobacteria
        handle = Entrez.esearch(db='Taxonomy', term=term)
        record = Entrez.read(handle)

        temp_id = record['IdList'][0]
        fname = term + '.zip'
    
        cmd = './datasets download genome taxon ' + temp_id + ' --refseq --reference --filename ' + fname
    
        os.system(cmd)


# Unzip all proteobacteria data
def unzip_all(directory):

    cmd = 'unzip '

    for filename in os.listdir(directory):

        # Unzip all .zip files
        if filename.endswith('.zip'):

            temp_dir = filename.split('.')[0]

            os.system(cmd + filename + ' -d ' + temp_dir)


# Search for all other proteobacteria
get_terms(search_terms)

# Search for all gammaproteobacteria
get_terms(gamma_terms)

unzip_all('.')
