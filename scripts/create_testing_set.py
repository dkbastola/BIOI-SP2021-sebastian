#!/usr/bin/python3
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

species_dir = 'testing/'
codon_dict = {}
t_file = open('testing_set.csv', 'w')


# Return DNA seq in fna file as a string
def get_dna(fna_path):
    seq = ''

    with open(fna_path, 'r') as f:

        header = f.readline().split('>')[1].split(',')[0].strip('\n')
        header = header.split(' ', 1)[1].rsplit(' ',1)[0]
        for line in f:
            seq += line.strip('\n')

    return seq, header 

def get_CDS(gb_path):

    record = SeqIO.read(gb_path, 'genbank')

    CDS_list = [f for f in record.features if f.type == "CDS"]
    return CDS_list

def create_codon_dict(cd):

    # Read file with codons
    with open('codonDic.txt', 'r') as f:
        codons = f.read().strip('\n').split(',')

    # Add codons to a dict
    for codon in codons:
        cd.update({codon: 0})

    return cd


def count_codons(cds_list, seq, cd):
    for cds in cds_list:
        temp_seq = cds.extract(seq)
        
        invalid_nts = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']

        # Check if the seq nts are in invalid_nts
        valid_nts = all(nt not in temp_seq for nt in invalid_nts)
        if valid_nts:
            # iterate through every codon

            for n in range(0, len(temp_seq) - len(temp_seq) % 3, 3):
                codon = temp_seq[n: n + 3]
                # count the codon
                cd.update({codon: cd[codon] + 1})
                print(cd)
        else:
            print(temp_seq)
            print(len(temp_seq) % 3)

        
    return cd

# Create codon dict
codon_dict = create_codon_dict(codon_dict)

# file header
header = 'species,' + ','.join(list(codon_dict.keys())) + '\n'
t_file.write(header)

for species in os.listdir(species_dir):

    current_dir = species_dir + species + '/ncbi_dataset/data/'

    # Check if current dir is a directory
    if os.path.isdir(current_dir):

        # Look for GCF
        for f in os.listdir(current_dir):
            if f.startswith('GCF'):
                Mtfna_path = ''
                annot_path = ''

                # Retrieve genome fna and documentation
                for Mt in os.listdir(current_dir + f):
                    if Mt.startswith('chrMT'):
                        Mtfna_path = current_dir + f + '/' + Mt
                        

                    if Mt.startswith('NC_'):
                        annot_path = current_dir + f + '/' + Mt
                        

                    if Mtfna_path and annot_path:
                        print(annot_path)

                        # Get dna seq and record name
                        fna, rname = get_dna(Mtfna_path)

                        # Get list of all CDS coordinates
                        CDSList = get_CDS(annot_path)

                        # count codons and write results
                        codon_dict = count_codons(CDSList, fna, codon_dict)

                        res_counts = rname + ',' + ','.join(list(map(str, codon_dict.values()))) + '\n'
                        t_file.write(res_counts)

                        # reset all dictionary counts to 0
                        codon_dict = dict.fromkeys(codon_dict, 0)



                        



                        




