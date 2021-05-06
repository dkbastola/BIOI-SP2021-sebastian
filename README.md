# bioiSeniorProject SPRING 2021

This repository contains all of the scripts used to generate datasets, train machine learning models, and test models used to analyze codon bias in mitochondria. 

Files description:
* `scripts/get_genomes.py`
* `scripts/create_training_set.py`
* `scripts/create_naive_model.py`
* `scripts/create_xgboost_model.py`
* `scripts/create_neural_model.py`
* `scripts/datasets`
* `scripts/get_mito.py`
* `scripts/species_list.txt`

# Requirements
* Python 3

Python libraries:
* numpy
* pandas
* matplotlib
* sklearn
* gffutils
* pyfaidx
* Bio

# Installation and setup
Make sure to have Python 3 installed on your machine, and then install all of the listed python libraries listed above. This can be achieved using pip3.

# Running Files
To recreate the full analysis described in the paper you must:
* Download all the proteobacterial genomes from NCBI
* Tabularize all proteobacterial genomes by codon counts
* Train the 3 models on the tabular dataset

To achieve the full analysis run the following files: 
<br/><br/> 
1. Download all genomes using NCBI's dataset tool: `./get_genomes.py`
* input: N/A
* output: unzipped directories containing .fna and .gff files for proteobacterial genomes  
**Note**: This process make take several hours and the files downloaded span more than 100+ GB.
<br/><br/>
2. Process all genomic files: `./create_training_set.py`
* input: .fna and .gff files fetched by `get_genomes.py`
* output: `training_set.csv` -- file that contains all codon counts as comma-separated values  
**Note**: This script may take 30 min or more to run.
<br/><br/>
3. Train machine learning models: `./create_naive_model.py`, `./create_xgboost_model.py`, and `./create_neural_model.py`
* input: `training_set.csv`
* output: Contingency table, model scores (precision, recall, etc...), mitochondrial predictions, and ROC plot.  
**Note**: The XGBoost model also displays variable weights.
<br/><br/>
# Alternative


