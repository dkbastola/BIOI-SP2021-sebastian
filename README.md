# bioiSeniorProject SPRING 2021

This repository contains all of the scripts used to generate datasets, train machine learning models, and test models used to analyze codon bias in mitochondria.  
<br/><br/>
<em>Codon counts are calculated by finding all CDS coordinates using .gff files, extracting those sequences within genomic FASTA files, and dividing each CDS sequence into codons.  
All codon counts are represented as .csv files for compatibility with machine learning algorithms.  
The classifiers built in this study will assign one of four proteobacterial classes: alpha, beta, delta, and epsilon. </em>

Files description:
* `scripts/get_genomes.py` - Retrieves all reference proteobacterial genomes from NCBI
* `scripts/create_training_set.py` - Creates proteobacterial codon count .csv file for training the models
* `scripts/create_naive_model.py` - Trains and tests the Naive Bayes model
* `scripts/create_xgboost_model.py` - Trains and tests the XGBoost model
* `scripts/create_neural_model.py` - Trains and create the Artificial Neural Network model (Multilayer Perceptron)
* `scripts/datasets` - Executable file from NCBI to download bulk datasets (used by `get_genomes.py` and `get_mito.py`)
* `scripts/get_mito.py` - Retrieves 37 reference mitochondrial genomes from NCBI
* `scripts/species_list.txt` - Contains all 37 mitochondrial species names
* `scripts/codonDict.txt` - Text file containing all 64 codons (used by `create_training_set.py` and `create_testing_set.py`)  

* `datasets/training_set.csv` - Codon counts from all proteobacterial genomes retrieved by `get_genomes.py` used to train all 3 models.
* `datasets/testing_set.csv` - Codon counts from 37 mitochondrial genomes retrieved by `get_mito.py` used to test all 3 models.

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
* Download all the proteobacterial genomes from NCBI (ensure you have a stable internet connection)
* Tabularize all proteobacterial genomes by codon counts to create a training dataset
* Repeat the last two steps for 37 mitochondrial genomes to create a testing dataset
* Train the 3 models on the proteobacterial dataset and test on the mitochondrial dataset

To achieve the full analysis run the following files: 
<br/><br/> 
1. Download all genomes using NCBI's dataset tool: `./get_genomes.py`
* input: N/A
* output: unzipped directories containing .fna and .gff files of proteobacterial genomes  
**Note**: This process make take several hours and the files downloaded span more than 100+ GB.
<br/><br/>
2. Process all genomic files: `./create_training_set.py`
* input: .fna and .gff files fetched by `get_genomes.py`
* output: `training_set.csv` -- file that contains all codon counts as comma-separated values  
**Note**: This script may take 30 min or more to run.
<br/><br/>
3. Download all mitochondrial genomes listed in `species_list.txt`: `./get_mito.py`
* input: `species_list.txt`
* output: unzipped directories containing .fna and .gff files of mitochondrial genomes
<br/><br/>
4. Process all mitochondrial genomic files: `./create_testing_set.py`
* input: .fna and .gff files fetched by `get_mito.py`
* output: `testing_set.csv` -- file that contains all codon counts as comma-separated values
<br/><br/>
5. Train and test  machine learning models: `./create_naive_model.py`, `./create_xgboost_model.py`, and `./create_neural_model.py`
* input: `training_set.csv`
* output: Contingency table, model scores (precision, recall, etc...), mitochondrial predictions, and ROC plot.  
**Note**: The XGBoost model also displays variable weights.
<br/><br/>
# Alternative Running Method
Since steps 1-4 are time consuming, they be skipped above by running `cp datasets/testing_set.csv datasets/training_set.csv scripts/` or by moving both files into the scripts folder using other methods.  
After moving the training and testing sets into `scripts`, you can run step 5.


