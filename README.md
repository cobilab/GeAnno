# GeAnno
*A supervised machine learning pipeline for ab initio gene finding*

## Table of Contents
- [Overview](#overview)
- [Quickstart](#quickstart)
- [Installation](#installation)
  - [Set up environment variable](#set-up-environment-variable)
  - [Install external tools](#install-external-tools)
  - [Install Python dependencies](#install-python-dependencies)
- [Running GeAnno](#running-geanno)
- [Pre-trained models](#pre-trained-models)
- [Training new models](#training-new-models)
  - [Option 1 - Using training scripts](#option-1---using-training-scripts-recommended-for-most-users)
  - [Option 2 - Using the experimental setup](#option-2---using-the-experimental-setup)
- [Repository structure](#repository-structure)

## Overview

**GeAnno** is a framework for identifying genomic regions of interest (genic vs. intergenic) using supervised machine learning.
It combines biologically motivated features with XGBoost to build accurate gene prediction tools.

## Quickstart

```bash
git clone https://github.com/Brums21/GeAnno.git
cd GeAnno

echo "export PLANT_DIR=$(pwd)" >> ~/.bashrc
source ~/.bashrc

chmod +x ./setup.sh
./setup.sh -t

python3 src/geanno.py -d example/dna/a_thaliana.fa -m models/models_genic_a_thaliana/model_undersampling_XGBoost_50.pkl
```

## Installation

### **Set up environment variable** 

```bash
echo "export PLANT_DIR=$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

This variable is used to locate datasets, binaries and outputs, and is needed for the normal functioning of this program.

### **Install external tools**

Run the setup script:

```bash
chmod +x setup.sh
```

- To see available options:

```bash
./setup.sh -h
```

- To install required dependencies (EMBOSS, JARVIS3, gto, etc.) and create a Python virtual environment where requirements are installed:
```bash
./setup.sh -t
```

- To also compile the C++ binaries:
```bash
./setup.sh -a
```
>Pre-compiled C++ binaries are already provided in bin/, so in most cases only -t is needed.

### **Install Python dependencies**

```bash
pip install -r requirements.txt
```

## **Running GeAnno**

### Running prediction on a DNA sequence:

```bash
python3 src/geanno.py -d <DNA_FILE> -m <MODEL>
```

You can use the following example and check if everything is working as intended:
```bash
python3 src/geanno.py -d example/dna/a_thaliana.fa -m models/models_genic_a_thaliana/model_undersampling_XGBoost_50.pkl
```

To see all options:
```bash
python3 src/geanno.py -h
```

## Pre-trained models

We provide all pre-trained models users can use to make their predictions.

- `models/models_genic_m_esculenta_PCA/model_undersampling_XGBoost_50.pkl`: trained on *M. esculenta* variants (undersampling, XGBoost), using PCA. Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_m_esculenta_1500_50_PCA.ipynb`

- `models/models_genic_a_thaliana/model_undersampling_XGBoost_50.pkl`: trained on *A. thaliana* species (undersampling, XGBoost). Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_a_thaliana_1500_50.ipynb`
- `models/models_genic_a_thaliana_PCA/model_undersampling_XGBoost_50.pkl`: trained on *A. thaliana* species (undersampling, XGBoost), using PCA. Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_a_thaliana_1500_50_PCA.ipynb`

- `models/models_genic_o_sativa/model_undersampling_XGBoost_50.pkl`: trained on *O. sativa* species (undersampling, XGBoost). Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_o_sativa_1500_50.ipynb`
- `models/models_genic_o_sativa_PCA/model_undersampling_XGBoost_50.pkl`: trained on *O. sativa* species (undersampling, XGBoost), using PCA. Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_o_sativa_1500_50_PCA.ipynb`

- `models/models_genic_genemark/model_undersampling_XGBoost_50.pkl`: trained on *D. melanogaster*, *C. elegans* and *A. thaliana* species (undersampling, XGBoost). Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_genemark_1500_50.ipynb`
- `models/models_genic_genemark_PCA/model_undersampling_XGBoost_50.pkl`: trained on *D. melanogaster*, *C. elegans* and *A. thaliana* species (undersampling, XGBoost), using PCA. Corresponds to the output of `experimental_setup/notebooks/classifiers/classifiers_genemark_1500_50_PCA.ipynb`

## **Training new models**

We offer two main ways to generate models in GeAnno:

### **Option 1** - Using training scripts (recommended for most users):

**1. Generate dataset with features**

From FASTA + GFF3 files downloaded in advance:

```bash
chmod +x train_models/generate_dataset_model.sh
train_models/generate_dataset_model.sh -f <FASTA_DIR> -d <GFF_DIR> -o dataset.csv
```

> **Note:** The annotation and corresponding DNA files should have the same basename.
> For example, `example.fa` must have a matching `example.gff3`.

You can obtain the list of parameters of this script by running 

```bash
train_models/generate_dataset_model.sh -h
```

**2. Train models**

Train Random Forest, XGBoost, or both:

```bash
python train_models/train_model.py -d dataset.csv -m RF XGBoost
```
Options include:
- `-cv`: Number of CV folds (default: 3)
- `-n_iter`: Iterations for random search (default: 10)
- `-n_jobs`: Parallel jobs for training (default: 4)
- `-max_workers`: Max parallel models trained (default: 1)
- `-ov`: Use oversampling instead of undersampling

Example (train only XGBoost with oversampling and 20 iterations):

```bash
python train_models/train_model.py -d dataset.csv -m XGBoost -ov --n_iter 20
```

### **Option 2** - Using the experimental setup:

To reproduce the experiments described in the notebooks, you can use the experimental setup workflow:

**1. Download necessary datasets:**

```bash
experimental_setup/dataset_scripts/download_datasets.sh -A
```

**2. Extract features**

```bash
experimental_setup/dataset_scripts/extract_features.sh
```

**3. Run the provided notebooks (in `experimental_setup/notebooks/`) to train and evaluate models.**

- Example: `classifiers/classifiers_a_thaliana_1500_50.ipynb` reproduces the *A. thaliana* model included in `models/models_genic_a_thaliana/model_undersampling_XGBoost_50.pkl`.

> A more detailed description of these scripts is available in [experimental_setup/README.md](experimental_setup/README.md).


## Repository structure

```text
GeAnno/
├── bin/                  # C++ binaries and installed tools
├── example/              # Example FASTA, GFF3, and model
├── experimental_setup/   # Scripts and notebooks for experiments
├── models/               # Trained models (pre-trained included)
├── src/                  # Python and C++ source code
└── train_models/         # Dataset generation and training scripts
```
