# Experimental setup

This directory contains the notebooks and scripts used to generate the datasets employed in the comparison with existing genome annotation tools in the references thesis.

To prepare the datasets, two steps are required:

---

### 1. Download the datasets 

Run:
```bash
$ dataset_scripts/download_datasets.sh -A
```

This will download all required **FASTA** and **GFF3** files directly from the [Ensembl](https://www.ensembl.org/) and [Ensembl Plants](https://plants.ensembl.org/) genome browsers. 

The files are saved under `experimental_setup/datasets/`, with the following directories created:

- `dataset_m_esculenta`
- `dataset_o_sativa`
- `dataset_a_thaliana`
- `dataset_genemark`

> **Note**: the M. esculenta dataset is particularly large, which makes both the download and subsequent feature extraction slow. If you only want the other datasets, we recommend running:
```bash
$ dataset_scripts/download_datasets.sh -a -o -g
```

You can also display the full list of available options and explanations with:
```bash
$ dataset_scripts/download_datasets.sh -h
```

### 2. Extract features

```bash
$ dataset_scripts/extract_features.sh
```

This script processes all datasets found in `experimental_setup/datasets/` and extracts their corresponding features for downstream analysis.




