#!/usr/bin/env bash
set -euo pipefail

: "${PLANT_DIR:?PLANT_DIR not set!}"

Help()
{
    echo "Download datasets for genome annotation experiments."
    echo
    echo "Syntax: script.sh [-h|a|o|m|g|A]"
    echo "options:"
    echo "h     Print this Help."
    echo "a     Downloads Arabidopsis thaliana dataset"
    echo "o     Downloads Oryza sativa dataset"
    echo "m     Downloads Manihot esculenta variants dataset"
    echo "g     Downloads Genemark test set dataset"
    echo "A     Downloads all datasets"
    echo
}

download_file() {
    local url="$1"
    local dest="$2"

    if [ ! -f "$dest" ]; then
        echo "Downloading $(basename "$dest")..."
        wget -q "$url" -O "$dest"
        # If compressed file, decompress it
        if [[ "$dest" == *.gz ]]; then
            gunzip -f "$dest"
        fi
    else
        echo "$(basename "$dest") already exists. Skipping..."
    fi
}

download_dataset() {
    local species="$1"
    local dna_dir="${PLANT_DIR}/experimental_setup/dataset/${species}/dna"
    local annot_dir="${PLANT_DIR}/datasets/${species}/annotation"
    shift

    mkdir -p "$dna_dir" "$annot_dir"

    declare -n dna_files="$1"
    declare -n gff_files="$2"

    echo "Downloading DNA for $species"
    for file in "${!dna_files[@]}"; do
        download_file "${dna_files[$file]}" "$dna_dir/$file"
    done

    echo "Downloading Annotations for $species"
    for file in "${!gff_files[@]}"; do
        download_file "${gff_files[$file]}" "$annot_dir/$file"
    done
}

# Manihot esculenta variants
declare -A DNA_m_esculenta=(
    ["AM560_V8.fa"]="http://www.cassava-mdb.com/static/download/genome_download/AM560_V8.1.fa"
    ["Baxi.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/Baxi.genome.fa"
    ["TME204.fa"]="http://www.cassava-mdb.com/static/download/genome_download/tme204.fasta"
    ["W14.fa"]="http://www.cassava-mdb.com/static/download/genome_download/w14_pseudo_10062018.fa"
    ["XX048.fa"]="http://www.cassava-mdb.com/static/download/genome_download/xx048_2n_review.fasta"
    ["a4011.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4011.genome.fasta"
    ["a4041.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4041.genome.fa"
    ["a4047.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4047.genome.fa"
    ["a4111.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4111.genome.fasta"
    ["a4144.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4144.genome.fasta"
    ["a4160.fa"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4160.genome.fasta"
)
declare -A GFF_m_esculenta=(
    ["AM560_V8.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/AM560_V8.1.gene.gff3"
    ["Baxi.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/Baxi.gff"
    ["TME204.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/tme204.gff3"
    ["W14.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/w14Mes.gff3"
    ["XX048.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/xx048.gene.gff"
    ["a4011.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4011.gff"
    ["a4041.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4041.gff"
    ["a4047.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4047.gff"
    ["a4111.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4111.gff"
    ["a4144.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4144.gff"
    ["a4160.gff3"]="http://www.cassava-mdb.com/static/download/genome_download/mp_new/a4160.gff"
)

# Oryza sativa
declare -A DNA_o_sativa=(
    ["o_sativa.fa"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
)
declare -A GFF_o_sativa=(
    ["o_sativa.gff3"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.62.gff3.gz"
)

# Arabidopsis thaliana
declare -A DNA_a_thaliana=(
    ["a_thaliana.fa"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
)
declare -A GFF_a_thaliana=(
    ["a_thaliana.gff3"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz"
)

#Genemark
declare -A DNA_genemark=(
    ["a_thaliana.fa"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    ["c_elegans.fa"]="https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
    ["d_melanogaster.fa"]="https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz"
)
declare -A GFF_genemark=(
    ["a_thaliana.gff3"]="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz"
    ["c_elegans.gff3"]="https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz"
    ["d_melanogaster.gff3"]="https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gff3.gz"
)



while getopts "haomgA" option; do
   case $option in
        h) Help; exit;;
        a) download_dataset "dataset_a_thaliana" DNA_a_thaliana GFF_a_thaliana ;;
        o) download_dataset "dataset_o_sativa" DNA_o_sativa GFF_o_sativa ;;
        m) download_dataset "dataset_m_esculenta" DNA_m_esculenta GFF_m_esculenta ;;
        g) download_dataset "dataset_genemark" DNA_genemark GFF_genemark ;;
        A)
            download_dataset "dataset_a_thaliana" DNA_a_thaliana GFF_a_thaliana
            download_dataset "dataset_o_sativa" DNA_o_sativa GFF_o_sativa
            download_dataset "dataset_m_esculenta" DNA_m_esculenta GFF_m_esculenta
            download_dataset "dataset_genemark" DNA_genemark GFF_genemark
            ;;
        \?) echo "Error: Invalid option."; echo "Use $0 -h for help."; exit 1;;
   esac
done

if [ $OPTIND -eq 1 ]; then
    echo "Error: No flag provided."; echo "Use $0 -h for help."; exit 1;;
fi