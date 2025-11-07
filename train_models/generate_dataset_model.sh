#!/bin/bash

revcomp() {
    echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev
}

get_seeded_random(){
  seed="$1";
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null;
}

process_window() {
    local fasta="$1"
    local chrom="$2"
    local start="$3"
    local end="$4"
    local strand="$5"
    local label="$6"
    local species="$7"
    local output_result="$8"
    local final_dir="$9"

    region="${chrom}:${start}-${end}"
    seq=$(samtools faidx "$fasta" "$region" | grep -v '^>' | tr -d '\n')
    [[ "$strand" == "-" ]] && seq=$(revcomp "$seq")

    [[ -z "$seq" ]] && return

    tmp_id=$(uuidgen)
    seq_file="${final_dir}/tmp/seq_${tmp_id}.seq"
    fasta_file="${final_dir}/tmp/seq_${tmp_id}.fa"
    orf_file="${final_dir}/tmp/orf_${tmp_id}.fa"

    echo ">seq" > "$fasta_file"
    echo "$seq" >> "$fasta_file"
    echo "$seq" > "$seq_file"

    getorf -sequence "$fasta_file" -outseq "$orf_file" -reverse no > /dev/null 2>&1

    longest_seq=$(awk '/^>/ { if (seq && length(seq) > length(longest)) longest=seq; seq=""; next }
                      { seq = seq $0 }
                      END { if (length(seq) > length(longest)) longest=seq; print longest }' "$orf_file")

    compression_output=$(JARVIS3 -t 4 -lr 0 -cm 1:1:0:0.9/0:0:0:0 "$seq_file")
    compression_nc=$(echo "$compression_output" | grep -oP '[0-9.]+(?= NC)' | head -1)

    {
        flock 200
        echo "\"$species\",\"$label\",\"$seq\",\"$longest_seq\",\"$compression_nc\"" >> "$output_result"
    } 200>"$output_result.lock"

    rm -f "$seq_file" "$fasta_file" "$orf_file"
}


generate_dataset(){
    local FASTA_DIR="$1"
    local GFF_DIR="$2"
    local WINDOW_SIZE="$3"
    local STEP="$4"
    local OUTPUT="$5"

    MAX_JOBS=20

    FINAL_DATASET_DIR="${PLANT_DIR}/dataset/final"

    mkdir -p "${FINAL_DATASET_DIR}"

    RESULTS_OUTPUT="${FINAL_DATASET_DIR}/results_${OUTPUT}"

    echo "\"species\",\"label\",\"seq\",\"longest_seq\",\"compression_nc\"" > "$RESULTS_OUTPUT"
    mkdir -p ${FINAL_DATASET_DIR}/tmp

    for FASTA in "$FASTA_DIR"/*; do
        SPECIES=$(basename "${FASTA%.fa}")
        GFF="$GFF_DIR/${SPECIES}.gff3"

        echo "[INFO] Processing $SPECIES"
        [[ -f "${FASTA}.fai" ]] || samtools faidx "$FASTA"

        mkdir -p ${FINAL_DATASET_DIR}/tmp_$SPECIES
        GENES_BED="${FINAL_DATASET_DIR}/tmp_$SPECIES/genes.bed"
        WINDOWS="${FINAL_DATASET_DIR}/tmp_$SPECIES/windows.bed"
        GENIC="${FINAL_DATASET_DIR}/tmp_$SPECIES/genic.bed"
        INTERGENIC="${FINAL_DATASET_DIR}/tmp_$SPECIES/intergenic.bed"
        COMBINED="${FINAL_DATASET_DIR}/tmp_$SPECIES/combined.bed"

        awk 'BEGIN{OFS="\t"} $3 == "gene" && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ {
            strand = ($7 ~ /^[+-]$/ ? $7 : "+");
            print $1, $4-1, $5, ".", ".", strand;
        }' "$GFF" > "$GENES_BED"

        awk -v w=$WINDOW_SIZE -v s=$STEP '{
            chr=$1; len=$2;
            for (i = 0; i + w < len; i += s) {
                printf "%s\t%d\t%d\t.\t.\t+\n", chr, i, i+w;
                printf "%s\t%d\t%d\t.\t.\t-\n", chr, i, i+w;
            }
        }' "${FASTA}.fai" > "$WINDOWS"

        bedtools intersect -s -a "$WINDOWS" -b "$GENES_BED" -f 1.0 -u | \
            awk 'BEGIN{OFS="\t"} {print $0, "gene"}'> "$GENIC"

        bedtools intersect -s -a "$WINDOWS" -b "$GENES_BED" -v | \
            awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' > "$INTERGENIC"

        cat "$GENIC" "$INTERGENIC" > "$COMBINED"

        running_jobs=0
        while IFS=$'\t' read -r chrom start end _ _ strand label; do
            (
                process_window "$FASTA" "$chrom" "$start" "$end" "$strand" "$label" "$SPECIES" "$RESULTS_OUTPUT" "$FINAL_DATASET_DIR"
            ) &
            ((running_jobs++))
            if (( running_jobs >= MAX_JOBS )); then
                wait -n
                ((running_jobs--))
            fi
        done < "$COMBINED"

        wait
        echo "[DONE] $SPECIES"
    done

    rm -rf "${FINAL_DATASET_DIR}/tmp" "${RESULTS_OUTPUT}.lock"

    ${PLANT_DIR}/bin/converter -i $RESULTS_OUTPUT -o ${FINAL_DATASET_DIR}/dataset_${OUTPUT} -t A,C,T,G,ATG,TAA,TAG,TGA,TATAAA,CAAT,CGCG,GT,AG,AAA,TTT,GGG,CGG,CAG

    rm "$RESULTS_OUTPUT"
}

Help() {
   echo
   echo "Generate dataset for genome annotation."
   echo
   echo "Syntax: $0 [-h] -f FASTA_DIR -d ANNOTATION_DIR [-w WINDOW_SIZE] [-s STEP] [-o OUTPUT]"
   echo
   echo "options:"
   echo "-h   Print this help and exit"
   echo "-f   Path to folder containing FASTA files (mandatory)"
   echo "-d   Path to folder containing GFF3 annotation files (mandatory)"
   echo "-w   Window size (default: 1500)"
   echo "-s   Step size (default: 50)"
   echo "-o   Output file name (default: dataset_<WINDOW>_<STEP>.csv)"
   echo
}

if [[ -z "$PLANT_DIR" ]]; then
    echo "Error: PLANT_DIR environment variable is not set."; exit 1
fi

for cmd in samtools getorf JARVIS3 bedtools; do
    if ! command -v $cmd >/dev/null 2>&1; then
        echo "Error: Required tool '$cmd' not found in PATH."; exit 1
    fi
done

WINDOW_SIZE="1500"
STEP="50"

while getopts "hf:d:w:s:o:" option; do
   case $option in
        h) Help; exit;;
        f) FASTA_DIR=$OPTARG;;
        d) ANNOTATION_FOLDER=$OPTARG;;
        w) WINDOW_SIZE=$OPTARG;;
        s) STEP=$OPTARG;;
        o) OUTPUT_FILE=$OPTARG;;
        \?) echo "Error: Invalid option."; echo "Use $0 -h for help."; exit 1;;
   esac
done

if [ $OPTIND -eq 1 ]; then
    echo "Error: No flag provided."
    echo "Use $0 -h for help."
    exit 1
fi

if [[ -z "$FASTA_DIR" || -z "$ANNOTATION_FOLDER" ]]; then
    echo "Error: Both -f (FASTA_DIR) and -d (ANNOTATION_DIR) are required."
    echo "Use $0 -h for help."
    exit 1
fi

if [[ -z "$OUTPUT_FILE" ]]; then
    OUTPUT_FILE="dataset_${WINDOW_SIZE}_${STEP}.csv"
fi

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "Error: FASTA_DIR ($FASTA_DIR) does not exist."; exit 1
fi

if [[ ! -d "$ANNOTATION_FOLDER" ]]; then
    echo "Error: ANNOTATION_DIR ($ANNOTATION_FOLDER) does not exist."; exit 1
fi

generate_dataset "$FASTA_DIR" "$ANNOTATION_FOLDER" "$WINDOW_SIZE" "$STEP" "$OUTPUT_FILE"