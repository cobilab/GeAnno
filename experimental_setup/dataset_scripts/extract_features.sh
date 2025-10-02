 #!/bin/bash

revcomp() {
    echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev
}

get_seeded_random()
{
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

DATASET_DIR="${PLANT_DIR}/experimental_setup/datasets/"

for dataset_dir in "${DATASET_DIR}/"*; do
    dataset=$(basename "$dataset_dir")

    FASTA_DIR="${dataset_dir}/dna/"
    GFF_DIR="${dataset_dir}/annotation/"
    WINDOW_SIZE=1500
    STRIDE=50

    OUTPUT="${dataset}_${WINDOW_SIZE}_${STRIDE}.csv"
    MAX_JOBS=20

    FINAL_DATASET_DIR="${dataset_dir}/final"

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

        awk -v w=$WINDOW_SIZE -v s=$STRIDE '{
            chr=$1; len=$2;
            for (i = 0; i + w < len; i += s) {
                printf "%s\t%d\t%d\t.\t.\t+\n", chr, i, i+w;
                printf "%s\t%d\t%d\t.\t.\t-\n", chr, i, i+w;
            }
        }' "${FASTA}.fai" > "$WINDOWS"

        bedtools intersect -s -a "$WINDOWS" -b "$GENES_BED" -f 1.0 -u | \
            awk 'BEGIN{OFS="\t"} {print $0, "gene"}'> "$GENIC"

        if [[ "$dataset" == "m_esculenta" ]]; then
            GENIC_COUNT=$(wc -l < "$GENIC")
            COUNT_INTERGENIC=$(( 2 * GENIC_COUNT ))

            bedtools intersect -s -a "$WINDOWS" -b "$GENES_BED" -v | \
                awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' | \
                shuf --random-source=<(get_seeded_random 42) -n $COUNT_INTERGENIC> "$INTERGENIC"

        else
            bedtools intersect -s -a "$WINDOWS" -b "$GENES_BED" -v | \
                awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' > "$INTERGENIC"
        fi

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

    if [[ "$dataset" == "m_esculenta" ]]; then
        ${PLANT_DIR}/bin/converter -i $RESULTS_OUTPUT -o ${FINAL_DATASET_DIR}/temporary_set_${OUTPUT} -t A,C,T,G,ATG,TAA,TAG,TGA,TATAAA,CAAT,CGCG,GT,AG,AAA,TTT,GGG,CGG,CAG

        p=0.40
        gawk -v p="$p" -v seed=42 'BEGIN{srand(seed)} { if (rand() < p) print }' ${FINAL_DATASET_DIR}/temporary_set_${OUTPUT} > ${FINAL_DATASET_DIR}/dataset_${OUTPUT}

        rm "${FINAL_DATASET_DIR}/temporary_set_${OUTPUT}"
    else
        ${PLANT_DIR}/bin/converter -i $RESULTS_OUTPUT -o ${FINAL_DATASET_DIR}/dataset_${OUTPUT} -t A,C,T,G,ATG,TAA,TAG,TGA,TATAAA,CAAT,CGCG,GT,AG,AAA,TTT,GGG,CGG,CAG
    fi

    rm "$RESULTS_OUTPUT"
done

