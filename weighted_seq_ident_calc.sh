# Chris Benson
# 07/24/2024

#!/bin/bash

# Dynamically determine the directory of the current script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# Set the path to the 'bin' directory relative to the script's location.
BIN_DIR="$SCRIPT_DIR/bin"

# Default values.
threads=3
window_size=1000000
slide_size=500000
include_mean_line="no"
clean=false
use_snps=false
ymax=""

# Function to print usage/help.
usage() {
    cat << EOF
Usage: $0 -ref ref.fa -query asm1.fa asm2.fa asm3.fa [OPTIONS]

Plots gap compressed sequence identity (or SNP frequency (-use_snps)) across chromosomes using weighted averages to calculate sequence identity per window.
Gap Compressed sequence Identity is a measure of Average Nucleotide Identity (ANI).
Supports 1v1, 1v2, 1v3, etc.

Arguments:
  -ref              Reference genome file (required)
  -query            One or more query genome files (required)
  -threads [INT]    Number of threads (optional, default: 3)
  -window_size [INT] Window size (optional, default: 1000000)
  -slide_size [INT] Slide size (optional, default: 500000)
  -include_mean_line <yes|no> Include an average density line in the R plot (optional, default: no)
  -clean            Remove intermediate files after completion (optional)
  -use_snps         Use SNP frequency instead of gap compressed sequence identity (optional, default: off, 2x slower)
  -ymax [FLOAT]     Maximum y-axis value for the R plot (optional, default: use the max observed y-value)

Examples:
  bash $0 -ref ref.fa -query query1.fa query2.fa -threads 40 -window_size 1000000 -slide_size 500000 -include_mean_line yes -clean -use_snps -ymax 0.01

EOF
    exit 1
}

# Parse arguments.
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -ref)
        ref="$2"
        shift # past argument.
        shift # past value.
        ;;
        -query)
        shift # past argument.
        while [[ $# -gt 0 && "$1" != -* ]]; do
            queries+=("$1")
            shift # past value.
        done
        ;;
        -threads)
        threads="$2"
        shift # past argument.
        shift # past value.
        ;;
        -window_size)
        window_size="$2"
        shift # past argument.
        shift # past value.
        ;;
        -slide_size)
        slide_size="$2"
        shift # past argument.
        shift # past value.
        ;;
        -include_mean_line)
        include_mean_line="$2"
        shift # past argument.
        shift # past value.
        ;;
        -clean)
        clean=true
        shift # past argument.
        ;;
        -use_snps)
        use_snps=true
        shift # past argument.
        ;;
        -ymax)
        ymax="$2"
        shift # past argument.
        shift # past value.
        ;;
        *)
        usage
        ;;
    esac
done

# Check if reference and query genomes are provided.
if [ -z "$ref" ] || [ ${#queries[@]} -eq 0 ]; then
    usage
fi

ref_basename=$(basename "${ref%.fa}")
window_bed="${ref_basename}_${window_size}bp_window.bed"
fai_file="${ref_basename}.fa.fai"

# Create output files in the current directory.
output_pafs=()
output_beds=()
weighted_de_tsvs=()
weighted_de_avg_beds=()
reports=()
vcf_files=()
snp_avg_beds=()
snp_avg_reports=()

for query in "${queries[@]}"; do
    query_basename=$(basename "${query%.fa}")
    output_pafs+=("${ref_basename}.${query_basename}.paf")
    output_beds+=("${ref_basename}.${query_basename}.paf.bed")
    weighted_de_tsvs+=("${ref_basename}.${query_basename}.weighted_de.tsv")
    weighted_de_avg_beds+=("${ref_basename}.${query_basename}.weighted_de_avg.bed")
    reports+=("${ref_basename}.${query_basename}.report")
    vcf_files+=("${ref_basename}.${query_basename}.vcf")
    snp_avg_beds+=("${ref_basename}.${query_basename}.snp_avg.bed")
    snp_avg_reports+=("${ref_basename}.${query_basename}.snp_avg.report")
done

if [ "$use_snps" = true ]; then
    echo "Step 1: Align queries to reference and call SNPs"
    for i in "${!queries[@]}"; do
        if [ -f "${vcf_files[$i]}" ]; then
            echo "VCF file ${vcf_files[$i]} already exists, skipping SNP calling..."
        else
            echo "Running minimap2 and SNP calling for ${queries[$i]}..."
            minimap2 --secondary=no -ax asm5 -t "$threads" "$ref" "${queries[$i]}" 2>/dev/null | samtools view -@ "$threads" -bS 2>/dev/null | samtools sort -@ "$threads" -o "${ref_basename}.${query_basename}.sorted.bam" 2>/dev/null && \
            bcftools mpileup -q "$threads" --ff UNMAP,SECONDARY,QCFAIL,DUP --threads "$threads" -Ou -f "$ref" "${ref_basename}.${query_basename}.sorted.bam" 2>/dev/null | \
            bcftools call --threads "$threads" -mv -Ov -o "${vcf_files[$i]}" 2>/dev/null
        fi
    done
else
    echo "Step 1: Align queries to reference"
    for i in "${!queries[@]}"; do
        if [ -f "${output_pafs[$i]}" ]; then
            echo "PAF file ${output_pafs[$i]} already exists, skipping minimap2..."
        else
            echo "Running minimap2 for ${queries[$i]}..."
            minimap2 -c -t "$threads" --secondary=no -x asm5 "$ref" "${queries[$i]}" > "${output_pafs[$i]}" 2>/dev/null
        fi
    done
fi

echo "Step 2: Create windows file for the reference"
if [ -f "$window_bed" ]; then
    echo "Window BED file $window_bed already exists, skipping bedtools makewindows..."
else
    samtools faidx "$ref" 2>/dev/null
    bedtools makewindows -g <(cut -f 1,2 "$fai_file") -w "$window_size" -s "$slide_size" > "$window_bed" 2>/dev/null
fi

if [ "$use_snps" = true ]; then
    echo "Step 3-6: Calculating SNP fraction per window, chromosome, and whole-genome"
    for i in "${!queries[@]}"; do
        if [ -f "${snp_avg_beds[$i]}" ]; then
            echo "SNP avg BED file ${snp_avg_beds[$i]} already exists, skipping SNP fraction calculation..."
        else
            python "$BIN_DIR/window_vcf.py" -window_bed "$window_bed" -vcf "${vcf_files[$i]}" -id "${queries[$i]%.fa}" -output_window "${snp_avg_beds[$i]}" -output_report "${snp_avg_reports[$i]}"
        fi
    done
else
    echo "Step 3: Convert PAF to bed"
    for i in "${!queries[@]}"; do
        if [ -f "${output_beds[$i]}" ]; then
            echo "BED file ${output_beds[$i]} already exists, skipping PAF to BED conversion..."
        else
            cat "${output_pafs[$i]}" | cut -f 6,8,9,21 > "${output_beds[$i]}"
        fi
    done

    echo "Step 4: Calculate weighted de scores"
    for i in "${!queries[@]}"; do
        if [ -f "${weighted_de_tsvs[$i]}" ]; then
            echo "Weighted DE TSV file ${weighted_de_tsvs[$i]} already exists, skipping weighted_de_scores.py..."
        else
            python "$BIN_DIR/weighted_de_scores.py" -window_bed "$window_bed" -minimap_bed "${output_beds[$i]}" -output "${weighted_de_tsvs[$i]}" --threads "$threads"
        fi
    done

    echo "Step 5: Generate de scores per window"
    for i in "${!queries[@]}"; do
        if [ -f "${weighted_de_avg_beds[$i]}" ]; then
            echo "Weighted DE avg BED file ${weighted_de_avg_beds[$i]} already exists, skipping weighted_de_average.py..."
        else
            python "$BIN_DIR/weighted_de_average.py" "${weighted_de_tsvs[$i]}" "${queries[$i]%.fa}" > "${weighted_de_avg_beds[$i]}"
        fi
    done

    echo "Step 6: Calculate weighted average for each chromosome and whole-genome"
    for i in "${!queries[@]}"; do
        if [ -f "${reports[$i]}" ]; then
            echo "Report file ${reports[$i]} already exists, skipping weighted_de_chroms.py..."
        else
            python "$BIN_DIR/weighted_de_chroms.py" -bed "${output_beds[$i]}" -fai "$fai_file" > "${reports[$i]}"
        fi
    done
fi

echo "Step 7: Merge alignments by overlap"
if [ "$use_snps" = true ]; then
    python "$BIN_DIR/window_merger.py" "${snp_avg_beds[@]}" -o "merged.snp_avg.bed"
else
    python "$BIN_DIR/window_merger.py" "${weighted_de_avg_beds[@]}" -o "merged.weighted_de_avg.bed"
fi

echo "Step 8: Build the R plot"
if [ "$use_snps" = true ]; then
    Rscript "$BIN_DIR/divergence_plotter.R" -in "merged.snp_avg.bed" -include_mean_line "$include_mean_line" $([ -n "$ymax" ] && echo "-ymax $ymax")
else
    Rscript "$BIN_DIR/divergence_plotter.R" -in "merged.weighted_de_avg.bed" -include_mean_line "$include_mean_line" $([ -n "$ymax" ] && echo "-ymax $ymax")
fi

if [ "$clean" = true ]; then
    echo "Cleaning up intermediate files..."
    rm "$fai_file"
    rm "$window_bed"
    for i in "${!queries[@]}"; do
        if [ "$use_snps" = true ]; then
#            rm -f "${ref_basename}.${query_basename}.sorted.bam"
#            rm -f "${vcf_files[$i]}"
            rm -f "${snp_avg_beds[$i]}"
#            rm -f "${snp_avg_reports[$i]}"
        else
#            rm -f "${output_pafs[$i]}"
            rm -f "${output_beds[$i]}"
            rm -f "${weighted_de_tsvs[$i]}"
            rm -f "${weighted_de_avg_beds[$i]}"
        fi
    done
fi

echo "Pipeline completed successfully."

# [END]
