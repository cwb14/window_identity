#!/bin/bash

# Dynamically determine the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# Set the path to the 'bin' directory relative to the script's location
BIN_DIR="$SCRIPT_DIR/bin"

# Default values
threads=3
window_size=1000000
slide_size=500000
include_mean_line=false
clean=false

# Function to print usage/help
usage() {
    cat << EOF
Usage: $0 -ref ref.fa -query asm1.fa asm2.fa asm3.fa [OPTIONS]

Plots gap compressed sequence identity across chromosomes using weighted averages to calculate sequence identity per window.
Supports 1v1, 1v2, 1v3, etc.

Arguments:
  -ref              Reference genome file (required)
  -query            One or more query genome files (required)
  -threads [INT]    Number of threads for minimap2 (optional, default: 3)
  -window_size [INT] Window size for bedtools (optional, default: 1000000)
  -slide_size [INT] Slide size for bedtools (optional, default: 500000)
  include_mean_line Boolean flag to include mean line in the R plot (optional)
  -clean            Boolean flag to remove intermediate files after completion (optional)

Examples:
  bash $0 -ref ref.fa -query query1.fa query2.fa -threads 40 -window_size 1000000 -slide_size 500000 include_mean_line -clean

EOF
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -ref)
        ref="$2"
        shift # past argument
        shift # past value
        ;;
        -query)
        shift # past argument
        while [[ $# -gt 0 && "$1" != -* ]]; do
            queries+=("$1")
            shift # past value
        done
        ;;
        -threads)
        threads="$2"
        shift # past argument
        shift # past value
        ;;
        -window_size)
        window_size="$2"
        shift # past argument
        shift # past value
        ;;
        -slide_size)
        slide_size="$2"
        shift # past argument
        shift # past value
        ;;
        include_mean_line)
        include_mean_line=true
        shift # past argument
        ;;
        -clean)
        clean=true
        shift # past argument
        ;;
        *)
        usage
        ;;
    esac
done

# Check if reference and query genomes are provided
if [ -z "$ref" ] || [ ${#queries[@]} -eq 0 ]; then
    usage
fi

ref_basename=$(basename "${ref%.fa}")
window_bed="${ref_basename}_${window_size}bp_window.bed"
fai_file="${ref_basename}.fa.fai"

# Create output files in the current directory
output_pafs=()
output_beds=()
weighted_de_tsvs=()
weighted_de_avg_beds=()
reports=()

for query in "${queries[@]}"; do
    query_basename=$(basename "${query%.fa}")
    output_pafs+=("${ref_basename}.${query_basename}.paf")
    output_beds+=("${ref_basename}.${query_basename}.paf.bed")
    weighted_de_tsvs+=("${ref_basename}.${query_basename}.weighted_de.tsv")
    weighted_de_avg_beds+=("${ref_basename}.${query_basename}.weighted_de_avg.bed")
    reports+=("${ref_basename}.${query_basename}.report")
done

echo "Step 1: Align queries to reference"
# Align queries to reference using minimap2
for i in "${!queries[@]}"; do
    if [ -f "${output_pafs[$i]}" ]; then
        echo "PAF file ${output_pafs[$i]} already exists, skipping minimap2..."
    else
        echo "Running minimap2 for ${queries[$i]}..."
        minimap2 -c -t "$threads" --secondary=no -x asm5 "$ref" "${queries[$i]}" > "${output_pafs[$i]}" 2>/dev/null
    fi
done

echo "Step 2: Create windows file for the reference"
# Create a window file for the reference genome
if [ -f "$window_bed" ]; then
    echo "Window BED file $window_bed already exists, skipping bedtools makewindows..."
else
    samtools faidx "$ref"
    bedtools makewindows -g <(cut -f 1,2 "$fai_file") -w "$window_size" -s "$slide_size" > "$window_bed"
fi

echo "Step 3: Convert PAF to bed"
# Convert PAF to BED format
for i in "${!queries[@]}"; do
    if [ -f "${output_beds[$i]}" ]; then
        echo "BED file ${output_beds[$i]} already exists, skipping PAF to BED conversion..."
    else
        cat "${output_pafs[$i]}" | cut -f 6,8,9,21 > "${output_beds[$i]}"
    fi
done

echo "Step 4: Calculate weighted de scores"
# Calculate weighted DE scores
for i in "${!queries[@]}"; do
    if [ -f "${weighted_de_tsvs[$i]}" ]; then
        echo "Weighted DE TSV file ${weighted_de_tsvs[$i]} already exists, skipping weighted_de_scores.py..."
    else
        python "$BIN_DIR/weighted_de_scores.py" -window_bed "$window_bed" -minimap_bed "${output_beds[$i]}" -output "${weighted_de_tsvs[$i]}" &
    fi
done
wait

echo "Step 5: Generate de scores per window"
# Generate DE scores per window
for i in "${!queries[@]}"; do
    if [ -f "${weighted_de_avg_beds[$i]}" ]; then
        echo "Weighted DE avg BED file ${weighted_de_avg_beds[$i]} already exists, skipping weighted_de_average.py..."
    else
        python "$BIN_DIR/weighted_de_average.py" "${weighted_de_tsvs[$i]}" "${queries[$i]%.fa}" > "${weighted_de_avg_beds[$i]}"
    fi
done

echo "Step 6: Calculate weighted average for each chromosome and whole-genome"
# Calculate weighted average for each chromosome and whole-genome
for i in "${!queries[@]}"; do
    if [ -f "${reports[$i]}" ]; then
        echo "Report file ${reports[$i]} already exists, skipping weighted_de_chroms.py..."
    else
        python "$BIN_DIR/weighted_de_chroms.py" -bed "${output_beds[$i]}" -fai "$fai_file" > "${reports[$i]}"
    fi
done

echo "Step 7: Merge alignments by overlap"
# Merge alignments by overlap
python "$BIN_DIR/window_merger.py" "${weighted_de_avg_beds[@]}" -o "merged.weighted_de_avg.bed"

echo "Step 8: Build the R plot"
# Build the R plot
Rscript "$BIN_DIR/divergence_plotter.R" "merged.weighted_de_avg.bed" $([ "$include_mean_line" = true ] && echo "include_mean_line")

# Clean up intermediate files if -clean flag is provided
if [ "$clean" = true ]; then
    echo "Cleaning up intermediate files..."
    rm "$fai_file"
    rm "$window_bed"
    for i in "${!queries[@]}"; do
        rm "${output_pafs[$i]}"
        rm "${output_beds[$i]}"
        rm "${weighted_de_tsvs[$i]}"
        rm "${weighted_de_avg_beds[$i]}"
    done
fi

echo "Pipeline completed successfully."
