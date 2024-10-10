#!/bin/bash

# K2P requires we know the number of the total positions compared. 
# An alignment is two sided (A->B & B-A).
# You'd think the total number of compared positions would be identical regardless of whether you looked at it from A->B or B->A, but seems not so (maybe due to differences in gaps).
# This script looks at A->B & B->A and picks the smaller of the two to use as the total number of positions compared during the k2p calculation.
# It relies on 'paftools.js bedcov', which reports 'target bases overlapping regions'. 'Overlapping region' is the genome bed. 'target bases' are ones that were mapped during the full genome alignment.
# The 'paftools.js bedcov' 'target bases overlapping regions' value should be similar to the total positions compared.
# 'paftools.js stat' reports a similar value with 'Number of mapped bases', but this value appears to be REF dependent. Ie, its either A->B or B->A but not both.
# Im not sure that this is the best way to calculate total positions compared, but seems like a good approximate.
# If all aligned sequences had coverage=1, I could simply subtract the sum of indels from the sum of matches and missmatches. Coverage is != 1. 

# Usage: bash script.sh [--keep_intermediates]

# Handle '--keep_intermediates' flag
if [[ "$1" == "--keep_intermediates" ]]; then
    KEEP_INTERMEDIATES=true
else
    KEEP_INTERMEDIATES=false
fi

# Dynamically determine the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
BIN_DIR="$SCRIPT_DIR/"

# Initialize arrays to keep track of created files
created_files=()
aln_size_files=()

# Step 1: Create genome bed files and genome size files
mapfile -t genomes < genome_list.txt

echo "Step 1: Creating genome bed files..."
for genome in "${genomes[@]}"; do
    echo "Processing genome: $genome"
    bioawk -c fastx '{ print $name, length($seq) }' < "${genome}_mod.fa" \
        | awk '{print $1, 1, $2}' OFS="\t" > "${genome}_mod.bed"
    created_files+=("${genome}_mod.bed")
done

# Step 2: Partitioned PAF files and convert to BED
echo "Step 2: Partitioning PAF files and converting to BED..."
for genome in "${genomes[@]}"; do
    echo "Running paf_to_bed.py for genome: $genome"
    python "$BIN_DIR/paf_to_bed.py" -ref "$genome" -paf alignment_adjust.paf

    # Collect generated bed files
    for bed_file in ${genome}.*.bed; do
        if [[ "$bed_file" != "${genome}_mod.bed" ]]; then
            created_files+=("$bed_file")
        fi
    done
done

# Step 3: Run paftools for all reciprocal pairs
echo "Step 3: Running paftools.js bedcov for all reciprocal pairs..."
for genome in "${genomes[@]}"; do
    for bed_file in ${genome}.*.bed; do
        if [[ "$bed_file" != "${genome}_mod.bed" ]]; then
            output_file="${bed_file}.aln.size"
            echo "Calculating alignment size for $bed_file"
            paftools.js bedcov "${genome}_mod.bed" "$bed_file" 2>&1 \
                | awk '/target bases overlapping regions:/ {print $(NF-1)}' > "$output_file"
            created_files+=("$output_file")
            aln_size_files+=("$output_file")
        fi
    done
done

# Step 4: Identify reciprocal .aln.size files and remove the larger one
echo "Step 4: Identifying reciprocal .aln.size file pairs and removing the larger one..."

# Create an associative array to track processed pairs
declare -A processed_pairs

# Iterate over all .aln.size files
for aln_file in *.aln.size; do
    # Extract genome IDs from filename (format: GenomeA.GenomeB.bed.aln.size)
    if [[ "$aln_file" =~ ^([^\.]+)\.([^\.]+)\.bed\.aln\.size$ ]]; then
        genome1="${BASH_REMATCH[1]}"
        genome2="${BASH_REMATCH[2]}"
        
        # Create a sorted key to identify reciprocal pairs uniquely
        if [[ "$genome1" < "$genome2" ]]; then
            pair_key="${genome1}.${genome2}"
            reciprocal_file="${genome2}.${genome1}.bed.aln.size"
        else
            pair_key="${genome2}.${genome1}"
            reciprocal_file="${genome1}.${genome2}.bed.aln.size"
        fi
        
        # Check if this pair has already been processed
        if [[ -z "${processed_pairs["$pair_key"]}" ]]; then
            # Mark this pair as processed
            processed_pairs["$pair_key"]=1
            
            # Check if reciprocal file exists
            if [[ -f "$reciprocal_file" ]]; then
                # Read the numerical values from both files
                num1=$(<"$aln_file")
                num1=$(echo "$num1" | tr -d '[:space:]')
                num2=$(<"$reciprocal_file")
                num2=$(echo "$num2" | tr -d '[:space:]')
                
                # Validate that both num1 and num2 are integers
                if ! [[ "$num1" =~ ^[0-9]+$ && "$num2" =~ ^[0-9]+$ ]]; then
                    echo "Warning: Non-integer value found in $aln_file or $reciprocal_file. Skipping this pair."
                    continue
                fi
                
                echo "Comparing $aln_file ($num1) and $reciprocal_file ($num2)"
                
                # Compare the two sizes and remove the larger file
                if (( num1 > num2 )); then
                    echo "Removing larger file: $aln_file"
                    rm -f "$aln_file"
                    # Remove from tracking arrays
                    aln_size_files=("${aln_size_files[@]/$aln_file}")
                    created_files=("${created_files[@]/$aln_file}")
                elif (( num1 < num2 )); then
                    echo "Removing larger file: $reciprocal_file"
                    rm -f "$reciprocal_file"
                    # Remove from tracking arrays
                    aln_size_files=("${aln_size_files[@]/$reciprocal_file}")
                    created_files=("${created_files[@]/$reciprocal_file}")
                else
                    echo "Both files have equal sizes. Keeping both: $aln_file and $reciprocal_file"
                fi
            else
                echo "No reciprocal file found for $aln_file. Skipping..."
            fi
        fi
    else
        echo "Warning: Filename $aln_file does not match the expected pattern. Skipping..."
    fi
done


# Step 5: Cleanup
if [[ "$KEEP_INTERMEDIATES" == false ]]; then
    echo "Step 5: Cleaning up intermediate files..."
    for file in "${created_files[@]}"; do
        if [[ ! " ${aln_size_files[@]} " =~ " ${file} " ]]; then
            rm -f "$file"
        fi
    done
    echo "Cleanup completed."
else
    echo "Intermediate files are kept as per '--keep_intermediates' flag."
fi

echo "Script execution completed."
