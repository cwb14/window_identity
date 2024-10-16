#!/bin/bash

# Dynamically determine the directory of the current script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"

# Default values
WINDOW_SIZE=1000000
SLIDE_SIZE=500000
MUTATION_RATE=1.3e-8
PROTEIN="$SCRIPT_DIR/Amborella_trichopoda.AMTR1.0.pep.all.fa"
THREADS=4
PROCESSES=4
INCLUDE_MEAN_LINE="no"
YMAX=""
NAMES=""
REF=""
QUERY_GENOMES=()

# Parsing command-line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
    -ref)
        REF="$2"
        shift # past argument
        shift # past value
        ;;
    -query)
        shift # past argument
        while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
            QUERY_GENOMES+=("$1")
            shift
        done
        ;;
    -window_size)
        WINDOW_SIZE="$2"
        shift
        shift
        ;;
    -slide_size)
        SLIDE_SIZE="$2"
        shift
        shift
        ;;
    -peptide)
        PROTEIN="$2"
        shift
        shift
        ;;
    -threads)
        THREADS="$2"
        shift
        shift
        ;;
    -processes)
        PROCESSES="$2"
        shift
        shift
        ;;
    -mutation_rate)
        MUTATION_RATE="$2"
        shift
        shift
        ;;
    -names)
        NAMES="$2"
        shift
        shift
        ;;
    -include_mean_line)
        INCLUDE_MEAN_LINE="$2"
        shift
        shift
        ;;
    -ymax)
        YMAX="$2"
        shift
        shift
        ;;
    *)
        echo "Unknown option $1"
        exit 1
        ;;
    esac
done

# Check for required arguments
if [[ -z "$REF" ]]; then
    echo "Error: Reference genome (-ref) is required."
    exit 1
fi

if [[ ${#QUERY_GENOMES[@]} -eq 0 ]]; then
    echo "Error: At least one query genome (-query) is required."
    exit 1
fi

# Collect all genomes (reference and query)
ALL_GENOMES=("$REF" "${QUERY_GENOMES[@]}")

# Extract genome IDs
declare -A GENOME_IDS
for genome in "${ALL_GENOMES[@]}"; do
    filename=$(basename "$genome")
    id="${filename%%.*}"
    GENOME_IDS["$genome"]="$id"
done

# Step 1 - Cleanup the input genome file
modified_fastas_exist=true
for genome in "${ALL_GENOMES[@]}"; do
    mod_fasta="${GENOME_IDS[$genome]}_mod.fa"
    if [[ ! -s "$mod_fasta" ]]; then
        modified_fastas_exist=false
        break
    fi
done

if [[ ! -s "jcvi_list.txt" ]]; then
    modified_fastas_exist=false
fi

if [[ "$modified_fastas_exist" = false ]]; then
    echo "Step 1 - Cleanup the input genome file."
    echo "Running: python $BIN_DIR/fasta_renamer_diploid.py -genomes ${ALL_GENOMES[@]}"
    python "$BIN_DIR/fasta_renamer_diploid.py" -genomes "${ALL_GENOMES[@]}"
else
    echo "Step 1 - Modified fasta files and jcvi_list.txt exist. Skipping."
fi

# Step 2 - Align proteins to each genomic input
miniprot_outputs_exist=true
for genome in "${ALL_GENOMES[@]}"; do
    id="${GENOME_IDS[$genome]}"
    gff_file="${id}.gff"
    if [[ ! -s "$gff_file" ]]; then
        miniprot_outputs_exist=false
        break
    fi
done

if [[ "$miniprot_outputs_exist" = false ]]; then
    echo "Step 2 - Align proteins to each genomic input."
    for genome in "${ALL_GENOMES[@]}"; do
        id="${GENOME_IDS[$genome]}"
        mod_fasta="${id}_mod.fa"
        gff_file="${id}.gff"
        if [[ ! -s "$gff_file" ]]; then
            echo "Running: miniprot -It $THREADS --gff $mod_fasta $PROTEIN -P $id --outn=10 > $gff_file"
            miniprot -It "$THREADS" --gff "$mod_fasta" "$PROTEIN" -P "$id" --outn=10 >"$gff_file"
        else
            echo "GFF file $gff_file exists. Skipping miniprot for $id."
        fi
    done
else
    echo "Step 2 - GFF files exist. Skipping."
fi

# Step 3 - Convert GFF into bed for each genomic input
for genome in "${ALL_GENOMES[@]}"; do
    id="${GENOME_IDS[$genome]}"
    gff_file="${id}.gff"
    bed_file="${id}.bed"
    if [[ ! -s "$bed_file" ]]; then
        echo "Converting $gff_file to $bed_file"
        cat "$gff_file" | grep 'mRNA' | awk -F '\t' '{split($9,a,";"); split(a[1],id,"="); print $1 "\t" $4 "\t" $5 "\t" id[2]}' | grep "$id" | grep -v 'mapped' | awk 'NF==4' | bedtools sort -i - >"$bed_file"
    else
        echo "BED file $bed_file exists. Skipping."
    fi
done

# Step 4 - Generate pseudo-CDS for each genomic input
for genome in "${ALL_GENOMES[@]}"; do
    id="${GENOME_IDS[$genome]}"
    mod_fasta="${id}_mod.fa"
    bed_file="${id}.bed"
    cds_file="${id}.cds"
    if [[ ! -s "$cds_file" ]]; then
        echo "Generating pseudo-CDS for $id"
        bedtools getfasta -fi "$mod_fasta" -bed "$bed_file" -name | bioawk -c fastx '{ sub(/::.*/, "", $name); print ">"$name"\n"$seq }' >"$cds_file"
    else
        echo "CDS file $cds_file exists. Skipping."
    fi
done

# Step 5 - Align pseudo-CDS
expected_anchor_files=()
while read -r line; do
    ID1=$(echo "$line" | awk '{print $1}')
    ID2=$(echo "$line" | awk '{print $2}')
    anchor_file="${ID1}.${ID2}.anchors"
    expected_anchor_files+=("$anchor_file")
done < jcvi_list.txt

anchors_exist=true
for anchor_file in "${expected_anchor_files[@]}"; do
    if [[ ! -s "$anchor_file" ]]; then
        anchors_exist=false
        break
    fi
done

if [[ "$anchors_exist" = false ]]; then
    echo "Step 5 - Align pseudo-CDS"
    echo "Running: python $BIN_DIR/jcvi_diploid.py -p $PROCESSES --cpus $THREADS"
    python "$BIN_DIR/jcvi_diploid.py" -p "$PROCESSES" --cpus "$THREADS"
    rm -f *.nsq *.nin *.nhr *.ndb *.nto *.not *.ntf *.njs *.des *.sds *.tis *.ssp *.bck *.suf *.prj

    # Initialize current processes and threads
    current_p="$PROCESSES"
    current_cpus="$THREADS"

    for i in {1..200}; do
        anchors_exist=true
        for anchor_file in "${expected_anchor_files[@]}"; do
            if [[ ! -s "$anchor_file" ]]; then
                anchors_exist=false
                break
            fi
        done

        if [[ "$anchors_exist" = true ]]; then
            echo "All anchor files exist after $i attempt(s)."
            break
        fi

        echo "Retrying jcvi_diploid_retry.py - attempt $i"

        # Decrement processes and threads, ensuring they don't go below 1
        if [[ $i -ne 1 ]]; then
            current_p=$((current_p / 2))
            current_cpus=$((current_cpus / 2))
            # Ensure minimum value of 1
            if [[ $current_p -lt 1 ]]; then
                current_p=1
            fi
            if [[ $current_cpus -lt 1 ]]; then
                current_cpus=1
            fi
        fi

        echo "Using -p $current_p --cpus $current_cpus"

        python "$BIN_DIR/jcvi_diploid_retry.py" -p "$current_p" --cpus "$current_cpus"
        rm -f *.nsq *.nin *.nhr *.ndb *.nto *.not *.ntf *.njs *.des *.sds *.tis *.ssp *.bck *.suf *.prj
    done

    # Final check after all attempts
    anchors_exist=true
    for anchor_file in "${expected_anchor_files[@]}"; do
        if [[ ! -s "$anchor_file" ]]; then
            anchors_exist=false
            break
        fi
    done

    if [[ "$anchors_exist" = false ]]; then
        echo "Error: Unable to generate all anchor files after multiple attempts."
        exit 1
    fi
else
    echo "Step 5 - Anchor files exist. Skipping."
fi


# Steps 6-8 for each pair in jcvi_list.txt
while read -r line; do
    ID1=$(echo "$line" | awk '{print $1}')
    ID2=$(echo "$line" | awk '{print $2}')
    anchor_file="${ID1}.${ID2}.anchors"
    clean_anchor_file="${ID1}.${ID2}.clean.anchors"

    # Step 6 - Clean the anchor files
    if [[ ! -s "$clean_anchor_file" ]]; then
        echo "Cleaning $anchor_file to $clean_anchor_file"
        python "$BIN_DIR/anchor_builder.py" "$anchor_file" >"$clean_anchor_file"
    else
        echo "Cleaned anchor file $clean_anchor_file exists. Skipping."
    fi

    # Step 7 - Convert the cleaned anchor files to coordinate format
    raw_coords_file="${ID1}.${ID2}.anchors.raw.coords"
    if [[ ! -s "$raw_coords_file" ]]; then
        echo "Converting $clean_anchor_file to $raw_coords_file"
        python "$BIN_DIR/gene_coords_extractor_all4.py" -mcscan "$clean_anchor_file" | sort | uniq >"$raw_coords_file"
    else
        echo "Raw coords file $raw_coords_file exists. Skipping."
    fi

    # Step 8 - Consolidate overlapping anchors
    coords_file="${ID1}.${ID2}.anchors.coords"
    if [[ ! -s "$coords_file" ]]; then
        echo "Consolidating $raw_coords_file to $coords_file"
        python "$BIN_DIR/anchor_coord_consolidator.py" "$raw_coords_file" >"$coords_file"
    else
        echo "Coords file $coords_file exists. Skipping."
    fi

done <jcvi_list.txt

# Step 9 - Merge the pairwise anchor files
coords_files=()
while read -r line; do
    ID1=$(echo "$line" | awk '{print $1}')
    ID2=$(echo "$line" | awk '{print $2}')
    coords_file="${ID1}.${ID2}.anchors.coords"
    coords_files+=("$coords_file")
done <jcvi_list.txt

echo "Step 9 - Merge the pairwise anchor files"
cat "${coords_files[@]}" >all.anchors.coords

# Step 10 - Align the genomic anchors
if [[ ! -s "alignment.paf" ]]; then
    echo "Step 10 - Align the genomic anchors"
    python "$BIN_DIR/synmap.py" --timer 8h -t "$THREADS" -p "$PROCESSES" -c all.anchors.coords >/dev/null 2>&1 # I could add timer as a command line parameter.
    # I can use 'paftools.js view' to visualize aln quality and optimize parameters.
    # The minimap2 parameters likely need optimized since k2p seems off.
else
    echo "Step 10 - alignment.paf exists. Skipping."
fi

# Step 11 - Adjust the PAF coordinates, calculate the weighted divergence, and build the tree
if [[ ! -s "alignment_adjust.paf" ]]; then
    echo "Adjusting alignment.paf to alignment_adjust.paf"
    python "$BIN_DIR/synmap_adjust.py" -paf alignment.paf | \
    awk '{for(i=1;i<=NF;i++) if($i ~ /^cg:Z:/){c=substr($i,6);s=0; while(match(c,/([0-9]+)M/,a)){s+=a[1];c=substr(c,RSTART+RLENGTH)}; print $0 "\t" s}}' | \
    python "$BIN_DIR/calculate_k2p.py" | sed 's/-0.000000/0.000000/g' > alignment_adjust.paf
    # 'calculate_k2p.py' includes '-d' and '-r' flags that may help filter low confidence variants and improve k2p.
else
    echo "alignment_adjust.paf exists. Skipping."
fi

if [[ ! -s "alignment_adjust.tsv" ]]; then
    echo "Calculating weighted divergence"
    python "$BIN_DIR/weighted_paf.py" -paf alignment_adjust.paf
else
    echo "alignment_adjust.tsv exists. Skipping."
fi

if [[ ! -s "alignment_genome.tsv" ]]; then
    echo "Calculating weighted averages"
    python "$BIN_DIR/weighted_average3.py" -input alignment_de.tsv -chrom_out de_chrom.tsv -genome_out de_genome.tsv
    python "$BIN_DIR/weighted_average3.py" -input alignment_k2p.tsv -chrom_out k2p_chrom.tsv -genome_out k2p_genome.tsv
else
    echo "alignment_genome.tsv exists. Skipping."
fi

if [[ ! -s "ANI_matrix.tsv" ]]; then
    echo "Building the matrix"
    python "$BIN_DIR/matrix_builder.py" -in de_genome.tsv > de_matrix.tsv
    python "$BIN_DIR/matrix_builder.py" -in k2p_genome.tsv > k2p_matrix.tsv
else
    echo "ANI_matrix.tsv exists. Skipping."
fi

# Step 12 - Subset the alignment file according to the reference
REF_ID="${GENOME_IDS[$REF]}"
REF_ID="${REF_ID%%.*}"

subset_pafs_exist=true
for genome in "${QUERY_GENOMES[@]}"; do
    QUERY_ID="${GENOME_IDS[$genome]}"
    bed_file="${REF_ID}.${QUERY_ID}.bed"
    if [[ ! -s "$bed_file" ]]; then
        subset_pafs_exist=false
        break
    fi
done

if [[ "$subset_pafs_exist" = false ]]; then
    echo "Step 12 - Subsetting alignment file according to the reference"
    python "$BIN_DIR/paf_to_bed.py" -ref "$REF_ID" -paf alignment_adjust.paf
else
    echo "Subset PAF files exist. Skipping."
fi

# Step 13 - Build the window file
if [[ ! -s "${REF_ID}_mod.fa.fai" ]]; then
    echo "Indexing ${REF_ID}_mod.fa"
    samtools faidx "${REF_ID}_mod.fa"
fi

window_file="${REF_ID}_mod.window"
if [[ ! -s "$window_file" ]]; then
    echo "Building window file for $REF_ID"
    bedtools makewindows -g <(cut -f 1,2 "${REF_ID}_mod.fa.fai") -w "$WINDOW_SIZE" -s "$SLIDE_SIZE" >"$window_file"
else
    echo "Window file $window_file exists. Skipping."
fi

# Step 14 - Calculate weights for each window in all 'REF.QUERY.paf' files
for genome in "${QUERY_GENOMES[@]}"; do
    QUERY_ID="${GENOME_IDS[$genome]}"
    bed_file="${REF_ID}.${QUERY_ID}.bed"
    output_tsv="${REF_ID}.${QUERY_ID}.tsv"
    if [[ ! -s "$output_tsv" ]]; then
        echo "Calculating weights for $bed_file"
        python "$BIN_DIR/weighted_de_scores.py" -window_bed "${REF_ID}_mod.window" -minimap_bed "$bed_file" -output "$output_tsv" --threads "$THREADS"
    else
        echo "Output TSV $output_tsv exists. Skipping."
    fi
done

# Step 15 - Calculate the weighted average for each window
for genome in "${QUERY_GENOMES[@]}"; do
    QUERY_ID="${GENOME_IDS[$genome]}"
    input_tsv="${REF_ID}.${QUERY_ID}.tsv"
    avg_tsv="${REF_ID}.${QUERY_ID}.avg.tsv"
    if [[ ! -s "$avg_tsv" ]]; then
        echo "Calculating weighted average for $input_tsv"
        python "$BIN_DIR/weighted_de_average.py" "$input_tsv" "${QUERY_ID}" >"$avg_tsv"
    else
        echo "Average TSV $avg_tsv exists. Skipping."
    fi
done

# Step 16 - Generate per-chromosome divergence report
for genome in "${QUERY_GENOMES[@]}"; do
    QUERY_ID="${GENOME_IDS[$genome]}"
    bed_file="${REF_ID}.${QUERY_ID}.bed"
    chrom_report="${REF_ID}.${QUERY_ID}.chrom.tsv" # Assuming the script generates output with this name
    if [[ ! -s "$chrom_report" ]]; then
        echo "Generating per-chromosome divergence report for $bed_file"
        python "$BIN_DIR/weighted_de_chroms.py" -bed "$bed_file" -fai "${REF_ID}_mod.fa.fai" > "$chrom_report"
    else
        echo "Chromosome divergence report $chrom_report exists. Skipping."
    fi
done

# Step 17 - Merge the window weighted averages for plotting
avg_tsv_files=()
for genome in "${QUERY_GENOMES[@]}"; do
    QUERY_ID="${GENOME_IDS[$genome]}"
    avg_tsv="${REF_ID}.${QUERY_ID}.avg.tsv"
    avg_tsv_files+=("$avg_tsv")
done

merged_file="merged.weighted_de_avg.bed"
if [[ ! -s "$merged_file" ]]; then
    echo "Merging weighted averages for plotting"
    python "$BIN_DIR/window_merger.py" "${avg_tsv_files[@]}" -o "$merged_file"
else
    echo "Merged weighted averages file $merged_file exists. Skipping."
fi

# Step 18 - Plot divergence across chromosomes
echo "Plotting divergence across chromosomes"
plot_cmd="Rscript $BIN_DIR/divergence_plotter.R -in merged.weighted_de_avg.bed -include_mean_line \"$INCLUDE_MEAN_LINE\""
if [[ -n "$YMAX" ]]; then
    plot_cmd+=" -ymax $YMAX"
fi
echo "Running: $plot_cmd"
eval "$plot_cmd"

# Build the tree from the cleaned matrix.
echo "Step 19 - Building phylogenetic tree using UPGMA"
upgma_cmd="Rscript $BIN_DIR/upgma.R --method upgma"
if [[ -n "$MUTATION_RATE" ]]; then
    upgma_cmd+=" --mutation_rate $MUTATION_RATE"
fi
if [[ -n "$NAMES" ]]; then
    upgma_cmd+=" --name_key $NAMES"
fi
echo "Running: $upgma_cmd"
eval "$upgma_cmd"
