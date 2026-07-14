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
# Divergence mapping option for synmap.py
X_TYPE="asm10"
# Protein liftover (Steps 2-4) and anchoring (Step 5).
OUTN=10
OUTS=0.95
OUTC=0.1
TESORTER="yes"
CSCORE=0.99
# Synteny block consolidation (see Steps 6-9). Defaults match synLTR/module1.py.
MIN_BLOCK_SIZE=15000
STITCH_GAPS="yes"

# Help menu function
usage() {
    cat << EOF
Usage: $(basename "$0") -ref REF_GENOME -query GENOME1 [GENOME2 ...] [options]

FASTA inputs (-ref, -query, -peptide) may be plain or compressed with gzip/bgzip
(.gz/.bgz/.bgzf), bzip2 (.bz2), or xz (.xz/.lzma). The genome ID is the filename
minus those suffixes and the FASTA suffix, so 'Psupina.fa.gz' -> 'Psupina'.

Required:
  -ref GENOME                  Reference genome
  -query GENOME [GENOME ...]   Query genome(s)

Options:
  -window_size SIZE            Window size (default: $WINDOW_SIZE)
  -slide_size SIZE             Slide size (default: $SLIDE_SIZE)
  -peptide FILE                Protein file (default: $PROTEIN)
  -threads N                   Number of threads (default: $THREADS)
  -processes N                 Number of processes (default: $PROCESSES)
  -mutation_rate RATE          Mutation rate (default: $MUTATION_RATE)
  -names STRING                Names key file
  -include_mean_line yes|no    Include mean line in plot (default: $INCLUDE_MEAN_LINE)
  -ymax YMAX                   Y-axis maximum for plot
  -x asm5|asm10|asm20          Sequence divergence mapping for synmap.py (default: $X_TYPE)
  -outn N                      miniprot --outn: max alignments reported per protein
                               (default: $OUTN). Raise for polyploids.
  -outs FLOAT                  miniprot --outs: keep alignments scoring >= this fraction of
                               the best hit for that protein (default: $OUTS).
  -outc FLOAT                  miniprot --outc: min fraction of the protein that must align
                               (default: $OUTC).
  -tesorter yes|no             Strip TE-derived peptides from the reference proteome with a
                               two-pass TEsorter + blastp screen (default: $TESORTER). TE
                               proteins seed false anchors genome-wide; leave this on unless
                               you are certain the proteome is already TE-free. Slow.
                               A proteome with no detectable TEs (curated reference
                               annotations are usually already TE-filtered) passes through
                               unchanged -- the screen reports 0 removed and continues.
  -cscore FLOAT                jcvi --cscore (default: $CSCORE). ~0.99 is RBH-like; lower it
                               for polyploids.
  -min_block_size N            Anchor blocks with BOTH sides >= N are kept as-is; smaller
                               blocks are merged into overlapping neighbours (default: $MIN_BLOCK_SIZE).
                               Larger values merge more aggressively, yielding fewer/bigger
                               blocks and a faster minimap2 step in Step 10.
  -stitch_gaps yes|no          Fill the gap between consecutive syntenic blocks with a
                               synthetic block, so inter-anchor intervals are aligned rather
                               than dropped (default: $STITCH_GAPS). Guarded against inversions
                               and rearrangements.
  -h, --help                   Show this help message and exit
EOF
}

# Parsing command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
    -h|--help)
        usage
        exit 0
        ;;
    -ref)
        REF="$2"
        shift; shift
        ;;
    -query)
        shift
        while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
            QUERY_GENOMES+=("$1")
            shift
        done
        ;;
    -window_size)
        WINDOW_SIZE="$2"
        shift; shift
        ;;
    -slide_size)
        SLIDE_SIZE="$2"
        shift; shift
        ;;
    -peptide)
        PROTEIN="$2"
        shift; shift
        ;;
    -threads)
        THREADS="$2"
        shift; shift
        ;;
    -processes)
        PROCESSES="$2"
        shift; shift
        ;;
    -mutation_rate)
        MUTATION_RATE="$2"
        shift; shift
        ;;
    -names)
        NAMES="$2"
        shift; shift
        ;;
    -include_mean_line)
        INCLUDE_MEAN_LINE="$2"
        shift; shift
        ;;
    -ymax)
        YMAX="$2"
        shift; shift
        ;;
    -x)
        X_TYPE="$2"
        shift; shift
        ;;
    -outn)
        OUTN="$2"
        shift; shift
        ;;
    -outs)
        OUTS="$2"
        shift; shift
        ;;
    -outc)
        OUTC="$2"
        shift; shift
        ;;
    -tesorter)
        TESORTER="$2"
        shift; shift
        ;;
    -cscore)
        CSCORE="$2"
        shift; shift
        ;;
    -min_block_size)
        MIN_BLOCK_SIZE="$2"
        shift; shift
        ;;
    -stitch_gaps)
        STITCH_GAPS="$2"
        shift; shift
        ;;
    *)
        echo "Error: Unknown option $1"
        usage
        exit 1
        ;;
    esac
done

# Check for required arguments
if [[ -z "$REF" ]]; then
    echo "Error: Reference genome (-ref) is required."
    usage
    exit 1
fi

if [[ ${#QUERY_GENOMES[@]} -eq 0 ]]; then
    echo "Error: At least one query genome (-query) is required."
    usage
    exit 1
fi

if [[ ! -s "$PROTEIN" ]]; then
    echo "Error: reference proteome not found or empty: $PROTEIN"
    echo "       Supply one with -peptide."
    exit 1
fi

if ! [[ "$MIN_BLOCK_SIZE" =~ ^[0-9]+$ ]]; then
    echo "Error: -min_block_size must be a non-negative integer (got '$MIN_BLOCK_SIZE')."
    exit 1
fi

if ! [[ "$OUTN" =~ ^[0-9]+$ ]]; then
    echo "Error: -outn must be a non-negative integer (got '$OUTN')."
    exit 1
fi

case "$STITCH_GAPS" in
    yes|no) ;;
    *) echo "Error: -stitch_gaps must be 'yes' or 'no' (got '$STITCH_GAPS')."; exit 1 ;;
esac

case "$TESORTER" in
    yes|no) ;;
    *) echo "Error: -tesorter must be 'yes' or 'no' (got '$TESORTER')."; exit 1 ;;
esac

# Fail fast on missing tools rather than part-way through a long run.
required_tools=(python miniprot bedtools bioawk samtools cd-hit diamond)
if [[ "$TESORTER" == "yes" ]]; then
    required_tools+=(TEsorter seqkit blastp makeblastdb)
fi
missing_tools=()
for tool in "${required_tools[@]}"; do
    command -v "$tool" >/dev/null 2>&1 || missing_tools+=("$tool")
done
if [[ ${#missing_tools[@]} -gt 0 ]]; then
    echo "Error: missing required tool(s): ${missing_tools[*]}"
    echo "       Install them into the active environment and re-run."
    exit 1
fi

# Assembled once and reused by liftover.py in Steps 2-4.
LIFTOVER_OPTS=(--outn "$OUTN" --outs "$OUTS" --outc "$OUTC" --cdhit)
if [[ "$TESORTER" == "yes" ]]; then
    LIFTOVER_OPTS+=(--TEsorter)
fi

# Assembled once and reused by the consolidator in Step 8.
CONSOLIDATOR_OPTS=(-t "$MIN_BLOCK_SIZE")
if [[ "$STITCH_GAPS" == "yes" ]]; then
    CONSOLIDATOR_OPTS+=(--stitch-gaps)
fi

# Collect all genomes (reference and query)
ALL_GENOMES=("$REF" "${QUERY_GENOMES[@]}")

# Derive a genome ID from a path: basename, minus one compression suffix, minus one
# FASTA suffix. MUST stay in lockstep with fasta_stem() in bin/fastaio.py -- the shell
# names the files ('{id}_mod.fa', '{id}.pep', the anchor chain) that the Python steps
# write, so a disagreement means the pipeline looks for outputs that exist under
# another name. Only known suffixes are stripped, so 'Poa.annua.v2.fa.gz' keeps its
# dots and yields 'Poa.annua.v2'.
fasta_stem() {
    local f
    f=$(basename "$1")
    case "$f" in
        *.gz|*.bgz|*.bgzf|*.bz2|*.xz|*.lzma) f="${f%.*}" ;;
    esac
    case "$f" in
        *.fa|*.fas|*.fasta|*.fna|*.ffn|*.faa|*.mfa|*.pep|*.seq) f="${f%.*}" ;;
    esac
    printf '%s' "$f"
}

# Validate the genomes before doing anything expensive. gzip/bzip2/xz inputs are read
# directly by Step 1; anything else compressed must be decompressed by hand first.
for genome in "${ALL_GENOMES[@]}"; do
    if [[ ! -s "$genome" ]]; then
        echo "Error: genome not found or empty: $genome"
        exit 1
    fi
    case "$genome" in
        *.zst|*.zstd|*.bz|*.lz4|*.Z|*.zip)
            echo "Error: unsupported compression for $genome."
            echo "       Supported: plain, .gz/.bgz/.bgzf, .bz2, .xz/.lzma. Decompress it and re-run."
            exit 1
            ;;
    esac
done

# Extract genome IDs
declare -A GENOME_IDS
for genome in "${ALL_GENOMES[@]}"; do
    GENOME_IDS["$genome"]=$(fasta_stem "$genome")
done

# Two inputs whose IDs collide would silently overwrite each other's outputs.
mapfile -t _sorted_ids < <(printf '%s\n' "${GENOME_IDS[@]}" | sort)
mapfile -t _dupe_ids < <(printf '%s\n' "${_sorted_ids[@]}" | uniq -d)
if [[ ${#_dupe_ids[@]} -gt 0 ]]; then
    echo "Error: duplicate genome ID(s): ${_dupe_ids[*]}"
    echo "       Genome IDs come from the filename; give the inputs distinct names."
    exit 1
fi

# A genome whose sequences would all be dropped by the renamer yields an empty *_mod.fa, which
# only surfaces much later as an opaque miniprot "failed to open/build the index". Catch it here,
# before cd-hit/TEsorter/miniprot spend CPU on a run that is already doomed. Reads headers only.
if ! python "$BIN_DIR/fasta_renamer_diploid.py" --preflight -genomes "${ALL_GENOMES[@]}"; then
    echo "Aborting before Step 1: no work has been done. Fix the input(s) above and re-run." >&2
    exit 1
fi


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
    if ! python "$BIN_DIR/fasta_renamer_diploid.py" -genomes "${ALL_GENOMES[@]}"; then
        echo "Error: Step 1 (fasta_renamer_diploid.py) failed. Aborting before the liftover." >&2
        exit 1
    fi
else
    echo "Step 1 - Modified fasta files and jcvi_list.txt exist. Skipping."
fi

# Steps 2-4 - Liftover the reference proteome onto each genome (liftover.py).
#
# Replaces the old raw-miniprot -> awk GFF-to-BED -> bedtools pseudo-CDS chain, matching
# synLTR/module1.py. liftover.py wraps miniprot with score (--outs) and coverage (--outc)
# filters, de-duplicates the reference proteome with cd-hit, and strips TE-derived peptides
# via a two-pass TEsorter + blastp screen. TE proteins seed false anchors genome-wide, so
# removing them is what makes the anchors trustworthy. Output is .pep (not pseudo-CDS), so
# Step 5 anchors on protein rather than nucleotide.
#
# liftover.py names its outputs after the genome file, so '{id}_mod.fa' yields '{id}_mod.pep'
# and '{id}_mod.bed'. jcvi_list.txt and the whole anchor chain key off the bare '{id}', so the
# two are renamed below. '{id}_mod.gff' is deliberately left in place: liftover.py reuses an
# existing GFF and skips the expensive miniprot run on resume.
liftover_outputs_exist=true
for genome in "${ALL_GENOMES[@]}"; do
    id="${GENOME_IDS[$genome]}"
    if [[ ! -s "${id}.pep" || ! -s "${id}.bed" ]]; then
        liftover_outputs_exist=false
        break
    fi
done

if [[ "$liftover_outputs_exist" = false ]]; then
    echo "Steps 2-4 - Liftover $PROTEIN onto each genome (tesorter=$TESORTER)"
    mod_fastas=()
    for genome in "${ALL_GENOMES[@]}"; do
        mod_fastas+=("${GENOME_IDS[$genome]}_mod.fa")
    done

    echo "Running: python $BIN_DIR/liftover.py --genome ${mod_fastas[*]} --reference $PROTEIN ${LIFTOVER_OPTS[*]} --threads $THREADS --outputs gff pep bed"
    python "$BIN_DIR/liftover.py" \
        --genome "${mod_fastas[@]}" \
        --reference "$PROTEIN" \
        "${LIFTOVER_OPTS[@]}" \
        --threads "$THREADS" \
        --outputs gff pep bed

    # Rename to the bare genome ID that jcvi_list.txt and the anchor chain expect.
    for genome in "${ALL_GENOMES[@]}"; do
        id="${GENOME_IDS[$genome]}"
        for ext in pep bed; do
            if [[ -s "${id}_mod.${ext}" ]]; then
                mv -f "${id}_mod.${ext}" "${id}.${ext}"
            fi
        done
        if [[ ! -s "${id}.pep" || ! -s "${id}.bed" ]]; then
            echo "Error: liftover produced no ${id}.pep / ${id}.bed."
            echo "       Check that $PROTEIN aligns to ${id}_mod.fa."
            exit 1
        fi
    done
else
    echo "Steps 2-4 - Liftover outputs (.pep/.bed) exist. Skipping."
fi

# Step 5 - Anchor proteins between each genome pair
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
    # --prot anchors on the lifted-over peptides via diamond_blastp (--dbtype prot) rather than
    # on nucleotide pseudo-CDS. '*.dmnd' is the diamond DB jcvi builds; it is cleared alongside
    # the legacy BLAST/LAST DBs so a failed attempt cannot leave a partial DB behind.
    echo "Step 5 - Anchor proteins (jcvi --prot, --cscore $CSCORE)"
    echo "Running: python $BIN_DIR/jcvi_diploid.py -p $PROCESSES --cpus $THREADS --prot --cscore $CSCORE"
    python "$BIN_DIR/jcvi_diploid.py" -p "$PROCESSES" --cpus "$THREADS" --prot --cscore "$CSCORE"
    rm -f *.nsq *.nin *.nhr *.ndb *.nto *.not *.ntf *.njs *.des *.sds *.tis *.ssp *.bck *.suf *.prj *.dmnd

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

        python "$BIN_DIR/jcvi_diploid_retry.py" -p "$current_p" --cpus "$current_cpus" --prot --cscore "$CSCORE"
        rm -f *.nsq *.nin *.nhr *.ndb *.nto *.not *.ntf *.njs *.des *.sds *.tis *.ssp *.bck *.suf *.prj *.dmnd
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
#
# Synteny chain, kept in lockstep with synLTR/module1.py step (5):
#     anchor_builder -> gene_coords_extractor -> subtracter -> subtracter -> consolidator
#
# Order matters. The subtracter runs FIRST so that containments are dropped and partial
# overlaps trimmed; only then does the consolidator merge and stitch. Stitching a gap is
# only meaningful between blocks that no longer overlap, so consolidating first (the old
# order) made the inter-anchor gaps ill-defined.
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

    # Step 7 - Convert the cleaned anchor files to coordinate format.
    # 'sort | uniq' is load-bearing: anchor_coord_subtracter treats each of two identical
    # lines as containing the other and removes BOTH, silently deleting the block. Exact
    # duplicates must never reach it.
    raw_coords_file="${ID1}.${ID2}.anchors.raw.coords"
    if [[ ! -s "$raw_coords_file" ]]; then
        echo "Converting $clean_anchor_file to $raw_coords_file"
        python "$BIN_DIR/gene_coords_extractor_all4.py" -mcscan "$clean_anchor_file" | sort | uniq >"$raw_coords_file"
    else
        echo "Raw coords file $raw_coords_file exists. Skipping."
    fi

    # Step 8 - Polish overlapping anchors (two passes), then consolidate and stitch gaps.
    polished_file="${ID1}.${ID2}.anchors.coords.polished"
    polished2_file="${ID1}.${ID2}.anchors.coords.polished2"
    coords_file="${ID1}.${ID2}.anchors.coords"
    if [[ ! -s "$coords_file" ]]; then
        echo "Polishing $raw_coords_file (pass 1 -> $polished_file)"
        python "$BIN_DIR/anchor_coord_subtracter.py" "$raw_coords_file" "$polished_file"

        # Second pass: the partial-overlap trim mutates records in place while iterating,
        # so one pass is not guaranteed to reach a fixed point.
        echo "Polishing $polished_file (pass 2 -> $polished2_file)"
        python "$BIN_DIR/anchor_coord_subtracter.py" "$polished_file" "$polished2_file"

        echo "Consolidating $polished2_file to $coords_file (-t $MIN_BLOCK_SIZE, stitch_gaps=$STITCH_GAPS)"
        python "$BIN_DIR/anchor_coord_consolidator.py" "${CONSOLIDATOR_OPTS[@]}" "$polished2_file" >"$coords_file"
    else
        echo "Coords file $coords_file exists. Skipping."
    fi

done <jcvi_list.txt

# Step 9 - Merge the pairwise anchor files.
# Each pair is already polished and consolidated, so the merge is the final product.
# The subtracter groups by (ref_chrom, query_chrom, strand), and those keys are unique to a
# pair, so polishing per-pair and merging is equivalent to polishing the merged file -- but
# it keeps the O(n^2) polish on small per-pair inputs.
POLISHED="all.anchors.coords.polished"
if [[ ! -s "$POLISHED" ]]; then
    coords_files=()
    while read -r line; do
        ID1=$(echo "$line" | awk '{print $1}')
        ID2=$(echo "$line" | awk '{print $2}')
        coords_files+=("${ID1}.${ID2}.anchors.coords")
    done <jcvi_list.txt

    echo "Step 9 - Merge the pairwise anchor files"
    cat "${coords_files[@]}" >"$POLISHED"
else
    echo "Step 9 - Polished coords file $POLISHED exists. Skipping."
fi

# Step 10 - Align the genomic anchors
if [[ ! -s "alignment.paf" ]]; then
    echo "Step 10 - Align the genomic anchors"
    python -u "$BIN_DIR/synmap_split.py" --timer 10m -t "$THREADS" -p "$PROCESSES" --preset "$X_TYPE" -c all.anchors.coords.polished
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
    # The four added columns on the rightmost end are (query_genome_length, target_genome_length, number_matches_and_mismatches, and k2p).
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
# GENOME_IDS already holds the stem; do not strip again, or a dotted ID loses its tail.
REF_ID="${GENOME_IDS[$REF]}"

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
