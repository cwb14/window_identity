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
# Step 10 alignment. ALIGNER=minimap2 is the historical default and is unchanged.
ALIGNER="minimap2"
# block    = span each syntenic block end to end (historical).
# genepair = span adjacent gene pairs within a block; skips the consolidator, whose 15kb
#            merge and gap-stitching would undo the finer granularity.
PARTITION="block"
# Segments more length-skewed than this are skipped and logged; 0 disables. Borrowed from
# riparian.py's --max-len-ratio default (its rescue logic is deliberately NOT borrowed --
# it exists to keep skewed blocks for plot completeness, the opposite of what we want).
MAX_LEN_RATIO=5
# Pairwise Ks on the syntenic anchors (Steps 20-21), via ParaAT -> KaKs_Calculator.
KAKS="yes"
KAKS_METHOD="YN"
KS_RATE=1.5e-8
KS_MAX=2.0

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
  -aligner minimap2|last|wfa   Step 10 aligner (default: $ALIGNER). 'wfa' is a global
                               wavefront aligner (BiWFA, no heuristic) that returns one
                               end-to-end record per segment; 'last' is lastal/last-split.
  -partition block|genepair    Step 10 segment granularity (default: $PARTITION).
                               'genepair' spans adjacent gene pairs and skips the
                               consolidator.
  -max_len_ratio N             Skip segments whose length ratio exceeds N; 0 disables
                               (default: $MAX_LEN_RATIO)
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
  -kaks yes|no                 Estimate pairwise Ks on the syntenic anchors with ParaAT +
                               KaKs_Calculator, then build the Ks density plot, distance
                               matrix, and tree (default: $KAKS). Adds the in-frame CDS to
                               the liftover and fetches/builds both tools on first use.
  -kaks_method METHOD          KaKs_Calculator method: NG, LWL, LPB, MLWL, MLPB, GY, YN,
                               MYN, MS, or MA (default: $KAKS_METHOD). YN is the standard for
                               Ks distributions and is fast enough for whole-proteome anchor
                               sets; MA is more accurate but much slower. 'ALL' is not
                               supported -- it emits one row per method, which would mix
                               methods into a single median.
  -ks_rate RATE                Synonymous substitution rate for the Ks divergence-time tree
                               (default: $KS_RATE). Distinct from -mutation_rate, which is a
                               genome-wide nucleotide rate and calibrates the K2P tree.
  -ks_max FLOAT                Drop gene pairs with Ks >= this before taking the median
                               (default: $KS_MAX). Above ~2, Ks is saturated and unstable.
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
    -aligner)
        ALIGNER="$2"
        shift; shift
        ;;
    -partition)
        PARTITION="$2"
        shift; shift
        ;;
    -max_len_ratio)
        MAX_LEN_RATIO="$2"
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
    -kaks)
        KAKS="$2"
        shift; shift
        ;;
    -kaks_method)
        KAKS_METHOD="$2"
        shift; shift
        ;;
    -ks_rate)
        KS_RATE="$2"
        shift; shift
        ;;
    -ks_max)
        KS_MAX="$2"
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
case "$ALIGNER" in
    minimap2|last|wfa) ;;
    *) echo "Error: -aligner must be minimap2, last or wfa (got '$ALIGNER')."; exit 1 ;;
esac
case "$PARTITION" in
    block|genepair) ;;
    *) echo "Error: -partition must be block or genepair (got '$PARTITION')."; exit 1 ;;
esac

case "$TESORTER" in
    yes|no) ;;
    *) echo "Error: -tesorter must be 'yes' or 'no' (got '$TESORTER')."; exit 1 ;;
esac

case "$KAKS" in
    yes|no) ;;
    *) echo "Error: -kaks must be 'yes' or 'no' (got '$KAKS')."; exit 1 ;;
esac

# 'ALL' is deliberately excluded: KaKs_Calculator emits one row per method under ALL, and
# ks_summary.py would take the median across a mixture of methods.
case "$KAKS_METHOD" in
    NG|LWL|LPB|MLWL|MLPB|GY|YN|MYN|MS|MA) ;;
    *)
        echo "Error: -kaks_method must be one of NG LWL LPB MLWL MLPB GY YN MYN MS MA (got '$KAKS_METHOD')."
        exit 1
        ;;
esac

# Fail fast on missing tools rather than part-way through a long run.
required_tools=(python miniprot bedtools bioawk samtools cd-hit diamond)
if [[ "$TESORTER" == "yes" ]]; then
    required_tools+=(TEsorter seqkit blastp makeblastdb)
fi
if [[ "$KAKS" == "yes" ]]; then
    # ParaAT.pl is perl and shells out to an aligner. mafft is the one to use: ParaAT's
    # muscle command line is muscle-v3 syntax and silently breaks against muscle v5, and
    # clustalw2/t_coffee are rarely installed. git/make are only needed if the toolchain
    # has to be built, so setup_kaks_tools.sh checks for those itself.
    required_tools+=(perl mafft)
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

# 'inframe' runs cds_walker.py to emit frame-corrected CDS ('{id}_mod.cds.inframe') for the
# lifted-over genes. That is the nucleotide input ParaAT needs to back-translate the peptide
# alignments in Step 20; without it there is nothing to compute Ks from.
LIFTOVER_OUTPUTS=(gff pep bed)
LIFTOVER_EXTS=(pep bed)
if [[ "$KAKS" == "yes" ]]; then
    LIFTOVER_OUTPUTS+=(inframe)
    LIFTOVER_EXTS+=(cds.inframe)
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

# The liftover renames its per-genome peptides to '{id}.pep' (see Steps 2-4). If the
# reference proteome is itself named '{id}.pep' for one of the genomes -- '-ref Athal.fa
# -peptide Athal.pep' -- that rename overwrites the reference with the lifted-over peptides.
# The first run still succeeds, because the liftover reads the reference before the rename
# clobbers it, but the input is destroyed and any later run in that directory silently uses
# a self-referential proteome. Warn rather than abort: a directory where this already
# happened still resumes correctly from the cached '{id}_mod.gff'.
for genome in "${ALL_GENOMES[@]}"; do
    if [[ "$(basename "$PROTEIN")" == "${GENOME_IDS[$genome]}.pep" ]]; then
        echo "WARNING: the reference proteome '$PROTEIN' has the same name as the liftover"
        echo "         output for genome '${GENOME_IDS[$genome]}', and will be OVERWRITTEN by it."
        echo "         Keep a copy, or rename it (e.g. '${GENOME_IDS[$genome]}_ref.pep') and re-run."
        break
    fi
done

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
# The resume check must cover every extension requested, '.cds.inframe' included. A run
# directory left over from a -kaks no run already has .pep/.bed, so checking only those
# would skip the liftover forever and Step 20 would never find its CDS. Re-running is cheap:
# liftover.py reuses the existing '{id}_mod.gff' and the cached TEsorter screen, so neither
# miniprot nor TEsorter repeats -- only the peptide/CDS extraction does.
liftover_outputs_exist=true
for genome in "${ALL_GENOMES[@]}"; do
    id="${GENOME_IDS[$genome]}"
    for ext in "${LIFTOVER_EXTS[@]}"; do
        if [[ ! -s "${id}.${ext}" ]]; then
            liftover_outputs_exist=false
            break 2
        fi
    done
done

if [[ "$liftover_outputs_exist" = false ]]; then
    echo "Steps 2-4 - Liftover $PROTEIN onto each genome (tesorter=$TESORTER, kaks=$KAKS)"
    mod_fastas=()
    for genome in "${ALL_GENOMES[@]}"; do
        mod_fastas+=("${GENOME_IDS[$genome]}_mod.fa")
    done

    echo "Running: python $BIN_DIR/liftover.py --genome ${mod_fastas[*]} --reference $PROTEIN ${LIFTOVER_OPTS[*]} --threads $THREADS --outputs ${LIFTOVER_OUTPUTS[*]}"
    python "$BIN_DIR/liftover.py" \
        --genome "${mod_fastas[@]}" \
        --reference "$PROTEIN" \
        "${LIFTOVER_OPTS[@]}" \
        --threads "$THREADS" \
        --outputs "${LIFTOVER_OUTPUTS[@]}"

    # Rename to the bare genome ID that jcvi_list.txt and the anchor chain expect.
    for genome in "${ALL_GENOMES[@]}"; do
        id="${GENOME_IDS[$genome]}"
        for ext in "${LIFTOVER_EXTS[@]}"; do
            if [[ -s "${id}_mod.${ext}" ]]; then
                mv -f "${id}_mod.${ext}" "${id}.${ext}"
            fi
            if [[ ! -s "${id}.${ext}" ]]; then
                echo "Error: liftover produced no ${id}.${ext}."
                echo "       Check that $PROTEIN aligns to ${id}_mod.fa."
                exit 1
            fi
        done
    done
else
    echo "Steps 2-4 - Liftover outputs (${LIFTOVER_EXTS[*]}) exist. Skipping."
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
        echo "Converting $clean_anchor_file to $raw_coords_file (partition=$PARTITION)"
        if [[ "$PARTITION" == "genepair" ]]; then
            extractor="gene_coords_extractor_all4_pairs.py"
        else
            extractor="gene_coords_extractor_all4.py"
        fi
        python "$BIN_DIR/$extractor" -mcscan "$clean_anchor_file" | sort | uniq >"$raw_coords_file"
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

        if [[ "$PARTITION" == "genepair" ]]; then
            # The consolidator merges to >=MIN_BLOCK_SIZE and stitches gaps, which would
            # rebuild the very blocks gene-pair partitioning exists to avoid.
            echo "Partition=genepair: skipping the consolidator, using $polished2_file as-is"
            cp "$polished2_file" "$coords_file"
        else
            echo "Consolidating $polished2_file to $coords_file (-t $MIN_BLOCK_SIZE, stitch_gaps=$STITCH_GAPS)"
            python "$BIN_DIR/anchor_coord_consolidator.py" "${CONSOLIDATOR_OPTS[@]}" "$polished2_file" >"$coords_file"
        fi
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
    python -u "$BIN_DIR/synmap_split.py" --timer 10m -t "$THREADS" -p "$PROCESSES" \
        --preset "$X_TYPE" --aligner "$ALIGNER" --max-len-ratio "$MAX_LEN_RATIO" \
        -c all.anchors.coords.polished
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
# --prefix defaults to k2p_matrix, so this pass picks up k2p_matrix*.tsv only and the Ks
# matrices written in Step 21 are left to their own pass.
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

if [[ "$KAKS" != "yes" ]]; then
    echo "Steps 20-21 - Ks estimation disabled (-kaks no). Done."
    exit 0
fi

# Step 20 - Pairwise Ks on the syntenic anchors (ParaAT -> KaKs_Calculator).
#
# The anchors are already exactly what ParaAT wants. '{ID1}.{ID2}.clean.anchors' is two
# columns of gene IDs -- ParaAT's homolog format -- once MCscan's '###' block separators are
# stripped. And because liftover.py prefixes every gene ID with its genome ('Athal_mod000001'),
# the two genomes' .pep and .cds.inframe can simply be concatenated with no risk of an ID
# collision silently pairing a gene with itself.
#
# Ks is computed off the liftover of a single reference proteome onto every genome, so it is
# independent of Steps 10-19 (the minimap2/K2P branch) and reads only the Step 2-5 outputs.
echo "Step 20 - Pairwise Ks (ParaAT + KaKs_Calculator, method=$KAKS_METHOD)"

TOOLS_DIR="$SCRIPT_DIR/tools"
if ! bash "$BIN_DIR/setup_kaks_tools.sh" "$TOOLS_DIR"; then
    echo "Error: could not set up the ParaAT/KaKs_Calculator toolchain in $TOOLS_DIR." >&2
    exit 1
fi
# ParaAT.pl resolves Epal2nal.pl, mafft, and KaKs_Calculator by walking PATH and testing -x.
export PATH="$TOOLS_DIR/ParaAT:$TOOLS_DIR/shim:$PATH"
# Read by the shim. ParaAT hardcodes its KaKs command line, so -m cannot be passed through it.
export KAKS_METHOD

# ParaAT takes the processor count as a FILE, not a number, and re-reads it between batches
# so it can be retuned mid-run.
#
# Three constraints on the paths below, all forced by ParaAT.pl internals -- do not
# "tidy" them into absolute paths:
#   * ParaAT chdir()s into its output folder, then re-reads the processor file as
#     '../$ProcessFile'. So -p must be RELATIVE to the launch directory, and the -o folder
#     must sit exactly ONE level below it, or '../' misses and ParaAT silently falls back to
#     a default thread count.
#   * -h/-a/-n are read before that chdir, so they resolve against the launch directory.
PROC_FILE="paraat.proc"
echo "$THREADS" >"$PROC_FILE"

while read -r line; do
    ID1=$(echo "$line" | awk '{print $1}')
    ID2=$(echo "$line" | awk '{print $2}')
    kaks_tsv="${ID1}.${ID2}.kaks.tsv"

    if [[ -s "$kaks_tsv" ]]; then
        echo "Ks table $kaks_tsv exists. Skipping."
        continue
    fi

    homologs="${ID1}.${ID2}.homologs"
    pair_pep="${ID1}.${ID2}.paraat.pep"
    pair_cds="${ID1}.${ID2}.paraat.cds"
    outdir="kaks_${ID1}_${ID2}"

    grep -v '^###' "${ID1}.${ID2}.clean.anchors" | cut -f1,2 >"$homologs"
    cat "${ID1}.pep" "${ID2}.pep" >"$pair_pep"
    cat "${ID1}.cds.inframe" "${ID2}.cds.inframe" >"$pair_cds"

    echo "Running ParaAT on $(wc -l <"$homologs") anchor pairs from $ID1/$ID2 -> $outdir"
    # A partial ParaAT run leaves stale per-gene files that the next pass would fold into the
    # merged table, so start clean.
    rm -rf "$outdir"
    ParaAT.pl -h "$homologs" -a "$pair_pep" -n "$pair_cds" \
        -p "$PROC_FILE" -m mafft -f axt -k -o "$outdir"

    # ParaAT writes one .kaks per gene pair. Merge them into a single table, keeping one
    # header. Anchors whose peptide carries an internal stop (miniprot pseudogene calls) fail
    # in Epal2nal and are simply absent here, so the count is reported rather than assumed.
    first_kaks=$(find "$outdir" -name '*.kaks' -print -quit)
    if [[ -z "$first_kaks" ]]; then
        echo "Error: ParaAT/KaKs produced no .kaks files in $outdir." >&2
        echo "       Check $outdir/msg.msa and $outdir/msg.kaks." >&2
        exit 1
    fi
    # 'xargs -r' matters: with no matches, a bare 'xargs cat'/'xargs awk' can run the command
    # with no file operands, which makes it read stdin and block. FNR resets per file, so the
    # header skip stays correct even when xargs splits the list across several awk calls.
    head -1 "$first_kaks" >"$kaks_tsv"
    find "$outdir" -name '*.kaks' -print0 | xargs -0 -r awk 'FNR > 1 && NF {print}' >>"$kaks_tsv"
    echo "Wrote $kaks_tsv ($(($(wc -l <"$kaks_tsv") - 1)) of $(wc -l <"$homologs") anchor pairs scored)"

    # Keep the codon alignments as one file, then drop the directory. ParaAT emits one .axt
    # and one .kaks per gene pair, which is a few hundred thousand inodes across a run.
    find "$outdir" -name '*.cds_aln.axt' -print0 | xargs -0 -r cat >"${ID1}.${ID2}.axt"
    rm -rf "$outdir" "$pair_pep" "$pair_cds"
done <jcvi_list.txt

# Step 21 - Summarise Ks, then build the density plot, distance matrix, and tree.
if [[ ! -s "ks_genome.tsv" ]]; then
    echo "Step 21 - Summarising Ks (median per genome pair, 0 < Ks < $KS_MAX)"
    python "$BIN_DIR/ks_summary.py" \
        -list jcvi_list.txt \
        -max_ks "$KS_MAX" \
        -genome_out ks_genome.tsv \
        -long_out ks_all.tsv
else
    echo "Step 21 - ks_genome.tsv exists. Skipping."
fi

if [[ ! -s "ks_matrix.tsv" ]]; then
    echo "Building the Ks matrix"
    python "$BIN_DIR/matrix_builder.py" -in ks_genome.tsv >ks_matrix.tsv
else
    echo "ks_matrix.tsv exists. Skipping."
fi

echo "Plotting the Ks density distribution"
ks_plot_cmd="Rscript $BIN_DIR/ks_density_plotter.R -i ks_all.tsv -o ks_density.pdf --max_ks $KS_MAX"
echo "Running: $ks_plot_cmd"
eval "$ks_plot_cmd"

# Same UPGMA machinery as Step 19, pointed at ks_matrix*.tsv and calibrated with the
# synonymous rate rather than the genome-wide nucleotide rate.
echo "Building the Ks tree (UPGMA, synonymous rate $KS_RATE)"
ks_tree_cmd="Rscript $BIN_DIR/upgma.R --method upgma --prefix ks_matrix --xlab Ks --mutation_rate $KS_RATE"
if [[ -n "$NAMES" ]]; then
    ks_tree_cmd+=" --name_key $NAMES"
fi
echo "Running: $ks_tree_cmd"
eval "$ks_tree_cmd"
