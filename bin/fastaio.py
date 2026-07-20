#!/usr/bin/env python3
"""
Shared FASTA input helpers: transparent decompression and filename-stem derivation.

Genomes are routinely distributed gzip-compressed (NCBI, Ensembl, and every assembly
hub ship '.fa.gz'), so every entry point that accepts a user-supplied FASTA reads it
through open_fasta() rather than open().

The stem rules here are the single source of truth for how an input path becomes a
genome ID. weighted_seq_ident_calc2.sh mirrors fasta_stem() in shell; the two MUST
agree, or the pipeline writes '{id}_mod.fa' under one name and looks for it under
another.
"""
import bz2
import glob
import gzip
import lzma
import os
import shutil

# Compression suffixes we can decode with the standard library alone (no new deps).
# BGZF is a valid gzip stream, so '.bgz'/'.bgzf' fall through the gzip reader.
COMPRESSION_SUFFIXES = {
    ".gz": gzip.open,
    ".bgz": gzip.open,
    ".bgzf": gzip.open,
    ".bz2": bz2.open,
    ".xz": lzma.open,
    ".lzma": lzma.open,
}

# Suffixes stripped when deriving a genome ID. Kept explicit rather than using
# splitext() so that a dotted basename ('Poa.annua.v2.fa') keeps its dots.
FASTA_SUFFIXES = (
    ".fa", ".fas", ".fasta", ".fna", ".ffn", ".faa", ".mfa", ".pep", ".seq",
)


def compression_suffix(path):
    """Return the compression suffix of path (lowercased), or '' if uncompressed."""
    ext = os.path.splitext(path)[1].lower()
    return ext if ext in COMPRESSION_SUFFIXES else ""


def is_compressed(path):
    return bool(compression_suffix(path))


def open_fasta(path, mode="rt"):
    """
    Open a FASTA for reading, decompressing on the fly when the suffix says so.

    Returns a text-mode handle. The caller should iterate it rather than read it
    whole: these are genomes.
    """
    suffix = compression_suffix(path)
    if suffix:
        return COMPRESSION_SUFFIXES[suffix](path, mode)
    return open(path, mode)


def fasta_stem(path):
    """
    Basename minus one compression suffix minus one FASTA suffix.

        Pinfirma.fa        -> Pinfirma
        Pinfirma.fa.gz     -> Pinfirma
        Poa.annua.v2.fa.gz -> Poa.annua.v2

    Only known suffixes are stripped, so a version-dotted name survives intact.
    """
    name = os.path.basename(path)
    root, ext = os.path.splitext(name)
    if ext.lower() in COMPRESSION_SUFFIXES:
        name = root
    root, ext = os.path.splitext(name)
    if ext.lower() in FASTA_SUFFIXES:
        name = root
    return name


# Cache keyed by resolved work_dir. Pool workers each build their own copy once;
# the underlying files are small (one line per genome) and never change mid-run.
_GENOME_ID_CACHE = {}


def load_genome_ids(work_dir="."):
    """
    Canonical genome IDs for a pipeline working directory, as a frozenset.

    Source precedence, first non-empty wins:
      1. genome_list.txt   -- one bare ID per line, written by fasta_renamer_diploid.py
      2. jcvi_list.txt     -- tab-separated pairs; both columns are IDs
      3. *_mod.fa          -- glob the renamer's own output

    Returns an empty frozenset when none are present, which makes accession_of()
    fall back to the legacy split('_')[0]. That keeps these scripts usable standalone,
    outside a pipeline directory.
    """
    path = os.path.join(work_dir, "genome_list.txt")
    if os.path.exists(path):
        with open(path) as fh:
            ids = {line.strip() for line in fh if line.strip()}
        if ids:
            return frozenset(ids)

    path = os.path.join(work_dir, "jcvi_list.txt")
    if os.path.exists(path):
        ids = set()
        with open(path) as fh:
            for line in fh:
                ids.update(line.split())  # whitespace split: IDs may contain dots
        if ids:
            return frozenset(ids)

    ids = set()
    for path in glob.glob(os.path.join(work_dir, "*_mod.fa")):
        name = os.path.basename(path)
        ids.add(name[: -len("_mod.fa")])
    return frozenset(ids)


def genome_ids(work_dir="."):
    """Process-cached load_genome_ids(). Use this from hot paths."""
    key = os.path.abspath(work_dir)
    if key not in _GENOME_ID_CACHE:
        _GENOME_ID_CACHE[key] = load_genome_ids(work_dir)
    return _GENOME_ID_CACHE[key]


def accession_of(seqname, known_ids):
    """
    Recover the genome ID from a renamed sequence name.

        accession_of('annuaA_chr_chr1', {'annuaA_chr'}) -> 'annuaA_chr'

    fasta_renamer_diploid.py emits '{id}_chr{N}[LETTERS]' and '{id}_sca{N}', and IDs may
    themselves contain underscores, so the naive split('_')[0] is wrong. Longest-prefix
    match against the known IDs is used rather than a suffix regex because a name like
    'annuaA_chr_chr1' has two valid-looking '_chr' boundaries and only the ID list can
    say which is real.

    Falls back to the legacy split('_')[0] when known_ids is empty or nothing matches,
    so underscore-free IDs are unaffected and standalone use still works.
    """
    if seqname in known_ids:
        return seqname
    for gid in sorted(known_ids, key=len, reverse=True):
        if seqname.startswith(gid + "_"):
            return gid
    return seqname.split("_")[0]


def materialize_plain(path, dest_dir=".", log=None):
    """
    Guarantee an uncompressed copy on disk and return its path.

    Several tools in this pipeline (cd-hit, TEsorter, makeblastdb, samtools faidx,
    bedtools getfasta) cannot read a gzip stream. Rather than teach each call site,
    decompress once here. Uncompressed inputs are returned untouched -- no copy, no
    rewrite. Pass log=print (or a logger) to narrate; decompressing a genome is slow
    enough that a silent pause looks like a hang.
    """
    if not is_compressed(path):
        return path

    plain_name = os.path.basename(path)[: -len(compression_suffix(path))]
    plain_path = os.path.join(dest_dir, plain_name)

    # Reuse an existing decompression, but only if it is newer than the archive;
    # a stale or half-written copy is worse than paying to redo the work.
    if os.path.exists(plain_path) and os.path.getsize(plain_path) > 0:
        if os.path.getmtime(plain_path) >= os.path.getmtime(path):
            if log:
                log(f"Reusing existing decompressed copy: {plain_path}")
            return plain_path

    if log:
        log(f"Decompressing {path} -> {plain_path}")
    tmp_path = plain_path + ".partial"
    with open_fasta(path, "rb") as src, open(tmp_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 22)
    os.replace(tmp_path, plain_path)  # atomic: a killed run leaves no usable stub
    return plain_path
