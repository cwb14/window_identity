#!/usr/bin/env python3

import argparse
import sys

from gene_coords_extractor_all4 import (
    directionality as _directionality,
    get_gene_coords as _get_gene_coords,
)
import os

def modify_bed_file(file_name):
    """
    Read a BED file and return a list of stripped lines.
    """
    modified_lines = []
    with open(file_name, 'r') as f:
        for line in f:
            modified_lines.append(line.strip())
    return modified_lines


def read_bed_data(bed_files):
    """
    Merge data from multiple BED files into a single list of lines.
    """
    merged_data = []
    for bed_file in bed_files:
        if not os.path.exists(bed_file):
            sys.stderr.write(f"BED file not found: {bed_file}\n")
            continue
        bed_data = modify_bed_file(bed_file)
        merged_data.extend(bed_data)
    return merged_data


def get_gene_coords(merged_data):
    """geneID -> (chrom, start, end, strand).

    Delegates to the block extractor's parser so both partitioning modes share one
    strand-aware implementation; block orientation depends on column 6.
    """
    return _get_gene_coords(merged_data)


def read_gene_pairs(gene_pairs_filename):
    """
    Read a file where each line has two gene IDs (tab- or space-separated),
    with blocks separated by lines exactly equal to '###'.
    Returns a list of blocks, each block is a list of (gene1, gene2) tuples.
    """
    blocks = []
    current_block = []
    with open(gene_pairs_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '###':
                # end of current block
                if current_block:
                    blocks.append(current_block)
                    current_block = []
                continue
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                sys.stderr.write(f"Skipping malformed line: {line}\n")
                continue
            gene1, gene2 = parts
            current_block.append((gene1, gene2))
    # append any trailing block
    if current_block:
        blocks.append(current_block)
    return blocks


def process_pairs(blocks, gene_coords):
    """
    For each block, span adjacent pairs within that block only.
    """
    output_strings = []
    for block in blocks:
        if len(block) < 2:
            continue
        for i in range(len(block) - 1):
            geneA1, geneA2 = block[i]
            geneB1, geneB2 = block[i + 1]

            coords1 = gene_coords.get(geneA1)
            coords2 = gene_coords.get(geneA2)
            coords3 = gene_coords.get(geneB1)
            coords4 = gene_coords.get(geneB2)

            if None in (coords1, coords2, coords3, coords4):
                sys.stderr.write(
                    f"Missing coordinates for genes: {geneA1},{geneA2},{geneB1},{geneB2}\n"
                )
                continue

            # Orientation from relative gene strand, not gene order -- see
            # gene_coords_extractor_all4.directionality() for why the order rule was wrong.
            directionality = _directionality(coords1, coords2, coords3, coords4)
            if directionality is None:
                directionality = "+" if coords1[3] == coords2[3] else "-"

            # Build spanning coordinate ranges
            range1_start = min(coords1[1], coords3[1])
            range1_end   = max(coords1[2], coords3[2])
            range2_start = min(coords2[1], coords4[1])
            range2_end   = max(coords2[2], coords4[2])

            coord_str1 = f"{coords1[0]}:{range1_start}..{range1_end}"
            coord_str2 = f"{coords2[0]}:{range2_start}..{range2_end}"

            output_strings.append(f"{coord_str1}\t{coord_str2}\t{directionality}")

    return output_strings


def derive_bed_files(mcscan_filename):
    """
    Given an input filename like "species1.species2",
    derive the two BED file names: species1.bed and species2.bed.
    """
    name_parts = mcscan_filename.split('.')[:2]
    return [f"{name_parts[0]}.bed", f"{name_parts[1]}.bed"]


def main():
    parser = argparse.ArgumentParser(
        description='Extract gene coordinates and determine directionality between successive gene pairs (within blocks).'
    )
    parser.add_argument(
        '-mcscan', required=True,
        help='Input file containing gene-pairs, with blocks separated by lines "###".'
    )
    args = parser.parse_args()

    bed_files = derive_bed_files(args.mcscan)
    merged_data = read_bed_data(bed_files)
    gene_coords = get_gene_coords(merged_data)
    blocks = read_gene_pairs(args.mcscan)
    output_strings = process_pairs(blocks, gene_coords)

    for output_str in output_strings:
        print(output_str)


if __name__ == "__main__":
    main()
