#!/usr/bin/env python3

import argparse
import os
import sys

### coordinates between all pairs of genes. within and between blocks (although prolly would need to exclude between... if it works).
# cat Col0.Alt5.clean.anchors | grep -v '###' > Col0.Alt5.clean.anchors.mod
# python window_identity/bin/gene_coords_extractor_all4_pairs.py -mcscan Col0.Alt5.clean.anchors.mod > Col0.Alt5.clean.anchors.mod.coords

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
        bed_data = modify_bed_file(bed_file)
        merged_data.extend(bed_data)
    return merged_data


def get_gene_coords(merged_data):
    """
    Parse merged BED lines into a dict mapping geneID -> (chrom, start, end).
    """
    gene_coords = {}
    for line in merged_data:
        chrom, start, end, geneID = line.split("\t")
        gene_coords[geneID] = (chrom, int(start), int(end))
    return gene_coords


def read_gene_pairs(gene_pairs_filename):
    """
    Read a file where each line has two gene IDs (tab-separated).
    Returns a list of (gene1, gene2) tuples.
    """
    pairs = []
    with open(gene_pairs_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                sys.stderr.write(f"Skipping malformed line: {line}\n")
                continue
            gene1, gene2 = parts
            pairs.append((gene1, gene2))
    return pairs


def process_pairs(pairs, gene_coords):
    """
    For each adjacent pair of lines, compute the spanning coordinates
    for gene1->gene1_next and gene2->gene2_next, plus directionality.
    """
    output_strings = []
    if len(pairs) < 2:
        return output_strings

    for i in range(len(pairs) - 1):
        geneA1, geneA2 = pairs[i]
        geneB1, geneB2 = pairs[i + 1]

        coords1 = gene_coords.get(geneA1)
        coords2 = gene_coords.get(geneA2)
        coords3 = gene_coords.get(geneB1)
        coords4 = gene_coords.get(geneB2)

        if None in (coords1, coords2, coords3, coords4):
            sys.stderr.write(
                f"Missing coordinates for genes: {geneA1},{geneA2},{geneB1},{geneB2}\n"
            )
            continue

        # Determine directionality based on start positions
        direction1 = "+" if coords1[1] < coords3[1] else "-"
        direction2 = "+" if coords2[1] < coords4[1] else "-"
        directionality = "+" if direction1 == direction2 else "-"

        # Build spanning coordinate ranges
        range1_start = min(coords1[1], coords3[1])
        range1_end = max(coords1[2], coords3[2])
        range2_start = min(coords2[1], coords4[1])
        range2_end = max(coords2[2], coords4[2])

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
        description='Extract gene coordinates and determine directionality between successive gene pairs.'
    )
    parser.add_argument(
        '-mcscan', required=True,
        help='Input file containing one gene-pair per line (tab-separated).'
    )
    args = parser.parse_args()

    bed_files = derive_bed_files(args.mcscan)
    merged_data = read_bed_data(bed_files)
    gene_coords = get_gene_coords(merged_data)
    pairs = read_gene_pairs(args.mcscan)
    output_strings = process_pairs(pairs, gene_coords)

    for output_str in output_strings:
        print(output_str)

if __name__ == "__main__":
    main()
