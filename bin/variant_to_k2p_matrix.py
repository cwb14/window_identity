#!/usr/bin/env python3

# Estimates k2p from minimap2 full-genome alignments. 
# Outputs a k2p sample matrix.
# Script worked correctly with a small dummy dataset.

import os
import sys
import argparse
import math
from collections import defaultdict
import glob

def is_transition(ref_base, query_base):
    ref_base = ref_base.upper()
    query_base = query_base.upper()
    transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    valid_bases = {'A', 'C', 'G', 'T'}
    if ref_base == query_base:
        return None  # Not a variant
    if ref_base in valid_bases and query_base in valid_bases:
        if transitions.get(ref_base) == query_base:
            return 'transition'
        else:
            return 'transversion'
    else:
        return None  # Invalid nucleotide, skip

def compute_k2p_distance(p, q):
    if p + q >= 0.75:
        raise ValueError('p + q >= 0.75, distance undefined')
    a = 1 - 2 * p - q
    b = 1 - 2 * q
    if a <= 0 or b <= 0:
        raise ValueError('Negative argument to log')
    d = -0.5 * math.log(a * math.sqrt(b))
    return d

def main():
    parser = argparse.ArgumentParser(description='Compute K2P distances from variant data.')
    parser.add_argument('-in', dest='variant_file', required=True, help='Input variant file')
    args = parser.parse_args()

    variant_file = args.variant_file

    counts = defaultdict(lambda: {'Ti': 0, 'Tv': 0})
    genomes = set()

    # Process variant file
    with open(variant_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            REF_ID = fields[1]
            REF_SEQ = fields[6]
            QUERY_SEQ = fields[7]
            QUERY_ID = fields[8]

            REF_GENOME = REF_ID.split('_')[0]
            QUERY_GENOME = QUERY_ID.split('_')[0]

            genomes.update([REF_GENOME, QUERY_GENOME])

            result = is_transition(REF_SEQ, QUERY_SEQ)
            if result == 'transition':
                genome_pair = tuple(sorted([REF_GENOME, QUERY_GENOME]))
                counts[genome_pair]['Ti'] += 1
            elif result == 'transversion':
                genome_pair = tuple(sorted([REF_GENOME, QUERY_GENOME]))
                counts[genome_pair]['Tv'] += 1
            # else: ignore non-variant or invalid bases

    # Read length files
    total_lengths = {}
    length_files = glob.glob('*.bed.aln.size')
    for filename in length_files:
        base = os.path.basename(filename)
        parts = base.split('.')
        if len(parts) >= 3:
            genome1 = parts[0]
            genome2 = parts[1]
            genome_pair = tuple(sorted([genome1, genome2]))
            with open(filename, 'r') as f:
                length_line = f.readline().strip()
                try:
                    length = int(length_line)
                    total_lengths[genome_pair] = length
                except ValueError:
                    sys.stderr.write(f"Warning: Invalid length in file {filename}\n")
        else:
            sys.stderr.write(f"Warning: Unrecognized filename format: {filename}\n")

    genome_list = sorted(genomes)
    distance_matrix = defaultdict(dict)

    for i in genome_list:
        for j in genome_list:
            genome_pair = tuple(sorted([i, j]))
            if i == j:
                distance_matrix[i][j] = 0.0
            else:
                Ti = counts.get(genome_pair, {}).get('Ti', 0)
                Tv = counts.get(genome_pair, {}).get('Tv', 0)
                total_length = total_lengths.get(genome_pair, None)
                if total_length is None or total_length == 0:
                    sys.stderr.write(f"Warning: Total length for genomes {genome_pair} not found or zero.\n")
                    distance_matrix[i][j] = 'NA'
                    continue
                p = Ti / total_length
                q = Tv / total_length
                try:
                    d = compute_k2p_distance(p, q)
                    distance_matrix[i][j] = d
                except ValueError as e:
                    sys.stderr.write(f"Error computing K2P distance for genomes {genome_pair}: {e}\n")
                    distance_matrix[i][j] = 'NA'

    # Output the distance matrix
    print('\t' + '\t'.join(genome_list))
    for i in genome_list:
        line = [i]
        for j in genome_list:
            value = distance_matrix[i].get(j, 'NA')
            if isinstance(value, float):
                value = f"{value:.6f}"
            line.append(str(value))
        print('\t'.join(line))

if __name__ == '__main__':
    main()
