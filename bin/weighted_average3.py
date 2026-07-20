#!/usr/bin/env python3

import argparse
import sys
from collections import OrderedDict
import fastaio

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate weighted divergence averages for unique ID pairs at chromosome and genome levels.')
    parser.add_argument('-input', required=True, help='Path to the input TSV file.')
    parser.add_argument('-chrom_out', required=True, help='Path to the chromosome-level output TSV file.')
    parser.add_argument('-genome_out', required=True, help='Path to the genome-level output TSV file.')
    return parser.parse_args()

def extract_divergence(divergence_field):
    """
    Extracts the numeric divergence value from a string.
    Supports both 'de:f:<value>' and plain numeric values.
    """
    if divergence_field.startswith('de:f:'):
        try:
            return float(divergence_field.split(':')[-1])
        except (IndexError, ValueError):
            raise ValueError(f"Invalid divergence format: '{divergence_field}'")
    else:
        try:
            return float(divergence_field)
        except ValueError:
            raise ValueError(f"Invalid divergence format: '{divergence_field}'")


def process_file(input_file):
    """
    Processes the input file and calculates the sum of weights and weighted divergence for each unique chromosome ID pair.
    """
    # Use OrderedDict to preserve the order of first occurrence of each unique pair
    pair_data = OrderedDict()

    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            parts = line.split('\t')
            if len(parts) != 4:
                print(f"Warning: Line {line_num} does not have exactly 4 columns. Skipping.", file=sys.stderr)
                continue

            ID1, ID2, weight_str, divergence_str = parts

            try:
                weight = float(weight_str)
            except ValueError:
                print(f"Warning: Invalid weight '{weight_str}' on line {line_num}. Skipping.", file=sys.stderr)
                continue

            try:
                divergence = extract_divergence(divergence_str)
            except ValueError as e:
                print(f"Warning: {e} on line {line_num}. Skipping.", file=sys.stderr)
                continue

            # Create a sorted tuple of IDs to handle unordered pairs
            sorted_ids = tuple(sorted([ID1, ID2]))

            if sorted_ids not in pair_data:
                pair_data[sorted_ids] = {'sum_weights': 0.0, 'sum_weight_divergence': 0.0}

            pair_data[sorted_ids]['sum_weights'] += weight
            pair_data[sorted_ids]['sum_weight_divergence'] += weight * divergence

    return pair_data

def calculate_weighted_averages(pair_data):
    """
    Calculates the weighted divergence average for each ID pair.
    """
    results = []
    for (ID1, ID2), data in pair_data.items():
        sum_weights = data['sum_weights']
        sum_weight_divergence = data['sum_weight_divergence']
        if sum_weights == 0:
            weighted_avg = 0.0
            print(f"Warning: Sum of weights for pair ({ID1}, {ID2}) is zero. Setting weighted average to 0.", file=sys.stderr)
        else:
            weighted_avg = sum_weight_divergence / sum_weights
        results.append((ID1, ID2, sum_weights, weighted_avg))
    return results

def write_output(output_file, results):
    """
    Writes the calculated weighted averages to the output file.
    """
    with open(output_file, 'w') as out:
        for ID1, ID2, sum_weights, weighted_avg in results:
            # Format the weighted average with sufficient decimal places
            out.write(f"{ID1}\t{ID2}\t{sum_weights:.6f}\tde:f:{weighted_avg:.6f}\n")

# -------------------- Genome-Level Processing -------------------- #

def extract_genome_id(full_id):
    """
    Recover the genome ID, resolving against the canonical list so that IDs
    containing underscores are not truncated.
    """
    return fastaio.accession_of(full_id, fastaio.genome_ids())

def process_genome_file(input_file):
    """
    Processes the input file and calculates the sum of weights and weighted divergence for each unique genome ID pair.
    """
    genome_data = OrderedDict()

    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            parts = line.split('\t')
            if len(parts) != 4:
                print(f"Warning: Line {line_num} does not have exactly 4 columns. Skipping.", file=sys.stderr)
                continue

            ID1_full, ID2_full, weight_str, divergence_str = parts

            try:
                weight = float(weight_str)
            except ValueError:
                print(f"Warning: Invalid weight '{weight_str}' on line {line_num}. Skipping.", file=sys.stderr)
                continue

            try:
                divergence = extract_divergence(divergence_str)
            except ValueError as e:
                print(f"Warning: {e} on line {line_num}. Skipping.", file=sys.stderr)
                continue

            # Extract genome IDs
            genome1 = extract_genome_id(ID1_full)
            genome2 = extract_genome_id(ID2_full)

            # Create a sorted tuple of genome IDs to handle unordered pairs
            sorted_genomes = tuple(sorted([genome1, genome2]))

            if sorted_genomes not in genome_data:
                genome_data[sorted_genomes] = {'sum_weights': 0.0, 'sum_weight_divergence': 0.0}

            genome_data[sorted_genomes]['sum_weights'] += weight
            genome_data[sorted_genomes]['sum_weight_divergence'] += weight * divergence

    return genome_data

# -------------------- Main Execution -------------------- #

def main():
    args = parse_arguments()
    
    # Process chromosome-level data
    try:
        pair_data = process_file(args.input)
    except FileNotFoundError:
        print(f"Error: File '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error while processing the chromosome data: {e}", file=sys.stderr)
        sys.exit(1)

    chrom_results = calculate_weighted_averages(pair_data)
    write_output(args.chrom_out, chrom_results)

    # Process genome-level data
    try:
        genome_data = process_genome_file(args.input)
    except FileNotFoundError:
        print(f"Error: File '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error while processing the genome data: {e}", file=sys.stderr)
        sys.exit(1)

    genome_results = calculate_weighted_averages(genome_data)
    write_output(args.genome_out, genome_results)

if __name__ == "__main__":
    main()
