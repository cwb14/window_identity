#!/usr/bin/env python3

import argparse
import sys
from collections import OrderedDict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate weighted divergence averages for unique ID pairs.')
    parser.add_argument('-input', required=True, help='Path to the input TSV file.')
    return parser.parse_args()

def extract_divergence(divergence_field):
    """
    Extracts the numeric divergence value from a string formatted as 'de:f:<value>'.
    """
    try:
        return float(divergence_field.split(':')[-1])
    except (IndexError, ValueError):
        raise ValueError(f"Invalid divergence format: '{divergence_field}'")

def process_file(input_file):
    """
    Processes the input file and calculates the sum of weights and weighted divergence for each unique ID pair.
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

def main():
    args = parse_arguments()
    try:
        pair_data = process_file(args.input)
    except FileNotFoundError:
        print(f"Error: File '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error while processing the file: {e}", file=sys.stderr)
        sys.exit(1)

    results = calculate_weighted_averages(pair_data)

    # Output the results
    for ID1, ID2, sum_weights, weighted_avg in results:
        # Format the weighted average with sufficient decimal places
        print(f"{ID1}\t{ID2}\t{sum_weights:.6f}\tde:f:{weighted_avg:.6f}")

if __name__ == "__main__":
    main()
