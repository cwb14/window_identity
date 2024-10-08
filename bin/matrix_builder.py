#!/usr/bin/env python3

import argparse
import sys
from collections import OrderedDict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Build a pairwise divergence matrix from input TSV file.')
    parser.add_argument('-in', '--input', required=True, help='Input TSV file path')
    return parser.parse_args()

def read_input(file_path):
    ids_order = []
    divergences = {}
    
    with open(file_path, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # Skip empty lines or comments
            parts = line.split('\t')
            if len(parts) < 4:
                continue  # Skip malformed lines
            id1, id2, weight, divergence = parts
            # Extract the numeric part of divergence
            try:
                div_value = float(divergence.split(':')[-1])
            except ValueError:
                div_value = 0.0  # Default to 0.0 if parsing fails
            
            # Add IDs to order list if not already present
            for id_ in [id1, id2]:
                if id_ not in ids_order:
                    ids_order.append(id_)
            
            # Initialize nested dictionaries
            if id1 not in divergences:
                divergences[id1] = {}
            if id2 not in divergences:
                divergences[id2] = {}
            
            # Since the matrix is symmetric, assign both [id1][id2] and [id2][id1]
            divergences[id1][id2] = div_value
            divergences[id2][id1] = div_value
    
    return ids_order, divergences

def build_matrix(ids_order, divergences):
    # Initialize matrix with zeros for self comparisons
    matrix = OrderedDict()
    for id1 in ids_order:
        matrix[id1] = OrderedDict()
        for id2 in ids_order:
            if id1 == id2:
                matrix[id1][id2] = 0.0
            else:
                matrix[id1][id2] = divergences.get(id1, {}).get(id2, 0.0)
    return matrix

def print_matrix(ids_order, matrix):
    # Print header
    header = [''] + ids_order
    print('\t'.join(header))
    
    # Print each row
    for id1 in ids_order:
        row = [id1]
        for id2 in ids_order:
            row.append(f"{matrix[id1][id2]:.6f}")
        print('\t'.join(row))

def main():
    args = parse_arguments()
    ids_order, divergences = read_input(args.input)
    matrix = build_matrix(ids_order, divergences)
    print_matrix(ids_order, matrix)

if __name__ == '__main__':
    main()
