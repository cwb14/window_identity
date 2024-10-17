#!/usr/bin/env python3

import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate divergence and K2P weights from PAF alignments.")
    parser.add_argument('-paf', required=True, help='Input PAF file')
    return parser.parse_args()

def process_paf_line(fields):
    """
    Process a single PAF line and return the required output fields if applicable.
    
    Returns:
        tuple: (query_name, target_name, divergence_weight, divergence, k2p) or None
    """
    if len(fields) < 21:
        # Not enough columns
        return None
    
    query_name = fields[0]
    try:
        query_start = int(fields[2])
        query_end = int(fields[3])
        target_name = fields[5]
        target_start = int(fields[7])
        target_end = int(fields[8])
        
        # Get the query and target genome lengths from the fourth and third columns from the right, respectively
        query_genome_length = int(fields[-4])
        target_genome_length = int(fields[-3])
    except ValueError:
        # Invalid integer conversion
        return None
    
    # Calculate average_alignment_length
    alignment_length_query = query_end - query_start
    alignment_length_target = target_end - target_start
    average_alignment_length = (alignment_length_query + alignment_length_target) / 2
    
    # Calculate avg_genome_length using query_genome_length and target_genome_length
    avg_genome_length = (query_genome_length + target_genome_length) / 2
    
    # Calculate divergence_weight
    divergence_weight = average_alignment_length / avg_genome_length
    
    # Get divergence from column 21 (0-based index 20)
    divergence = fields[20]
    
    # Get K2P from the rightmost field
    k2p = fields[-1]
    
    return (query_name, target_name, divergence_weight, divergence, k2p)

def main():
    args = parse_arguments()
    
    try:
        with open(args.paf, 'r') as paf_file, \
             open("alignment_de.tsv", 'w') as divergence_file, \
             open("alignment_k2p.tsv", 'w') as k2p_file:
            
            for line in paf_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip empty lines and comments
                fields = line.split('\t')
                result = process_paf_line(fields)
                if result:
                    query_name, target_name, divergence_weight, divergence, k2p = result
                    # Format divergence_weight to 12 decimal places
                    divergence_weight_str = f"{divergence_weight:.12f}"
                    
                    # Write to alignment_de.tsv
                    divergence_file.write(f"{query_name}\t{target_name}\t{divergence_weight_str}\t{divergence}\n")
                    
                    # Write to alignment_k2p.tsv
                    k2p_file.write(f"{query_name}\t{target_name}\t{divergence_weight_str}\t{k2p}\n")
    
    except FileNotFoundError:
        sys.stderr.write(f"Error: File '{args.paf}' not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An error occurred: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END
