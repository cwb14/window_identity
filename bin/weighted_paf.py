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
        query_length = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        target_name = fields[5]
        target_length = int(fields[6])
        target_start = int(fields[7])
        target_end = int(fields[8])
    except ValueError:
        # Invalid integer conversion
        return None
    
    # Determine alignment type
    query_has_chr = 'chr' in query_name
    target_has_chr = 'chr' in target_name
    query_has_sca = 'sca' in query_name
    target_has_sca = 'sca' in target_name
    
    # Type 1: Both contain 'chr'
    if query_has_chr and target_has_chr:
        alignment_type = 1
    # Type 2: One contains 'chr' and the other contains 'sca'
    elif (query_has_chr and target_has_sca) or (query_has_sca and target_has_chr):
        alignment_type = 2
    # Type 3: Both contain 'sca' or neither contain 'chr' or 'sca' appropriately
    else:
        return None  # Ignore scaffold alignments and others
    
    # Calculate average_alignment_length
    alignment_length_query = query_end - query_start
    alignment_length_target = target_end - target_start
    average_alignment_length = (alignment_length_query + alignment_length_target) / 2
    
    if alignment_type == 1:
        average_chromosome_length = (query_length + target_length) / 2
        divergence_weight = average_alignment_length / average_chromosome_length
    elif alignment_type == 2:
        if query_has_chr:
            chromosome_length = query_length
        else:
            chromosome_length = target_length
        divergence_weight = average_alignment_length / chromosome_length
    else:
        return None  # Should not reach here
    
    # Get divergence from column21 (0-based index 20)
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
                    # Format divergence_weight to 5 decimal places
                    divergence_weight_str = f"{divergence_weight:.5f}"
                    
                    # Write to alignment_divergence.tsv
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
