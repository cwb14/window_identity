#!/usr/bin/env python3

import argparse
import sys
import fastaio

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate divergence and K2P weights from PAF alignments.")
    parser.add_argument('-paf', required=True, help='Input PAF file')
    return parser.parse_args()

def process_paf_line(fields):
    """
    Process a single PAF line and return the required output fields if applicable.

    Returns:
        dict or None: {
            'query_name': str,
            'target_name': str,
            'query_accession_name': str,
            'target_accession_name': str,
            'line_match_and_mismatch': float,
            'divergence': str,
            'k2p': str
        }
    """
    if len(fields) < 12:
        return None

    query_name = fields[0]
    target_name = fields[5]

    _ids = fastaio.genome_ids()
    query_accession_name = fastaio.accession_of(query_name, _ids)
    target_accession_name = fastaio.accession_of(target_name, _ids)

    # A trailing sd:Z: segment tag may be present (added by paf_emit for resume support).
    # Anchor the positional read against the end of the *untagged* row.
    tail = fields[:-1] if fields[-1].startswith("sd:Z:") else fields
    try:
        line_match_and_mismatch = float(tail[-2])
    except ValueError:
        return None

    # Get divergence from optional fields
    divergence = None
    for field in fields[12:-4]:  # Optional fields excluding the last four numeric fields
        if field.startswith('de:f:'):
            # Handle possible concatenation without tab
            if 'rl:i:' in field:
                divergence_value = field.split('de:f:')[1].split('rl:i:')[0]
                divergence = divergence_value.strip()
            else:
                divergence = field.split('de:f:')[1]
            break
    if divergence is None:
        divergence = '0'  # Default value if not found

    # Get K2P from the rightmost field
    k2p = fields[-1]

    return {
        'query_name': query_name,
        'target_name': target_name,
        'query_accession_name': query_accession_name,
        'target_accession_name': target_accession_name,
        'line_match_and_mismatch': line_match_and_mismatch,
        'divergence': divergence,
        'k2p': k2p
    }

def main():
    args = parse_arguments()

    try:
        with open(args.paf, 'r') as paf_file:
            lines = []
            total_match_and_mismatch_dict = {}

            for line in paf_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip empty lines and comments
                fields = line.split('\t')
                result = process_paf_line(fields)
                if result:
                    lines.append(result)
                    pair = tuple(sorted([result['query_accession_name'], result['target_accession_name']]))
                    total_match_and_mismatch = total_match_and_mismatch_dict.get(pair, 0)
                    total_match_and_mismatch += result['line_match_and_mismatch']
                    total_match_and_mismatch_dict[pair] = total_match_and_mismatch

        # Now process each line and write outputs
        with open("alignment_de.tsv", 'w') as divergence_file, \
             open("alignment_k2p.tsv", 'w') as k2p_file:

            for result in lines:
                pair = tuple(sorted([result['query_accession_name'], result['target_accession_name']]))
                total_match_and_mismatch = total_match_and_mismatch_dict[pair]
                divergence_weight = result['line_match_and_mismatch'] / total_match_and_mismatch
                divergence_weight_str = f"{divergence_weight:.12f}"
                query_name = result['query_name']
                target_name = result['target_name']
                divergence = result['divergence']
                k2p = result['k2p']

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
