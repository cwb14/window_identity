#!/usr/bin/env python3

import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert PAF to BED-like files based on a reference ID.")
    parser.add_argument('-ref', required=True, help='Reference ID to filter the PAF file.')
    parser.add_argument('-paf', required=True, help='Path to the input PAF file.')
    return parser.parse_args()

def extract_nonid(column):
    """
    Extracts the nonID part from a column by taking the substring before the first underscore.
    """
    return column.split('_')[0] if '_' in column else column

def process_paf(ref_id, paf_file):
    """
    Processes the PAF file and writes the extracted information to separate output files.
    """
    output_dict = {}

    with open(paf_file, 'r') as infile:
        for line_number, line in enumerate(infile, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # Skip empty lines and comments

            columns = line.split('\t')
            if len(columns) < 21:
                print(f"Warning: Line {line_number} has fewer than 21 columns. Skipping.")
                continue  # Skip lines that don't have enough columns

            # Check if ref_id is in column1 or column6
            if ref_id in columns[0]:
                # ref_id is in column1
                try:
                    extracted = [columns[0], columns[2], columns[3], columns[20]]
                    nonid = extract_nonid(columns[5])
                except IndexError:
                    print(f"Warning: Line {line_number} is malformed. Skipping.")
                    continue
            elif ref_id in columns[5]:
                # ref_id is in column6
                try:
                    extracted = [columns[5], columns[7], columns[8], columns[20]]
                    nonid = extract_nonid(columns[0])
                except IndexError:
                    print(f"Warning: Line {line_number} is malformed. Skipping.")
                    continue
            else:
                # ref_id not present in this line
                continue

            # Initialize the list for this nonid if not already present
            if nonid not in output_dict:
                output_dict[nonid] = []

            # Join the extracted columns with tabs and add to the list
            output_dict[nonid].append('\t'.join(extracted))

    # Write the output files
    for nonid, lines in output_dict.items():
        output_filename = f"{ref_id}.{nonid}.bed"
        with open(output_filename, 'w') as outfile:
            for extracted_line in lines:
                outfile.write(extracted_line + '\n')
        print(f"Written {len(lines)} lines to {output_filename}")

def main():
    args = parse_arguments()
    ref_id = args.ref
    paf_file = args.paf

    if not os.path.isfile(paf_file):
        print(f"Error: The file '{paf_file}' does not exist.")
        return

    process_paf(ref_id, paf_file)

if __name__ == "__main__":
    main()
