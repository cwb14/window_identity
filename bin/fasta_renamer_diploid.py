import re
import os
import subprocess
from argparse import ArgumentParser
from itertools import combinations

def rename_header(header, file_prefix):
#    chrom_regex = re.compile(r'(Chr|chr|Chro|chro|Chrom|chrom|Chromosome|chromosome)[ _-]*(\d+|[XYZW])')
    chrom_regex = re.compile(r'(Chr|chr|Chro|chro|Chrom|chrom|Chromosome|chromosome|CHROMOSOME|CHR|CHRO|CHROM)[ _-]*:? ?(\d+|[XYZW])')
    match = chrom_regex.search(header)
    if match:
        # Convert chromosome number to integer to remove leading zeros
        chrom_number = match.group(2)
        if chrom_number.isdigit():
            chrom_number = str(int(chrom_number))  # Remove leading zeros
        return f'>{file_prefix}_chr{chrom_number}', f'{file_prefix}_chr{chrom_number}'
    return None, None

def process_file(file_path, pass_files):
    file_prefix = os.path.splitext(os.path.basename(file_path))[0]
    with open(file_path, 'r') as fasta_file:
        lines = fasta_file.readlines()

    sequences = {}
    current_header = None
    for line in lines:
        if line.startswith('>'):
            current_header = line.strip()
            sequences[current_header] = []
        elif current_header:
            sequences[current_header].append(line.strip())

    sequence_lengths = {header: sum(len(seq) for seq in seqs) for header, seqs in sequences.items()}

    # Separate sequences into those with chr designations and those without
    chr_sequences = {header: seqs for header, seqs in sequences.items() if rename_header(header, file_prefix)[0]}
    sca_sequences = {header: seqs for header, seqs in sequences.items() if not rename_header(header, file_prefix)[0]}

    # Sort sca sequences by size
    sorted_sca_sequences = sorted(sca_sequences.items(), key=lambda item: -sequence_lengths[item[0]])

    new_lines = []
    used_headers = set()
    fallback_counter = 1
    header_mapping = []

    # Process chr sequences without sorting
    for header, seqs in chr_sequences.items():
        new_header, new_header_full = rename_header(header, file_prefix)
        old_header = header.strip()[1:]
        if new_header and new_header not in used_headers:
            used_headers.add(new_header)
            new_lines.append(f'{new_header}\n')
            header_mapping.append(f'{old_header}\t{new_header.strip(">")}\n')
        new_lines.append(''.join(seqs) + '\n')

    # Process sorted sca sequences with length check
    for header, seqs in sorted_sca_sequences:
        while f'>{file_prefix}_sca{fallback_counter}' in used_headers:
            fallback_counter += 1
        new_header_full = f'>{file_prefix}_sca{fallback_counter}'
        if len(new_header_full) > 14:
            print(f"Excluding sequence with header {new_header_full} due to length > 13 characters.")
            break
        old_header = header.strip()[1:]
        new_lines.append(f'{new_header_full}\n')
        used_headers.add(new_header_full)
        fallback_counter += 1
        header_mapping.append(f'{old_header}\t{new_header_full.strip(">")}\n')
        new_lines.append(''.join(seqs) + '\n')

    # Feed the new_lines output directly to bioawk via stdin
    bioawk_command = ["bioawk", "-c", "fastx", '{print ">"$name; print $seq}']
    process = subprocess.Popen(bioawk_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    fasta_content = ''.join(new_lines).encode('utf-8')
    bioawk_output, _ = process.communicate(fasta_content)

    # Write the bioawk processed output directly to the final file
    output_file_path = os.path.join(".", f"{os.path.splitext(file_path)[0]}_mod.fa")
    with open(output_file_path, 'wb') as output_file:
        output_file.write(bioawk_output)
    print(f"Processed {file_path}, output written to {output_file_path}")

    # Writing the header mapping file
    mapping_file_path = f"{file_path.rsplit('.', 1)[0]}_chrIDs.txt"
    with open(mapping_file_path, 'w') as mapping_file:
        mapping_file.writelines(header_mapping)
    print(f"Header mapping for {file_path} written to {mapping_file_path}")

    # Check if a pass file exists for this file prefix and process it if so
    pass_file = f'{file_prefix}.pass.list'
    if pass_file in pass_files:
        process_pass_file(pass_file, mapping_file_path)

def process_pass_file(pass_file, mapping_file_path):
    mappings = {}
    with open(mapping_file_path, 'r') as mapping_file:
        for line in mapping_file:
            parts = line.strip().rsplit('\t', 1)
            old = parts[0]
            new = parts[1]
            mappings[old] = new

    with open(pass_file, 'r') as pf:
        pass_lines = pf.readlines()

    new_pass_lines = []
    for line in pass_lines:
        columns = line.strip().split('\t')
        old_id_part = columns[0].split(':')[0]
        if old_id_part in mappings:
            columns[0] = columns[0].replace(old_id_part, mappings[old_id_part], 1)
        new_pass_lines.append('\t'.join(columns) + '\n')

    prefix = os.path.splitext(os.path.splitext(pass_file)[0])[0]
    output_pass_file = f"{prefix}_mod.pass.list"
    with open(output_pass_file, 'w') as opf:
        opf.writelines(new_pass_lines)
    print(f"Processed {pass_file}, output written to {output_pass_file}")

def generate_jcvi_list(genome_prefixes):
    with open('jcvi_list.txt', 'w') as jcvi_file:
        for pair in combinations(genome_prefixes, 2):
            jcvi_file.write('\t'.join(pair) + '\n')
    print("Generated jcvi_list.txt with all pairwise relationships.")

def write_prefix_list(filename, prefixes):
    with open(filename, 'w') as file:
        for prefix in prefixes:
            file.write(prefix + '\n')

def main(genomes, pass_files):
    genome_prefixes = [os.path.splitext(os.path.basename(genome))[0] for genome in genomes]
    pass_prefixes = [os.path.splitext(os.path.splitext(os.path.basename(pass_file))[0])[0] for pass_file in pass_files]

    for genome_file in genomes:
        process_file(genome_file, pass_files)
    
    generate_jcvi_list(genome_prefixes)
    write_prefix_list('genome_list.txt', genome_prefixes)
    write_prefix_list('pass_list.txt', pass_prefixes)

if __name__ == "__main__":
    parser = ArgumentParser(description="Rename FASTA headers and update pass files based on chromosome numbers.")
    parser.add_argument("-genomes", nargs='+', help="FASTA files to process.", required=True)
    parser.add_argument("-pass_files", nargs='*', help="Pass files to update based on header mappings.", default=[])
    args = parser.parse_args()

    main(args.genomes, args.pass_files)

# END
