#!/usr/bin/env python3

import os
import sys
import argparse
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import re
import tempfile  # Import tempfile for creating unique temp directories

# Function to parse syncoord and extract sequences
def extract_sequence(fasta_file, syncoord, reverse_complement=False):
    try:
        # Split the syncoord into accessionID_seqID and start..end
        accession_seqid, pos_range = syncoord.split(":")
        start, end = map(int, pos_range.split(".."))
        
        accession_id, seqid = accession_seqid.split("_")
        
        print(f"Parsed syncoord: accession_id={accession_id}, seqid={seqid}, start={start}, end={end}")

    except ValueError as e:
        print(f"Error parsing syncoord: {syncoord}. Expected format 'accessionID_seqID:start..end'. Error: {e}")
        return None

    print(f"Extracting sequence from {fasta_file} with syncoord {syncoord} (reverse_complement={reverse_complement})")

    # Read the genome file
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == f"{accession_id}_{seqid}":
                print(f"Found matching header: {record.id}")
                seq = record.seq[start-1:end]
                if reverse_complement:
                    print(f"Reverse complementing sequence: {record.id}")
                    seq = seq.reverse_complement()
                return SeqRecord(seq, id=record.id, description="")
    print(f"No matching header found for {syncoord} in {fasta_file}")
    return None

# Function to extract start position from file name
def extract_start_pos(file_name):
    # file_name is something like 'TaestD_chr4_100_200.fa'
    base_name = os.path.basename(file_name)
    base_name = base_name.replace('.fa', '')
    parts = base_name.split('_')
    start = int(parts[-2])
    return start

# Function to run minimap2 alignment and adjust coordinates
def run_minimap(seq1_file, seq2_file, output_file, kmer, threads):
    cmd = f"minimap2 -t {threads} --secondary=no -k {kmer} --cs=short -x asm5 -c {seq1_file} {seq2_file}"
    print(f"Running minimap2 with command: {cmd}")
    
    # Extract start positions from seq1_file and seq2_file
    start1 = extract_start_pos(seq1_file)
    start2 = extract_start_pos(seq2_file)
    
    # Run minimap2 and capture the output
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running minimap2: {result.stderr}")
        return
    
    # Process the output line by line
    adjusted_lines = []
    for line in result.stdout.strip().split('\n'):
        if line.strip() == '':
            continue
        cols = line.strip().split('\t')
        # PAF format has at least 12 columns
        if len(cols) < 12:
            print(f"Invalid PAF line: {line}")
            continue
        # Adjust columns 3 and 4 (0-based)
        cols[2] = str(int(cols[2]) + start1)
        cols[3] = str(int(cols[3]) + start1)
        # Adjust columns 8 and 9 (0-based)
        cols[7] = str(int(cols[7]) + start2)
        cols[8] = str(int(cols[8]) + start2)
        # Reconstruct the line
        adjusted_line = '\t'.join(cols)
        adjusted_lines.append(adjusted_line)
    
    # Write adjusted lines to output_file
    with open(output_file, 'a') as f:
        for line in adjusted_lines:
            f.write(line + '\n')

# Function to process each line of coords
def process_line(args):
    line, genome_dir, temp_dir_base, kmer, threads, output_file, debug_mode = args
    syncoord1, syncoord2, strand = line.strip().split("\t")
    
    print(f"Processing line: {line.strip()}")

    # Create a unique temp directory for this process_line call
    temp_dir = tempfile.mkdtemp(dir=temp_dir_base)
    print(f"Created temp directory: {temp_dir}")
    
    # Create temp files for the two syncoords
    seq1_file = os.path.join(temp_dir, syncoord1.replace(":", "_").replace("..", "_") + ".fa")
    seq2_file = os.path.join(temp_dir, syncoord2.replace(":", "_").replace("..", "_") + ".fa")
    
    print(f"Temp files: {seq1_file}, {seq2_file}")

    # Extract sequences
    genome1 = os.path.join(genome_dir, syncoord1.split(":")[0].split("_")[0] + "_mod.fa")
    genome2 = os.path.join(genome_dir, syncoord2.split(":")[0].split("_")[0] + "_mod.fa")
    
    print(f"Accessing genome files: {genome1}, {genome2}")

    seq1 = extract_sequence(genome1, syncoord1)
    seq2 = extract_sequence(genome2, syncoord2, reverse_complement=(strand == "-"))

    if seq1 is None or seq2 is None:
        print(f"Skipping line due to missing sequences for {syncoord1} or {syncoord2}")
        # Clean up temp directory
        if not debug_mode:
            os.rmdir(temp_dir)
        return

    # Write sequences to temp files
    SeqIO.write(seq1, seq1_file, "fasta")
    SeqIO.write(seq2, seq2_file, "fasta")
    print(f"Wrote temp sequence files: {seq1_file}, {seq2_file}")

    # Run minimap2 alignment and adjust coordinates
    run_minimap(seq1_file, seq2_file, output_file, kmer, threads)

    # Clean up temp files and temp directory if not in debug mode
    if not debug_mode:
        print(f"Removing temp files and directory: {temp_dir}")
        os.remove(seq1_file)
        os.remove(seq2_file)
        os.rmdir(temp_dir)
    else:
        print(f"Debug mode on, temp files and directory retained: {temp_dir}")

# Main function
def main():
    parser = argparse.ArgumentParser(description="Syntenic genomic sequence extraction and alignment using minimap2.")
    parser.add_argument("-pairedIDs", nargs='+', help="List of accession IDs to include for analysis.")
    parser.add_argument("-singleID", help="Single accession ID to include for analysis.")
    parser.add_argument("-k", "--kmer", type=int, default=28, help="Kmer size for minimap2.")
    parser.add_argument("-t", "--threads", type=int, default=3, help="Number of threads for minimap2.")
    parser.add_argument("-p", "--processes", type=int, default=10, help="Number of parallel processes.")
    parser.add_argument("-c", "--coords", required=True, help="Input file containing syntenic coordinates.")
    parser.add_argument("-d", "--debug", action="store_true", help="Debug mode: retain temp files.")

    args = parser.parse_args()

    # Check mutual exclusivity between pairedIDs and singleID
    if args.pairedIDs and args.singleID:
        parser.error("Cannot use both -pairedIDs and -singleID. Use one option.")

    # Set output file for minimap2 alignments
    output_file = "alignment.paf"
    print(f"Output file: {output_file}")

    # Load coords
    with open(args.coords, "r") as f:
        lines = f.readlines()

    print(f"Loaded {len(lines)} lines from {args.coords}")

    # Filter lines based on pairedIDs or singleID
    if args.pairedIDs or args.singleID:
        filtered_lines = []
        for line in lines:
            syncoord1, syncoord2, strand = line.strip().split("\t")
            accessions1 = syncoord1.split(":")[0].split("_")[0]  # Extract accessionID before underscore
            accessions2 = syncoord2.split(":")[0].split("_")[0]  # Extract accessionID before underscore

            if args.pairedIDs and (accessions1 in args.pairedIDs and accessions2 in args.pairedIDs):
                print(f"Adding line to filtered_lines (both in pairedIDs): {line.strip()}")
                filtered_lines.append(line)
            elif args.singleID and (accessions1 == args.singleID or accessions2 == args.singleID):
                print(f"Adding line to filtered_lines (one in singleID): {line.strip()}")
                filtered_lines.append(line)
    else:
        filtered_lines = lines

    print(f"Filtered down to {len(filtered_lines)} lines based on accessions of interest")

    # Prepare temp directory base
    temp_dir_base = "./minimap_temp"
    os.makedirs(temp_dir_base, exist_ok=True)

    # Process lines with multiprocessing
    args_list = [(line, ".", temp_dir_base, args.kmer, args.threads, output_file, args.debug) for line in filtered_lines]
    print(f"Processing {len(args_list)} lines with {args.processes} processes")
    with Pool(processes=args.processes) as pool:
        pool.map(process_line, args_list)

    # Clean up base temp directory if not in debug mode
    if not args.debug:
        print(f"Removing base temp directory: {temp_dir_base}")
        try:
            os.rmdir(temp_dir_base)
        except OSError:
            print(f"Could not remove {temp_dir_base}, it may not be empty.")

if __name__ == "__main__":
    main()
