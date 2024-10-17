import sys
import subprocess
import argparse

def parse_fai(fai_file):
    lengths = {}
    total_length = 0
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            seq_name = fields[0]
            seq_length = int(fields[1])
            lengths[seq_name] = seq_length
            total_length += seq_length  # Summing up the total length of the sequences
    return lengths, total_length

def ensure_fai(fasta_file):
    fai_file = fasta_file + ".fai"
    try:
        with open(fai_file, 'r'):
            pass  # FAI file exists
    except FileNotFoundError:
        print(f"FAI index not found for {fasta_file}. Generating using samtools...")
        subprocess.run(["samtools", "faidx", fasta_file])
    return fai_file

def infer_fasta_name(seq_name):
    # Infer fasta filename from the sequence name prefix
    fasta_name = seq_name.split('_')[0] + "_mod.fa"
    return fasta_name

def process_paf(paf_file):
    fai_cache = {}
    total_lengths_cache = {}

    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_name = fields[0]
            target_name = fields[5]

            # Infer fasta file names
            query_fasta = infer_fasta_name(query_name)
            target_fasta = infer_fasta_name(target_name)

            # Ensure FAI files are available and cache the parsed lengths and total genome lengths
            if query_fasta not in fai_cache:
                query_fai = ensure_fai(query_fasta)
                fai_cache[query_fasta], total_lengths_cache[query_fasta] = parse_fai(query_fai)

            if target_fasta not in fai_cache:
                target_fai = ensure_fai(target_fasta)
                fai_cache[target_fasta], total_lengths_cache[target_fasta] = parse_fai(target_fai)

            # Update query length (column 2) and target length (column 7)
            if query_name in fai_cache[query_fasta]:
                fields[1] = str(fai_cache[query_fasta][query_name])

            if target_name in fai_cache[target_fasta]:
                fields[6] = str(fai_cache[target_fasta][target_name])

            # Add query genome length and target genome length to the output
            query_genome_length = total_lengths_cache[query_fasta]
            target_genome_length = total_lengths_cache[target_fasta]
            fields.append(str(query_genome_length))
            fields.append(str(target_genome_length))

            # Output the modified PAF line with the new columns
            print('\t'.join(fields))

def main():
    parser = argparse.ArgumentParser(description="Adjust PAF file based on inferred FAI lengths.")
    parser.add_argument('-paf', required=True, help="Input PAF file")
    
    args = parser.parse_args()

    # Process the PAF file
    process_paf(args.paf)

if __name__ == "__main__":
    main()
    
# END
