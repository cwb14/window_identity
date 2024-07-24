# Chris Benson
# 07/24/2024

import csv
import argparse
import re
from collections import defaultdict

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Calculate SNP fraction within windows and chromosomes from a VCF file.")
    parser.add_argument('-window_bed', required=True, help="Path to the windows file.")
    parser.add_argument('-vcf', required=True, help="Path to the VCF file.")
    parser.add_argument('-id', required=True, help="ID to be used as the header for SNP fraction column in the output window file.")
    parser.add_argument('-output_window', required=True, help="Path to the output window file.")
    parser.add_argument('-output_report', required=True, help="Path to the output report file.")
    return parser.parse_args()

def read_vcf(vcf_file):
    """
    Reads the VCF file and extracts chromosome, position data, and chromosome sizes from the header.
    
    Parameters:
    vcf_file (str): Path to the VCF file.

    Returns:
    list: A list of tuples containing chromosome and position.
    dict: A dictionary of chromosome sizes.
    """
    vcf_data = []
    chrom_sizes = {}
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            # Extract chromosome sizes from the header.
            if line.startswith('##contig'):
                match = re.match(r'##contig=<ID=([^,]+),length=(\d+)>', line)
                if match:
                    chrom, length = match.groups()
                    chrom_sizes[chrom] = int(length)
            if line.startswith('#'):
                continue
            # Extract chromosome and position for each variant.
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            vcf_data.append((chrom, pos))
    return vcf_data, chrom_sizes

def read_windows(windows_file):
    """
    Reads the windows file and extracts chromosome, start, end, and calculates window size.
    
    Parameters:
    windows_file (str): Path to the windows file.

    Returns:
    list: A list of tuples containing chromosome, start, end, and window size.
    """
    windows_data = []
    with open(windows_file, 'r') as win:
        for line in win:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            win_size = end - start
            windows_data.append((chrom, start, end, win_size))
    return windows_data

def calculate_snp_fraction(vcf_data, windows_data):
    """
    Calculates the SNP count and fraction for each window.
    
    Parameters:
    vcf_data (list): List of tuples containing chromosome and position from VCF.
    windows_data (list): List of tuples containing chromosome, start, end, and window size from windows file.

    Returns:
    list: A list of tuples containing chromosome, start, end, and SNP fraction.
    """
    results = []
    for chrom, start, end, win_size in windows_data:
        # Ensure robust matching of chromosomes.
        snp_count = sum(1 for v_chrom, v_pos in vcf_data if v_chrom == chrom and start <= v_pos < end)
        snp_fract = round(snp_count / win_size, 5)
        results.append((chrom, start, end, snp_fract))
        # Debugging statement to show SNP count per window.
#        print(f"Window {chrom}:{start}-{end} has {snp_count} SNPs")
    return results

def calculate_chromosome_report(vcf_data, chrom_sizes):
    """
    Calculates the SNP fraction for each chromosome and for the whole genome.
    
    Parameters:
    vcf_data (list): List of tuples containing chromosome and position from VCF.
    chrom_sizes (dict): Dictionary containing chromosome sizes.

    Returns:
    list: A list of tuples containing chromosome and SNP fraction.
    """
    chrom_snp_counts = defaultdict(int)
    for chrom, pos in vcf_data:
        chrom_snp_counts[chrom] += 1
    
    report = []
    genome_snp_count = 0
    genome_size = 0
    for chrom, size in chrom_sizes.items():
        snp_count = chrom_snp_counts[chrom]
        snp_fract = round(snp_count / size, 5)
        report.append((chrom, snp_fract))
        genome_snp_count += snp_count
        genome_size += size

    genome_snp_fract = round(genome_snp_count / genome_size, 5)
    report.append(("WholeGenome", genome_snp_fract))
    return report

def write_output_window(output_file, results, id):
    """
    Writes the window results to the output file.
    
    Parameters:
    output_file (str): Path to the output file.
    results (list): List of tuples containing chromosome, start, end, and SNP fraction.
    id (str): The ID to be used in the header for SNP fraction column.
    """
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        # Write header.
        writer.writerow(["chromosome", "start_position", "end_position", id])
        # Write data.
        for row in results:
            writer.writerow(row)
    print(f'Window output written to {output_file}')

def write_output_report(output_file, report):
    """
    Writes the chromosome report to the output file.
    
    Parameters:
    output_file (str): Path to the output file.
    report (list): List of tuples containing chromosome and SNP fraction.
    """
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["chromosome", "avg_snp_fract"])
        for row in report:
            writer.writerow(row)
    print(f'Report output written to {output_file}')

def main():
    """
    Main function to execute the script logic.
    """
    args = parse_arguments()
    
    # Read VCF data and chromosome sizes.
    vcf_data, chrom_sizes = read_vcf(args.vcf)
    
    # Read window data and calculate SNP fractions for windows.
    windows_data = read_windows(args.window_bed)
    window_results = calculate_snp_fraction(vcf_data, windows_data)
    
    # Write window results to output file with the specified ID in the header.
    write_output_window(args.output_window, window_results, args.id)
    
    # Calculate and write chromosome report.
    chromosome_report = calculate_chromosome_report(vcf_data, chrom_sizes)
    write_output_report(args.output_report, chromosome_report)

if __name__ == "__main__":
    main()
