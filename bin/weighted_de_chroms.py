import pandas as pd
import argparse

def parse_fai(fai_file):
    """
    Parses the FAI file to extract chromosome sizes.
    
    Args:
    fai_file (str): Path to the .fai file.

    Returns:
    pd.DataFrame: DataFrame containing chromosome names and sizes.
    """
    # Read the .fai file
    fai_df = pd.read_csv(fai_file, sep='\t', header=None, names=['chromosome', 'chr_size', '2', '3', '4'])
    fai_df = fai_df[['chromosome', 'chr_size']]
    return fai_df

def parse_bed(bed_file):
    """
    Parses the BED file to extract alignment data and sequence identity values.
    
    Args:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame containing chromosome names, alignment starts and ends, and sequence identity values.
    """
    # Read the .bed file
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end', 'seq_ident'])
    # Extract numeric part of seq_ident
    bed_df['seq_ident'] = bed_df['seq_ident'].str.extract(r'(\d+\.\d+)$').astype(float)
    # Calculate alignment length
    bed_df['aln_len'] = (bed_df['end'] - bed_df['start']).abs()
    return bed_df

def compute_weighted_average(fai_df, bed_df):
    """
    Computes the weighted average sequence identity for each chromosome and for the whole genome.
    
    Args:
    fai_df (pd.DataFrame): DataFrame containing chromosome sizes.
    bed_df (pd.DataFrame): DataFrame containing alignment data and sequence identity values.

    Returns:
    (pd.DataFrame, float): DataFrame containing weighted average sequence identity for each chromosome, 
                           and weighted average sequence identity for the whole genome.
    """
    # Merge the fai and bed dataframes on chromosome
    merged_df = pd.merge(bed_df, fai_df, on='chromosome')
    # Calculate the weight for each alignment
    merged_df['weight'] = merged_df['aln_len'] / merged_df['chr_size']
    # Calculate the weighted seq_ident
    merged_df['weighted_seq_ident'] = merged_df['weight'] * merged_df['seq_ident']

    # Calculate weighted average seq_ident for each chromosome
    chrom_weighted_avg = merged_df.groupby('chromosome').apply(
        lambda x: x['weighted_seq_ident'].sum() / x['weight'].sum()
    ).reset_index(name='avg_weighted_de')

    # Calculate whole genome weighted average seq_ident
    genome_weighted_avg = merged_df['weighted_seq_ident'].sum() / merged_df['weight'].sum()

    return chrom_weighted_avg, genome_weighted_avg

def main():
    parser = argparse.ArgumentParser(
        description='Calculate weighted average sequence identity (seq_ident) for each chromosome and the whole genome.'
    )
    parser.add_argument('-bed', '--bedfile', required=True, help='Input BED file with alignment data')
    parser.add_argument('-fai', '--faifile', required=True, help='Input FAI file with chromosome sizes')
    args = parser.parse_args()

    # Parse the input files
    fai_df = parse_fai(args.faifile)
    bed_df = parse_bed(args.bedfile)

    # Compute weighted averages
    chrom_weighted_avg, genome_weighted_avg = compute_weighted_average(fai_df, bed_df)

    # Print the results with headers
    print("chromosome\tavg_weighted_de")
    for index, row in chrom_weighted_avg.iterrows():
        print(f"{row['chromosome']}\t{row['avg_weighted_de']:.4f}")
    print(f"WholeGenome\t{genome_weighted_avg:.4f}")

if __name__ == '__main__':
    main()
