import pandas as pd
import argparse
import sys

def merge_bed_files(bed_files):
    """
    Merge multiple BED files based on columns 1-3 (chromID, start, end).
    This script is not necissary for 1v1 comparisons since there are no additional alignments to merge.

    Parameters:
    bed_files (list): List of BED file paths to merge.

    Returns:
    pandas.DataFrame: Merged DataFrame with combined column 4 values.
    """
    # Create an empty dictionary to store data frames.
    data_frames = {}
    
    # Read each bed file into a data frame.
    for bed_file in bed_files:
        df = pd.read_csv(bed_file, sep='\t', header=0)
        # Use columns 1-3 as the index for merging.
        df.set_index(['chromosome', 'start_position', 'end_position'], inplace=True)
        # Add the dataframe to the dictionary with the filename as key.
        data_frames[bed_file] = df
    
    # Merge all data frames on their indices (columns 1-3).
    merged_df = pd.concat(data_frames.values(), axis=1)
    
    # Reset index to turn the index back into columns.
    merged_df.reset_index(inplace=True)
    
    return merged_df

def main():
    """
    Main function to parse arguments and merge BED files.
    """
    # Set up argument parser.
    parser = argparse.ArgumentParser(
        description='Merge multiple BED files by matching columns 1-3 and combining column 4 values.'
    )
    parser.add_argument(
        'bed_files', 
        metavar='BED_FILE', 
        type=str, 
        nargs='+', 
        help='List of BED files to merge'
    )
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        default='merged_bed.bed', 
        help='Output file name (default: merged_bed.bed)'
    )
    
    # Parse arguments.
    args = parser.parse_args()
    
    # Check if at least two files are provided.
    if len(args.bed_files) < 2:
        print("Error: Please provide at least two BED files to merge.", file=sys.stderr)
        sys.exit(1)
    
    # Merge the bed files.
    merged_df = merge_bed_files(args.bed_files)
    
    # Write the merged data frame to the specified output file.
    merged_df.to_csv(args.output, sep='\t', index=False)
    print(f"Merged BED file saved to {args.output}")

if __name__ == "__main__":
    main()
