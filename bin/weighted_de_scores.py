# Chris Benson
# July 11th, 2024.

import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_bed_file(filepath, has_score=False):
    """
    Read a BED file into a pandas DataFrame.
    
    Parameters:
    filepath (str): The path to the BED file.
    has_score (bool): Whether the BED file includes a score column.
    
    Returns:
    pd.DataFrame: The BED data.
    """
    if has_score:
        with open(filepath, 'r') as file:
            lines = file.readlines()
            data = []
            for line in lines:
                parts = line.split()
                if len(parts) == 4:
                    chrom, start, end, de_score = parts
                    data.append([chrom, int(start), int(end), de_score])
                elif len(parts) == 6:
                    chrom, start, end, de_score = parts[0], parts[2], parts[3], parts[4]
                    data.append([chrom, int(start), int(end), de_score])
        return pd.DataFrame(data, columns=['chrom', 'start', 'end', 'de_score'])
    else:
        return pd.read_csv(filepath, sep='\t', header=None, names=['chrom', 'start', 'end'])

def calculate_overlap(start1, end1, start2, end2):
    """
    Calculate the overlap between two intervals.
    
    Parameters:
    start1 (int): Start of the first interval.
    end1 (int): End of the first interval.
    start2 (int): Start of the second interval.
    end2 (int): End of the second interval.
    
    Returns:
    int: The size of the overlap.
    """
    return max(0, min(end1, end2) - max(start1, start2))

def process_chunk(chunk, bed2):
    results = []
    for idx1, row1 in chunk.iterrows():
        chrom1, start1, end1 = row1['chrom'], row1['start'], row1['end']
        output_row = [chrom1, start1, end1]
        
        overlaps = []
        
        for idx2, row2 in bed2.iterrows():
            chrom2, start2, end2, de_score = row2['chrom'], row2['start'], row2['end'], row2['de_score']
            if chrom1 == chrom2:
                overlap_size = calculate_overlap(start1, end1, start2, end2)
                if overlap_size > 0:
                    window_size = end1 - start1
                    overlap_fraction = overlap_size / window_size
                    overlaps.append((overlap_fraction, de_score))
        
        overlaps.sort(key=lambda x: (x[0], x[1]))
        
        for overlap_fraction, de_score in overlaps:
            output_row.append(overlap_fraction)
            output_row.append(de_score)
        
        results.append(output_row)
    return results

def main():
    # Command-line argument parsing.
    parser = argparse.ArgumentParser(description='Calculate overlaps between two BED files.')
    parser.add_argument('-window_bed', required=True, help='Path to the BED file with windows (e.g., windows.bed).')
    parser.add_argument('-minimap_bed', required=True, help='Path to the BED file with minimap results (e.g., minimap.bed).')
    parser.add_argument('-output', required=True, help='Path to the output file.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads to use for processing.')
    args = parser.parse_args()

    bed1_path = args.window_bed
    bed2_path = args.minimap_bed

    bed1 = read_bed_file(bed1_path)
    bed2 = read_bed_file(bed2_path, has_score=True)

    chunk_size = len(bed1) // args.threads
    chunks = [bed1.iloc[i:i + chunk_size] for i in range(0, len(bed1), chunk_size)]

    results = []

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_chunk, chunk, bed2): chunk for chunk in chunks}
        for future in as_completed(futures):
            results.extend(future.result())

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()
