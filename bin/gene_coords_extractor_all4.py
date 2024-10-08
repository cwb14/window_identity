import argparse
import os
import sys

def modify_bed_file(file_name):
    modified_lines = []
    with open(file_name, 'r') as f:
        for line in f:
            modified_lines.append(line.strip())
    return modified_lines

def read_bed_data(bed_files):
    merged_data = []
    for bed_file in bed_files:
        bed_data = modify_bed_file(bed_file)
        merged_data.extend(bed_data)
    return merged_data

def get_gene_coords(merged_data):
    gene_coords = {}
    for line in merged_data:
        chrom, start, end, geneID = line.split("\t")
        gene_coords[geneID] = (chrom, int(start), int(end))
    return gene_coords

def read_gene_ids_file(gene_ids_filename):
    clusters = []
    with open(gene_ids_filename, 'r') as f:
        cluster = []
        for line in f:
            if line.strip() == "###":
                if cluster:
                    clusters.append(cluster)
                    cluster = []
            else:
                genes = line.strip().split()
                cluster.append(genes)
    return clusters

def process_clusters(clusters, gene_coords):
    output_strings = []
    for cluster in clusters:
        for i in range(len(cluster) - 1):  
            coords_1 = gene_coords[cluster[i][0]]
            coords_2 = gene_coords[cluster[i][1]]
            coords_3 = gene_coords[cluster[i + 1][0]]
            coords_4 = gene_coords[cluster[i + 1][1]]

            direction_1 = "+" if coords_1[1] < coords_3[1] else "-"
            direction_2 = "+" if coords_2[1] < coords_4[1] else "-"
            directionality = "+" if direction_1 == direction_2 else "-"

            output_str = f"{coords_1[0]}:{min(coords_1[1], coords_3[1])}..{max(coords_1[2], coords_3[2])}\t"
            output_str += f"{coords_2[0]}:{min(coords_2[1], coords_4[1])}..{max(coords_2[2], coords_4[2])}\t"
            output_str += directionality
            output_strings.append(output_str)
    return output_strings

def derive_bed_files(mcscan_filename):
    name_parts = mcscan_filename.split('.')[0:2]
    bed_files = [f"{name_parts[0]}.bed", f"{name_parts[1]}.bed"]
    return bed_files

def main():
    parser = argparse.ArgumentParser(description='Extract gene coordinates and determine directionality.')
    parser.add_argument('-mcscan', required=True, help='Input file containing geneIDs.')
    args = parser.parse_args()
    
    bed_files = derive_bed_files(args.mcscan)
    merged_data = read_bed_data(bed_files)
    gene_coords = get_gene_coords(merged_data)
    clusters = read_gene_ids_file(args.mcscan)
    output_strings = process_clusters(clusters, gene_coords)
    
    for output_str in output_strings:
        print(output_str)

if __name__ == "__main__":
    main()
