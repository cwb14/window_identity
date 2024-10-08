import sys

def process_file(infile):
    """
    Process the anchors coordinate file to merge lines based on overlapping sequence ranges.
    The goal is to reduce I/O constraint during the minimap step.

    Parameters:
    infile (str): Path to the coordinate file
    """
    bins = {}
    merged_lines = []

    # Read the input file and process each line
    with open(infile, 'r') as f:
        for line in f:
            # Split the line into components
            parts = line.strip().split('\t')
            pair1 = parts[0].split(':')
            pair1_id = pair1[0]
            pair1_start, pair1_end = map(int, pair1[1].split('..'))

            pair2 = parts[1].split(':')
            pair2_id = pair2[0]
            pair2_start, pair2_end = map(int, pair2[1].split('..'))

            strand = parts[2]

            # Calculate lengths
            len1 = pair1_end - pair1_start
            len2 = pair2_end - pair2_start

            # Determine threshold status
            threshold_status = "pass" if len1 >= 1000000 and len2 >= 1000000 else "fail"

            # Create a bin key based on pair IDs and strand
            bin_key = f"{pair1_id}_{pair2_id}_{'plus' if strand == '+' else 'minus'}"
            new_line = [pair1_id, pair1_start, pair1_end, pair2_id, pair2_start, pair2_end, strand, len1, len2, bin_key, threshold_status]

            # Add the line to the appropriate bin
            if bin_key not in bins:
                bins[bin_key] = []
            bins[bin_key].append(new_line)

    # Merge lines within each bin
    for bin_key, lines in bins.items():
        while lines:
            current_line = lines.pop(0)
            merged = False

            # Check for overlaps and merge lines if possible
            for i, line in enumerate(lines):
                if (current_line[1] <= line[2] and current_line[2] >= line[1]) and (current_line[4] <= line[5] and current_line[5] >= line[4]):
                    if current_line[10] == "fail" or line[10] == "fail":
                        merged_line = [
                            current_line[0],
                            min(current_line[1], line[1]),
                            max(current_line[2], line[2]),
                            current_line[3],
                            min(current_line[4], line[4]),
                            max(current_line[5], line[5]),
                            current_line[6],
                            min(current_line[7], line[7]),
                            max(current_line[8], line[8]),
                            current_line[9],
                            "pass"  # Set the merged threshold_status to "pass"
                        ]
                        lines[i] = merged_line
                        merged = True
                        break

            if not merged:
                merged_lines.append(current_line)

    # Output the merged lines in the original format
    for line in merged_lines:
        print(f"{line[0]}:{line[1]}..{line[2]}\t{line[3]}:{line[4]}..{line[5]}\t{line[6]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python anchor_coord_consolidator.py all.recip.anchors.coords > outfile\n")
        print("This script processes the input file to merge lines based on overlapping coordinate ranges")
        print("Parameters:")
        print("infile : Path to the input file")
        print("outfile: Path to the output file (use redirection to specify the output file)\n")
        print("Example:")
        print("python script.py infile > outfile")
        sys.exit(1)

    infile = sys.argv[1]
    process_file(infile)
