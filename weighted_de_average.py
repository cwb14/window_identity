import argparse
import logging

# Configure logging.
logging.basicConfig(level=logging.WARNING, format='%(message)s')

def compute_weighted_average(line):
    # Split the input line into parts separated by whitespace
    parts = line.strip().split()
    category = parts[:3]  # The first three parts are the category (chromosome, start, end).
    measurements = parts[3:]  # The remaining parts are the measurements and their weights.

    weights = []
    values = []

    # Iterate through the measurements and weights.
    for i in range(0, len(measurements), 2):
        try:
            weight = float(measurements[i])  # Extract the weight.
            value = float(measurements[i+1].split(':')[-1])  # Extract the value (ignoring preceding characters).
            weights.append(weight)
            values.append(value)
        except (IndexError, ValueError) as e:
            logging.warning(f"Skipping invalid seq identity measurments: {measurements[i:i+2]}")
            logging.warning(f"Error processing line: {line}")
            continue

    if not weights or not values:
        logging.warning(f"No valid alignments found in line: {line}. Setting weighted identity to NA.")
        return f"{category[0]}\t{category[1]}\t{category[2]}\tNA"

    # Calculate the sum of weighted values and sum of weights.
    sum_weighted_values = sum(w * v for w, v in zip(weights, values))
    sum_weights = sum(weights)

    # Compute the weighted average.
    weighted_average = sum_weighted_values / sum_weights if sum_weights != 0 else 0

    # Return the result in the required format.
    return f"{category[0]}\t{category[1]}\t{category[2]}\t{weighted_average}"

def process_file(input_file, id):
    # Print the header line.
    print(f"chromosome\tstart_position\tend_position\t{id}")
    
    # Open the input file for reading.
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.strip():  # Ensure the line is not empty.
                result = compute_weighted_average(line)
                if result:  # Only print valid results.
                    print(result)

def main():
    parser = argparse.ArgumentParser(description='Compute a composite gap compressed sequence identiy measurement using weighted averages in windows.')
    parser.add_argument('input_file', help='The input file containing identity measurements and weights. The input file can be created with weighted_de_scores.py')
    parser.add_argument('id', help='The ID for the fourth column header.')
    
    # Parse the command-line arguments.
    args = parser.parse_args()
    
    # Process the input file.
    process_file(args.input_file, args.id)

if __name__ == '__main__':
    main()

"""
Example input file (weighted_de_scores.txt):
NIPT2T_chr1     0       1000000 0.000347        de:f:0.0086     0.000434        de:f:0.0161     0.000485        de:f:0.0062     0.000566        de:f:0.0229
NIPT2T_chr9     2500000 3500000

Columns 1-3 are chromosomes and window coordinates.
Additional columns are weights and cooresponding seq identity. 
For example, '0.000347' is the weight associated with the identiy measurment, 'de:f:0.0086'. 
'0.000434' is the weight associated with identity measrument, 'de:f:0.0161'. etc.

Example command to run the script:
python weighted_average.py weighted_de_scores.txt Y476h2 > weighted_average.bed

Expected output (weighted_average.bed):
chromosome      start_position  end_position    Y476h2
NIPT2T_chr1     0       1000000 0.000563
NIPT2T_chr9     2500000 3500000 1.0
"""
