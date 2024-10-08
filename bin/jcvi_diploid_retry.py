import argparse
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run jcvi.compara.catalog ortholog in parallel.')
parser.add_argument('-p', metavar='N', type=int, default=16, help='Number of processes to run in parallel (default: 16)')
parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use for each subprocess (default: 1)')
parser.add_argument('--blast', '-b', action='store_true', help='Include --align_soft blast in the command')
args = parser.parse_args()

# Ensure that both --cpus and -p are set to at least 1
args.cpus = max(1, args.cpus)
args.p = max(1, args.p)

# Define the paths to your list files
default_list_file_path = 'jcvi_list.txt'
missing_list_file_path = 'jcvi_list_missing.txt'

# Determine which list file to use
list_file_path = missing_list_file_path if os.path.exists(missing_list_file_path) else default_list_file_path

# Function to check if a file is missing or empty
def is_missing_or_empty(filename):
    return not os.path.exists(filename) or os.path.getsize(filename) == 0

# Function to run a command for a given line
def run_command(line):
    col1, col2 = line.strip().split()
    last_file = f"{col1}.{col2}.last"
    anchors_file = f"{col1}.{col2}.anchors"
    
    # Check if either 'last' or 'anchors' file is missing or empty
    if is_missing_or_empty(last_file) or is_missing_or_empty(anchors_file):
        command = [
            'python', '-m', 'jcvi.compara.catalog', 'ortholog',
            col1, col2, '--no_strip_names', f'--cpus={args.cpus}', '--notex', '--cscore=.99'
        ]
        if args.blast:
            command.extend(['--align_soft', 'blast'])
        subprocess.run(command)

# Read all lines from the selected list file into a list
with open(list_file_path, 'r') as file:
    lines = file.readlines()

# Use ThreadPoolExecutor to run commands in parallel
with ThreadPoolExecutor(max_workers=args.p) as executor:
    executor.map(run_command, lines)

# Once all lines have been processed, run the final script if needed
# subprocess.run(['python', 'jcvi_diploid.py'])
# END
