import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run jcvi.compara.catalog ortholog in parallel.')
parser.add_argument('-p', metavar='N', type=int, default=1, help='Number of processes to run in parallel (default: 1)')
parser.add_argument('--cpus', metavar='N', type=int, default=1, help='Number of CPUs to use (default: 1)')
parser.add_argument('--blast', action='store_true', help='Include --align_soft blast if specified')
args = parser.parse_args()

# Ensure that both --cpus and -p are set to at least 1
args.cpus = max(1, args.cpus)
args.p = max(1, args.p)

# Define the paths to your list files
list_file_path = 'jcvi_list.txt'
missing_list_file_path = 'jcvi_list_missing.txt'

# Determine which file to use
file_to_use = missing_list_file_path if os.path.isfile(missing_list_file_path) else list_file_path

# Function to run a command for a given line
def run_command(line):
    col1, col2 = line.strip().split()
    command = [
        'python', '-m', 'jcvi.compara.catalog', 'ortholog',
        col1, col2, '--no_strip_names', '--cpus={}'.format(args.cpus),
        '--notex', '--cscore=.99'
    ]
    if args.blast:
        command.extend(['--align_soft', 'blast'])
    subprocess.run(command)

# Read all lines from the file into a list
with open(file_to_use, 'r') as file:
    lines = file.readlines()

# Use ThreadPoolExecutor to run commands in parallel
with ThreadPoolExecutor(max_workers=args.p) as executor:
    executor.map(run_command, lines)

# Once all lines have been processed, run the final script
# subprocess.run(['python', 'jcvi_diploid.py'])

# [END]
