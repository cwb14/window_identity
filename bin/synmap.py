#!/usr/bin/env python3

import os
import sys
import argparse
import shutil
import subprocess
import tempfile
import threading
import time
from multiprocessing import Pool, Manager
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Global variables
MINIMAP2_BIN = None
PRESET = "asm10"
VERBOSE = False

# Logging helper
def log(msg):
    if VERBOSE:
        print(msg)

# Function to parse syncoord and extract sequences
def extract_sequence(fasta_file, syncoord, reverse_complement=False):
    try:
        accession_seqid, pos_range = syncoord.split(":")
        start, end = map(int, pos_range.split(".."))
        accession_id, seqid = accession_seqid.split("_")
        log(f"Parsed syncoord: accession_id={accession_id}, seqid={seqid}, start={start}, end={end}")
    except ValueError as e:
        log(f"Error parsing syncoord: {syncoord}. Expected format 'accessionID_seqID:start..end'. Error: {e}")
        return None

    log(f"Extracting sequence from {fasta_file} with syncoord {syncoord} (reverse_complement={reverse_complement})")

    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == f"{accession_id}_{seqid}":
                seq = record.seq[start-1:end]
                if reverse_complement:
                    seq = seq.reverse_complement()
                return SeqRecord(seq, id=record.id, description="")
    log(f"No matching header found for {syncoord} in {fasta_file}")
    return None

# Function to extract start position from file name
def extract_start_pos(file_name):
    base = os.path.basename(file_name).replace('.fa', '')
    parts = base.split('_')
    return int(parts[-2])

# Convert timer string to seconds
def convert_timer_to_seconds(timer_str):
    if not timer_str:
        return None
    unit = timer_str[-1]
    value = int(timer_str[:-1])
    if unit == 'd': return value * 86400
    if unit == 'h': return value * 3600
    if unit == 'm': return value * 60
    if unit == 's': return value
    raise ValueError(f"Invalid time unit in '{timer_str}'. Use 'd', 'h', 'm', or 's'.")

# Run minimap2 alignment and adjust coordinates
def run_minimap(seq1_file, seq2_file, output_file, kmer, threads, timer):
    cmd = (
        f"{MINIMAP2_BIN} -t {threads} --secondary=no "
        f"-k {kmer} --cs=short -x {PRESET} -c {seq1_file} {seq2_file}"
    )
    log(f"Running minimap2 with command: {cmd}")
    start1 = extract_start_pos(seq1_file)
    start2 = extract_start_pos(seq2_file)
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timer)
        if result.returncode != 0:
            log(f"Error running minimap2: {result.stderr}")
            return
        adjusted = []
        for line in result.stdout.strip().split('\n'):
            if not line: continue
            cols = line.split('\t')
            if len(cols) < 12:
                log(f"Invalid PAF line: {line}")
                continue
            cols[2] = str(int(cols[2]) + start2)
            cols[3] = str(int(cols[3]) + start2)
            cols[7] = str(int(cols[7]) + start1)
            cols[8] = str(int(cols[8]) + start1)
            adjusted.append('\t'.join(cols))
        with open(output_file, 'a') as f:
            for l in adjusted:
                f.write(l + '\n')
    except subprocess.TimeoutExpired:
        log(f"Process took longer than {timer} seconds. Skipping this alignment.")

# Process each line of coords
def process_line(args):
    line, genome_dir, temp_base, kmer, threads, output_file, debug_mode, timer, counter = args
    syncoord1, syncoord2, strand = line.strip().split("\t")
    log(f"Processing line: {line.strip()}")

    temp_dir = tempfile.mkdtemp(dir=temp_base)
    seq1_file = os.path.join(temp_dir, syncoord1.replace(":", "_").replace("..", "_") + ".fa")
    seq2_file = os.path.join(temp_dir, syncoord2.replace(":", "_").replace("..", "_") + ".fa")

    genome1 = os.path.join(genome_dir, syncoord1.split(":")[0].split("_")[0] + "_mod.fa")
    genome2 = os.path.join(genome_dir, syncoord2.split(":")[0].split("_")[0] + "_mod.fa")

    seq1 = extract_sequence(genome1, syncoord1)
    seq2 = extract_sequence(genome2, syncoord2, reverse_complement=(strand == "-"))

    if seq1 and seq2:
        SeqIO.write(seq1, seq1_file, "fasta")
        SeqIO.write(seq2, seq2_file, "fasta")
        run_minimap(seq1_file, seq2_file, output_file, kmer, threads, timer)
    # Clean up
    if not debug_mode:
        try:
            os.remove(seq1_file)
            os.remove(seq2_file)
            os.rmdir(temp_dir)
        except OSError:
            pass
    # Increment processed counter
    counter.value += 1

# Main entry point
def main():
    parser = argparse.ArgumentParser(
        description="Syntenic genomic sequence extraction and alignment using minimap2."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-pairedIDs", nargs='+', help="List of accession IDs to include for analysis.")
    group.add_argument("-singleID", help="Single accession ID to include for analysis.")
    parser.add_argument("-k", "--kmer", type=int, default=28,
                        help="Kmer size for minimap2.")
    parser.add_argument("-t", "--threads", type=int, default=3,
                        help="Number of threads for minimap2.")
    parser.add_argument("-p", "--processes", type=int, default=10,
                        help="Number of parallel processes.")
    parser.add_argument("-c", "--coords", required=True,
                        help="Input file containing syntenic coordinates.")
    parser.add_argument("--timer", type=str,
                        help="Set a timer for each minimap2 process. Format: [int][d/h/m/s].")
    parser.add_argument("--preset", choices=["asm5", "asm10", "asm20"], default="asm10",
                        help="minimap2 preset (default: asm10).")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose QC output.")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="Debug mode: retain temp files.")
    args = parser.parse_args()

    global PRESET, VERBOSE
    PRESET = args.preset
    VERBOSE = args.verbose

    output_file = "alignment.paf"
    log(f"Output file: {output_file}")

    # Download latest minimap2
    temp_minimap_dir = tempfile.mkdtemp(prefix="minimap2_download_")
    metadata = subprocess.check_output([
        "curl", "-s",
        "https://api.github.com/repos/lh3/minimap2/releases/latest"
    ], text=True)
    import json
    tag = json.loads(metadata)["tag_name"].lstrip("v")
    tarball = f"minimap2-{tag}_x64-linux.tar.bz2"
    url = f"https://github.com/lh3/minimap2/releases/download/v{tag}/{tarball}"
    log(f"Downloading minimap2 v{tag} from {url}")
    subprocess.run(
        f"curl -L {url} | tar -jxvf - -C {temp_minimap_dir}",
        shell=True, check=True
    )
    global MINIMAP2_BIN
    MINIMAP2_BIN = os.path.join(
        temp_minimap_dir,
        f"minimap2-{tag}_x64-linux",
        "minimap2"
    )
    log(f"Using minimap2 binary at: {MINIMAP2_BIN}")

    # Read and filter coords
    with open(args.coords, "r") as f:
        lines = f.readlines()
    filtered = []
    for line in lines:
        syn1, syn2, _ = line.strip().split("\t")
        acc1 = syn1.split(":")[0].split("_")[0]
        acc2 = syn2.split(":")[0].split("_")[0]
        if args.pairedIDs and (acc1 in args.pairedIDs and acc2 in args.pairedIDs):
            filtered.append(line)
        elif args.singleID and (acc1 == args.singleID or acc2 == args.singleID):
            filtered.append(line)
        elif not args.pairedIDs and not args.singleID:
            filtered = lines
    total = len(filtered)
    log(f"Total lines to process: {total}")

    # Prepare temp base
    temp_base = os.path.abspath("./minimap_temp")
    os.makedirs(temp_base, exist_ok=True)

    # Shared counter for progress
    manager = Manager()
    counter = manager.Value('i', 0)

    # Progress reporter thread
    stop_event = threading.Event()
    def reporter():
        iteration = 1
        while not stop_event.is_set() and counter.value < total:
            time.sleep(60)
            if counter.value >= total:
                break
            percent = int((counter.value / total) * 100)
            print(f"{iteration*1} minute: {percent}% complete")
            iteration += 1
    thread = threading.Thread(target=reporter)
    thread.daemon = True
    thread.start()

    # Launch pool of workers
    timer_sec = convert_timer_to_seconds(args.timer)
    tasks = [
        (line, ".", temp_base, args.kmer, args.threads,
         output_file, args.debug, timer_sec, counter)
        for line in filtered
    ]
    with Pool(processes=args.processes) as pool:
        pool.map(process_line, tasks)

    # Signal reporter to stop
    stop_event.set()

    # Cleanup
    if not args.debug:
        try:
            os.rmdir(temp_base)
        except OSError:
            pass
    log(f"Removing downloaded minimap2 directory: {temp_minimap_dir}")
    shutil.rmtree(temp_minimap_dir)

if __name__ == "__main__":
    main()
