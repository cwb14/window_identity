#!/usr/bin/env python3
"""Summarise per-gene Ks into one distance per genome pair.

Reads the merged KaKs_Calculator tables ('{ID1}.{ID2}.kaks.tsv') for every pair listed in
jcvi_list.txt and emits:

  -genome_out   ID1  ID2  n_pairs  median_Ks   (the 4-column layout matrix_builder.py eats)
  -long_out     pair  id1  id2  Ks             (one row per gene pair, for the density plot)

The genome-pair distance is the MEDIAN Ks over gene pairs. Ks has a long right tail from
paralogous and saturated anchors, which drags the mean upward; the median does not move.

Sites with Ks outside (0, max_ks) are dropped: Ks == 0 means the two CDS are identical at
synonymous sites (no information, and it would pull a distance toward zero), and Ks above
~2 is saturated, where the estimator is unstable and effectively unbounded.
"""

import argparse
import csv
import statistics
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarise per-gene Ks into one median distance per genome pair."
    )
    p.add_argument("-list", "--list", dest="list", default="jcvi_list.txt",
                   help="Genome-pair list, two IDs per line (default: jcvi_list.txt)")
    p.add_argument("-max_ks", "--max_ks", dest="max_ks", type=float, default=2.0,
                   help="Drop gene pairs with Ks >= this (saturation cutoff; default: 2.0)")
    p.add_argument("-genome_out", "--genome_out", dest="genome_out", default="ks_genome.tsv",
                   help="Per-genome-pair medians (default: ks_genome.tsv)")
    p.add_argument("-long_out", "--long_out", dest="long_out", default="ks_all.tsv",
                   help="Long-format per-gene Ks for plotting (default: ks_all.tsv)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Report per-pair retention counts")
    return p.parse_args()


def read_pairs(path):
    """Genome pairs from jcvi_list.txt. IDs may contain dots, so never split on '.'."""
    pairs = []
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) >= 2:
                pairs.append((fields[0], fields[1]))
    if not pairs:
        sys.exit(f"Error: no genome pairs read from {path}")
    return pairs


def read_ks(path, max_ks, verbose):
    """Ks column of a merged KaKs table. Read by header name -- KaKs_Calculator's column
    order shifts between methods."""
    kept, seen, dropped = [], 0, 0
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "Ks" not in reader.fieldnames:
            sys.exit(f"Error: {path} has no 'Ks' column (found: {reader.fieldnames})")
        for row in reader:
            seen += 1
            try:
                ks = float(row["Ks"])
            except (TypeError, ValueError):
                dropped += 1          # 'NA'/'nan' -- KaKs could not fit the model
                continue
            if not (0.0 < ks < max_ks):
                dropped += 1
                continue
            kept.append(ks)
    if verbose:
        print(f"  {path}: {len(kept)}/{seen} gene pairs kept "
              f"(0 < Ks < {max_ks}), {dropped} dropped", file=sys.stderr)
    return kept


def main():
    args = parse_args()
    pairs = read_pairs(args.list)

    with open(args.genome_out, "w") as gout, open(args.long_out, "w") as lout:
        lout.write("pair\tid1\tid2\tKs\n")
        for id1, id2 in pairs:
            kaks_file = f"{id1}.{id2}.kaks.tsv"
            ks_values = read_ks(kaks_file, args.max_ks, args.verbose)
            if not ks_values:
                sys.exit(f"Error: no usable Ks values in {kaks_file}. "
                         f"Every gene pair was NA or fell outside (0, {args.max_ks}).")

            median = statistics.median(ks_values)
            gout.write(f"{id1}\t{id2}\t{len(ks_values)}\t{median:.6f}\n")

            label = f"{id1} vs {id2}"
            for ks in ks_values:
                lout.write(f"{label}\t{id1}\t{id2}\t{ks:.6f}\n")

            print(f"{id1} vs {id2}: median Ks = {median:.4f} "
                  f"(n = {len(ks_values)} gene pairs)", file=sys.stderr)

    print(f"Wrote {args.genome_out} and {args.long_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
