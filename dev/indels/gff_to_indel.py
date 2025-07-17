#!/usr/bin/env python3
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Print GFF lines whose positions match indel records (with optional wiggle) "
            "and whose chromosome matches indel column 2 or column 5. "
            "Appends matching indel lines to the output."
        )
    )
    parser.add_argument(
        "-indel", required=True,
        help="Path to indel.tsv file"
    )
    parser.add_argument(
        "-gff", required=True,
        help="Path to GFF file"
    )
    wiggle_group = parser.add_mutually_exclusive_group()
    wiggle_group.add_argument(
        "--wiggle", type=int, default=None,
        help="Fixed wiggle in bases (± N bp) for matching positions"
    )
    wiggle_group.add_argument(
        "--wiggle-fraction", type=float, default=None,
        help="Fraction of feature length to use as dynamic wiggle (e.g. 0.1 for 10% of feature length)"
    )
    return parser.parse_args()

def build_indel_index(indel_path):
    """
    Build a dict mapping each chromosome name to a list of
    (start, end, full_indel_line) tuples, indexing by both cols 2 and 5.
    """
    d = defaultdict(list)
    with open(indel_path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                chrom1 = fields[1]
                start = int(fields[2])
                end = int(fields[3])
                chrom2 = fields[4]
            except ValueError:
                continue
            d[chrom1].append((start, end, line))
            if chrom2 != chrom1:
                d[chrom2].append((start, end, line))
    return d

def main():
    args = parse_args()
    indel_index = build_indel_index(args.indel)

    with open(args.gff) as gff:
        for raw in gff:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            chrom = cols[0]
            try:
                g_start = int(cols[3])
                g_end   = int(cols[4])
            except ValueError:
                continue

            if chrom not in indel_index:
                continue

            # determine threshold
            if args.wiggle is not None:
                threshold = args.wiggle
            elif args.wiggle_fraction is not None:
                feature_len = g_end - g_start
                threshold = args.wiggle_fraction * feature_len
            else:
                threshold = 0

            matches = []
            for i_start, i_end, indel_line in indel_index[chrom]:
                if abs(i_start - g_start) <= threshold and abs(i_end - g_end) <= threshold:
                    matches.append(indel_line)

            if matches:
                print(f"{line}\t{' | '.join(matches)}")

if __name__ == "__main__":
    main()
