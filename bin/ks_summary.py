#!/usr/bin/env python3
"""Summarise per-gene Ks into one distance per genome pair.

Reads the merged KaKs_Calculator tables ('{ID1}.{ID2}.kaks.tsv') for every pair listed in
jcvi_list.txt and emits:

  -genome_out   ID1  ID2  n_pairs  median_Ks   (the 4-column layout matrix_builder.py eats)
  -long_out     pair  id1  id2  Ks             (one row per gene pair, for the density plot)

The genome-pair distance is the MEDIAN Ks over gene pairs. Ks has a long right tail from
paralogous and saturated anchors, which drags the mean upward; the median does not move.

Gene pairs with 0 <= Ks < max_ks are kept. Above ~2 Ks is saturated, where the estimator is
unstable and effectively unbounded.

Ks 'NA' is read as 0 and kept. KaKs_Calculator 3.0 prints any Ks estimate below 1e-6 as 'NA'
(base.cpp, parseOutput). That is a display convention, not a failure flag:
  * no method ever assigns Ks the tool's internal NA sentinel;
  * with zero synonymous differences, YN00 cannot fit F84/K80 and falls back to Jukes-Cantor,
    which gives exactly 0 (DistanceF84);
  * saturation yields large finite values (capped at 99), not 'NA'.
Dropping every 'NA' used to discard the most similar gene pairs -- 51% of PaA/PiA anchors -- and
biased recent divergences upward.

A per-gene Ks of 0 is still only "below detection" (roughly < 1/S synonymous sites), and Ka/Ks
is genuinely undefined there. That is a reason to be cautious about a single gene, not to drop
it from a genome-wide distribution.
"""

import argparse
import csv
import statistics
import sys


# KaKs_Calculator 3.0 writes its result rows with NO header line, so a merged table has none
# either: the shell's merge step takes 'head -1' of the first .kaks file as the header, and
# that line is already data. DictReader would then consume a data row as the header and never
# find a 'Ks' column. These are the column names the tool itself uses internally, verbatim
# from KaKs_Calculator-3.0/bin/KaKs.cpp:21-25, and are applied when no header is present.
KAKS_COLUMNS = [
    "Sequence", "Method", "Ka", "Ks", "Ka/Ks",
    "P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
    "Substitutions", "Syn-Subs", "Nonsyn-Subs",
    "Fold-Syn-Subs(0:2:4)", "Fold-Nonsyn-Subs(0:2:4)",
    "Divergence-Distance", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)",
    "GC(1:2:3)", "ML-Score", "AICc", "Akaike-Weight", "Model",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarise per-gene Ks into one median distance per genome pair."
    )
    p.add_argument("-list", "--list", dest="list", default="jcvi_list.txt",
                   help="Genome-pair list, two IDs per line (default: jcvi_list.txt)")
    p.add_argument("-max_ks", "--max_ks", dest="max_ks", type=float, default=2.0,
                   help="Drop gene pairs with Ks >= this (saturation cutoff; default: 2.0). "
                        "Ks = 0 is kept")
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


def parse_ks(row):
    """Ks of one KaKs_Calculator row: 'NA' is 0.0; anything unparseable is None."""
    ks = row["Ks"]
    if ks == "NA":
        return 0.0
    try:
        return float(ks)
    except (TypeError, ValueError):
        return None


def iter_ks(path):
    """(Sequence, Ks-or-None) for every row of a merged KaKs table. Read by header name --
    KaKs_Calculator's column order shifts between methods."""
    with open(path, newline="") as fh:
        first = fh.readline()
        if not first.strip():
            sys.exit(f"Error: {path} is empty.")
        fh.seek(0)
        fields = first.rstrip("\n").split("\t")
        if "Ks" in fields:
            reader = csv.DictReader(fh, delimiter="\t")
        elif len(fields) == len(KAKS_COLUMNS):
            # Headerless KaKs_Calculator 3.0 output -- name the columns ourselves and treat
            # every line, including the first, as data.
            reader = csv.DictReader(fh, delimiter="\t", fieldnames=KAKS_COLUMNS)
        else:
            sys.exit(f"Error: {path} has no 'Ks' column and its {len(fields)} columns do not "
                     f"match the {len(KAKS_COLUMNS)}-column KaKs_Calculator layout. "
                     f"First fields: {fields[:5]}")
        for row in reader:
            yield row["Sequence"], parse_ks(row)


def read_ks(path, max_ks, verbose):
    """Ks values with 0 <= Ks < max_ks; unparseable rows and saturated pairs are dropped."""
    kept, seen, bad, saturated = [], 0, 0, 0
    for _, ks in iter_ks(path):
        seen += 1
        if ks is None:
            bad += 1
        elif not (0.0 <= ks < max_ks):
            saturated += 1
        else:
            kept.append(ks)
    if bad:
        print(f"WARNING: {path}: skipped {bad} row(s) with an unparseable Ks", file=sys.stderr)
    if verbose:
        print(f"  {path}: {len(kept)}/{seen} gene pairs kept (0 <= Ks < {max_ks}; "
              f"{sum(k == 0 for k in kept)} with Ks = 0), {saturated} saturated dropped",
              file=sys.stderr)
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
                         f"Every gene pair was unparseable or had Ks >= {args.max_ks}.")

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
