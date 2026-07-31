#!/usr/bin/env python3
"""Chain collinear PAF records into blocks.

The aligner emits one record per alignment, so a chromosome arm arrives as tens of
thousands of fragments even when it is perfectly collinear. AnchorWave reports the
same sequence as a handful of blocks (51 genome-wide on the Pinfirma/PannuaA pair,
the largest spanning 214 Mb) purely because its chaining charges no penalty for a
large gap between anchors. This recovers that view without touching the aligner:
it is a read-only pass over the PAF that groups records which are already
collinear, and reports the block structure alongside per-block identity.

It does not modify or re-align anything. The PAF stays the record of what aligned;
this says how those records line up.

Reads the PAF streaming and keeps only the 12 mandatory columns, so a multi-GB
file costs a few tens of MB.
"""

import argparse
import sys
from collections import defaultdict


def parse_paf(path, min_len=0):
    """(t_name, q_name, strand) -> [(t_start, t_end, q_start, q_end, nmatch, alnlen)]"""
    groups = defaultdict(list)
    n = 0
    with open(path) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if len(f) < 12:
                continue
            try:
                q_start, q_end = int(f[2]), int(f[3])
                t_start, t_end = int(f[7]), int(f[8])
                nmatch, alnlen = int(f[9]), int(f[10])
            except ValueError:
                continue
            if t_end - t_start < min_len:
                continue
            n += 1
            groups[(f[5], f[0], f[4])].append(
                (t_start, t_end, q_start, q_end, nmatch, alnlen))
    return groups, n


def chain_group(recs, max_gap, slack, strand, gap_factor=1.0):
    """Split one (target, query, strand) group into collinear runs.

    Collinear means the query advances in step with the target -- forward for '+',
    backward for '-'. A run breaks when either axis jumps more than max_gap, when
    the query moves the wrong way by more than `slack` (which absorbs the small
    overlaps produced by segments that share an anchor gene), or when the jump is
    large relative to how much the chain has actually aligned.

    That last rule is what keeps the output honest. With only an absolute gap
    limit, two 5 kb records 18 Mb apart chain into an "18.3 Mb block" -- and on
    real data that produced blocks totalling 2,030 Mb over a 1.12 Gb genome, i.e.
    mostly spurious. Requiring a chain to have aligned roughly as much as it wants
    to skip means a long block has to be supported by dense alignment, which is
    exactly what makes AnchorWave's megabase blocks real rather than an artifact
    of a permissive gap penalty.
    """
    recs.sort(key=lambda r: (r[0], r[2]))
    blocks = []
    cur = None
    for r in recs:
        t_start, t_end, q_start, q_end, nmatch, alnlen = r
        if cur is None:
            cur = dict(t_lo=t_start, t_hi=t_end, q_lo=q_start, q_hi=q_end,
                       n=1, matched=nmatch, alnlen=alnlen)
            continue
        t_gap = t_start - cur["t_hi"]
        if strand == "+":
            q_gap = q_start - cur["q_hi"]
        else:
            q_gap = cur["q_lo"] - q_end
        earned = gap_factor * cur["alnlen"] if gap_factor > 0 else float("inf")
        if (t_gap > max_gap or q_gap > max_gap
                or t_gap < -slack or q_gap < -slack
                or t_gap > earned or q_gap > earned):
            blocks.append(cur)
            cur = dict(t_lo=t_start, t_hi=t_end, q_lo=q_start, q_hi=q_end,
                       n=1, matched=nmatch, alnlen=alnlen)
            continue
        cur["t_hi"] = max(cur["t_hi"], t_end)
        cur["t_lo"] = min(cur["t_lo"], t_start)
        cur["q_hi"] = max(cur["q_hi"], q_end)
        cur["q_lo"] = min(cur["q_lo"], q_start)
        cur["n"] += 1
        cur["matched"] += nmatch
        cur["alnlen"] += alnlen
    if cur is not None:
        blocks.append(cur)
    return blocks


def main():
    p = argparse.ArgumentParser(
        description="Chain collinear PAF records into blocks (read-only).")
    p.add_argument("paf", help="Input PAF (e.g. alignment_adjust.paf)")
    p.add_argument("-o", "--output", default="-",
                   help="Output TSV, '-' for stdout (default: %(default)s)")
    p.add_argument("--max-gap", type=int, default=1_000_000,
                   help="Break a chain when either axis jumps more than this "
                        "(default: %(default)s). Raise it for more "
                        "chromosome-scale blocks; AnchorWave effectively uses no "
                        "bp limit at all, which is why its blocks reach 214 Mb.")
    p.add_argument("--max-gap-factor", type=float, default=1.0,
                   help="A chain may only jump a gap up to this multiple of the "
                        "sequence it has already aligned; 0 disables (default: "
                        "%(default)s). Without it, two small distant records chain "
                        "into a bogus megabase block.")
    p.add_argument("--min-aln-frac", type=float, default=0.0,
                   help="Drop blocks whose aligned length is below this fraction of "
                        "their reference span (default: %(default)s).")
    p.add_argument("--slack", type=int, default=1000,
                   help="Backward movement tolerated before a chain breaks, to "
                        "absorb overlaps between segments sharing an anchor gene "
                        "(default: %(default)s).")
    p.add_argument("--min-record-len", type=int, default=0,
                   help="Ignore records shorter than this on the target axis.")
    p.add_argument("--min-block-len", type=int, default=0,
                   help="Do not report blocks shorter than this on the target axis.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    groups, n_rec = parse_paf(args.paf, args.min_record_len)
    if args.verbose:
        print(f"[chain] {n_rec:,} records in {len(groups):,} "
              f"(target, query, strand) groups", file=sys.stderr)

    rows = []
    for (t_name, q_name, strand), recs in groups.items():
        for b in chain_group(recs, args.max_gap, args.slack, strand,
                             args.max_gap_factor):
            span = b["t_hi"] - b["t_lo"]
            if span < args.min_block_len:
                continue
            if args.min_aln_frac > 0 and span > 0 and \
                    b["alnlen"] / span < args.min_aln_frac:
                continue
            rows.append((t_name, b["t_lo"], b["t_hi"], q_name, b["q_lo"], b["q_hi"],
                         strand, b["n"], span, b["q_hi"] - b["q_lo"],
                         b["matched"], b["alnlen"]))
    rows.sort(key=lambda r: (r[0], r[1]))

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        out.write("\t".join([
            "ref_name", "ref_start", "ref_end", "qry_name", "qry_start", "qry_end",
            "strand", "n_records", "ref_span", "qry_span", "matched", "aln_len",
            "identity"]) + "\n")
        for r in rows:
            ident = r[10] / r[11] if r[11] else 0.0
            out.write("\t".join(str(x) for x in r) + f"\t{ident:.6f}\n")
    finally:
        if out is not sys.stdout:
            out.close()

    if args.verbose and rows:
        spans = sorted((r[8] for r in rows), reverse=True)
        tot = sum(spans)
        cum = 0
        n50 = spans[-1]
        for s in spans:
            cum += s
            if cum >= tot / 2:
                n50 = s
                break
        print(f"[chain] {len(rows):,} blocks, {n_rec:,} records "
              f"({n_rec/len(rows):.1f} per block); block N50 {n50:,} bp, "
              f"largest {spans[0]:,} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
