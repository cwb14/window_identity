#!/usr/bin/env python3
import sys
import argparse

def process_file(infile: str, threshold: int, stitch_gaps: bool,
                 max_stitch_ratio: float = 5.0, stitch_flank_factor: float = 1.0):
    """
    Process the anchors coordinate file to merge lines based on overlapping
    or touching (butt heads) sequence ranges, per (pair1_id, pair2_id, strand, block) bin.

    A pair is considered "pass" if both lengths >= threshold; otherwise "fail".
    Only lines where at least one side is "fail" are eligible to merge, so a merge
    never crosses a syntenic block.

    If --stitch-gaps is set, after merging we insert synthetic lines to fill the gaps
    between consecutive syntenic BLOCKS of the same (pair1_id, pair2_id, strand), whenever
    both sequences have a positive gap. A block's extent is the span of its records, so this
    behaves the same under -partition block (a block is one record) and genepair (a block is
    many gene-pair segments, and only its outer edges bound a real gap). Touching blocks are
    not stitched. Guards: no opposite-strand block in the gap, consecutive in both genomes,
    gap sizes within max_stitch_ratio of each other, and each gap no longer than
    stitch_flank_factor x the smaller neighbouring block on that genome.
    """
    # ---------- parsing ----------
    bins = {}

    def parse_range(token: str):
        # e.g., "9311v2_chr1:18291..25404" -> ("9311v2_chr1", 18291, 25404)
        id_part, coords = token.split(':', 1)
        start_str, end_str = coords.split('..', 1)
        return id_part, int(start_str), int(end_str)

    with open(infile, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines quietly

            pair1_id, pair1_start, pair1_end = parse_range(parts[0])
            pair2_id, pair2_start, pair2_end = parse_range(parts[1])
            strand = parts[2].strip()
            # Column 4 is the syntenic block id from the coord extractors. It goes
            # into the bin key so a merge can never span two blocks. Files written
            # before this column existed have three columns; they bin as before.
            block = parts[3].strip() if len(parts) > 3 else ""

            # Calculate lengths (half-open style; matches original behavior)
            len1 = pair1_end - pair1_start
            len2 = pair2_end - pair2_start

            threshold_status = "pass" if (len1 >= threshold and len2 >= threshold) else "fail"
            bin_key = f"{pair1_id}_{pair2_id}_{'plus' if strand == '+' else 'minus'}_{block}"

            new_line = [
                pair1_id, pair1_start, pair1_end,   # 0..2
                pair2_id, pair2_start, pair2_end,   # 3..5
                strand,                             # 6
                len1, len2,                         # 7..8
                bin_key,                            # 9
                threshold_status,                   # 10
                block                               # 11
            ]
            bins.setdefault(bin_key, []).append(new_line)

    # ---------- merging ----------
    def ranges_touch_or_overlap(a_start, a_end, b_start, b_end):
        # Allow touching: end == start
        return (a_start <= b_end) and (a_end >= b_start)

    final_bins = {}

    for bin_key, lines in bins.items():
        # Sort to make merging deterministic and efficient
        lines.sort(key=lambda x: (x[1], x[4]))  # by pair1_start then pair2_start
        merged_list = []

        while lines:
            current_line = lines.pop(0)
            merged = False

            for i, line in enumerate(lines):
                cond1 = ranges_touch_or_overlap(current_line[1], current_line[2], line[1], line[2])
                cond2 = ranges_touch_or_overlap(current_line[4], current_line[5], line[4], line[5])

                # Merge only while the record being placed is still under the
                # threshold. The old condition also fired when the *other* record
                # was "fail", and unconditionally stamped the result "pass" -- so a
                # record that had already reached the threshold went on absorbing
                # every small neighbour. Since consecutive gene-pair intervals share
                # an anchor gene they always touch, and the median interval is far
                # below the threshold, that merge ran transitively down a whole
                # chromosome arm (4.12 Gb of segment span over a 1.12 Gb genome,
                # single segments up to 159 Mb).
                #
                # Recomputing the status instead of assuming "pass" lets a merged
                # record keep growing only until it genuinely reaches the threshold,
                # which bounds a segment at roughly threshold + one interval. It
                # terminates because every merge strictly increases both lengths.
                if cond1 and cond2 and current_line[10] == "fail":
                    new_len1 = max(current_line[2], line[2]) - min(current_line[1], line[1])
                    new_len2 = max(current_line[5], line[5]) - min(current_line[4], line[4])
                    merged_line = [
                        current_line[0],
                        min(current_line[1], line[1]),
                        max(current_line[2], line[2]),
                        current_line[3],
                        min(current_line[4], line[4]),
                        max(current_line[5], line[5]),
                        current_line[6],
                        new_len1, new_len2,
                        current_line[9],
                        "pass" if (new_len1 >= threshold and new_len2 >= threshold) else "fail",
                        current_line[11]
                    ]
                    lines[i] = merged_line
                    merged = True
                    break

            if not merged:
                merged_list.append(current_line)

        # Keep per-bin results for optional stitching
        final_bins[bin_key] = merged_list

    # ---------- optional stitching ----------
    if stitch_gaps:
        from collections import defaultdict

        # Stitching works on BLOCK EXTENTS, not on single records. Under -partition block a
        # block IS one record, so this changes nothing there. Under genepair a block is many
        # gene-pair segments; only the block's outer edges bound a real gap, and the segments
        # bordering one are too small (and too locally scrambled) to judge it by themselves.
        extents = {}
        for lines in final_bins.values():
            for rec in lines:
                # A record with no block id (legacy 3-column coords) is its own block.
                key = (rec[0], rec[3], rec[11] or f"__rec{id(rec)}")
                e = extents.get(key)
                if e is None:
                    e = extents[key] = {"p1": rec[0], "p2": rec[3], "blk": rec[11],
                                        "a0": rec[1], "a1": rec[2], "b0": rec[4], "b1": rec[5],
                                        "plus": 0, "minus": 0, "bin": rec[9]}
                e["a0"] = min(e["a0"], rec[1]); e["a1"] = max(e["a1"], rec[2])
                e["b0"] = min(e["b0"], rec[4]); e["b1"] = max(e["b1"], rec[5])
                e["plus" if rec[6] == "+" else "minus"] += rec[2] - rec[1]
        for e in extents.values():
            # bp-weighted majority, matching how riparian judges a block's orientation.
            e["strand"] = "+" if e["plus"] >= e["minus"] else "-"

        by_pair = defaultdict(list)
        for e in extents.values():
            by_pair[(e["p1"], e["p2"])].append(e)

        def has_opposite_strand_between(prev, nxt):
            gap1_start, gap1_end = prev["a1"], nxt["a0"]
            gap2_start, gap2_end = prev["b1"], nxt["b0"]
            if not (gap1_end > gap1_start and gap2_end > gap2_start):
                return False  # nothing to check
            for cand in by_pair[(prev["p1"], prev["p2"])]:
                if cand is prev or cand is nxt or cand["strand"] == prev["strand"]:
                    continue
                # an opposite-strand block overlapping BOTH gap windows "occupies" the gap
                if (cand["a0"] <= gap1_end and cand["a1"] >= gap1_start) and \
                   (cand["b0"] <= gap2_end and cand["b1"] >= gap2_start):
                    return True
            return False

        groups = defaultdict(list)
        for e in extents.values():
            groups[(e["p1"], e["p2"], e["strand"])].append(e)

        for group in groups.values():
            if len(group) < 2:
                continue
            by_a = sorted(group, key=lambda e: (e["a0"], e["b0"]))
            by_b = sorted(group, key=lambda e: (e["b0"], e["a0"]))
            pos_in_b = {id(e): i for i, e in enumerate(by_b)}

            for prev, nxt in zip(by_a, by_a[1:]):
                # Positive gaps on both sequences?
                gap1_start, gap1_end = prev["a1"], nxt["a0"]
                gap2_start, gap2_end = prev["b1"], nxt["b0"]
                if not (gap1_end > gap1_start and gap2_end > gap2_start):
                    continue  # no stitch for touching/overlapping

                # GUARD #1: an inversion sits in the gap
                if has_opposite_strand_between(prev, nxt):
                    continue

                # GUARD #2: the two blocks must be consecutive in BOTH genomes
                if abs(pos_in_b[id(prev)] - pos_in_b[id(nxt)]) != 1:
                    continue

                # GUARD #3: the two sides of the gap must be of comparable size. A 50 Mb
                # against 5 Mb stitch asserts a correspondence the anchors never showed, and
                # the record is unalignable end to end.
                g1 = gap1_end - gap1_start
                g2 = gap2_end - gap2_start
                if max_stitch_ratio > 0 and max(g1, g2) > max_stitch_ratio * min(g1, g2):
                    continue

                # GUARD #4: the gap must be backed by synteny on both sides -- no longer than
                # stitch_flank_factor x the smaller neighbouring BLOCK, on each genome.
                # Without it, two tiny spurious blocks on a non-homologous chromosome pair
                # were stitched across 75 Mb, and the fake block then displaced real (rescued)
                # blocks in the riparian plot. Same principle as paf_chain_blocks.py
                # --max-gap-factor.
                if stitch_flank_factor > 0 and (
                        g1 > stitch_flank_factor * min(prev["a1"] - prev["a0"],
                                                       nxt["a1"] - nxt["a0"]) or
                        g2 > stitch_flank_factor * min(prev["b1"] - prev["b0"],
                                                       nxt["b1"] - nxt["b0"])):
                    continue

                # All guards pass: insert the synthetic record into the left block's bin, so it
                # prints among that block's records and the order of existing records is kept.
                final_bins[prev["bin"]].append([
                    prev["p1"], gap1_start, gap1_end,
                    prev["p2"], gap2_start, gap2_end,
                    prev["strand"],
                    g1, g2,
                    prev["bin"],
                    "stitched",  # internal marker; output format ignores this
                    prev["blk"]  # block id of the left neighbour; nothing downstream reads it
                ])

    # ---------- output ----------
    # Flatten bins in insertion order; within bin keep sorted order for determinism
    for bin_key in final_bins:
        for line in sorted(final_bins[bin_key], key=lambda x: (x[1], x[4])):
            row = (f"{line[0]}:{line[1]}..{line[2]}\t"
                   f"{line[3]}:{line[4]}..{line[5]}\t{line[6]}")
            if line[11]:
                row += f"\t{line[11]}"
            print(row)

def main():
    parser = argparse.ArgumentParser(
        description="Merge anchor coordinate lines by overlapping/touching ranges within (pair1_id, pair2_id, strand) bins, with optional gap stitching."
    )
    parser.add_argument(
        "infile",
        help="Path to the input anchors coordinate file (e.g., all.recip.anchors.coords)"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=1_000_000,
        help="Length threshold for each side to be considered 'pass' (default: 1000000)"
    )
    parser.add_argument(
        "--stitch-gaps",
        action="store_true",
        help="After merging, insert synthetic lines to fill positive gaps between consecutive records within each (pair1, pair2, strand) bin."
    )
    parser.add_argument(
        "--max-stitch-ratio",
        type=float,
        default=5.0,
        help="Do not stitch a gap pair whose two sides differ by more than this "
             "factor; 0 disables the check (default: %(default)s, the same limit Step 10 "
             "and riparian apply to skewed segments). A 50 Mb vs 5 Mb stitch asserts a "
             "correspondence the anchors never showed, and the resulting segment is "
             "unalignable end-to-end."
    )
    parser.add_argument(
        "--stitch-flank-factor",
        type=float,
        default=1.0,
        help="Do not stitch a gap longer than this multiple of the smaller neighbouring "
             "block, on either genome; 0 disables the check (default: %(default)s). Keeps "
             "two tiny spurious blocks from being bridged across tens of Mb."
    )
    args = parser.parse_args()

    if args.threshold < 0:
        print("Threshold must be non-negative.", file=sys.stderr)
        sys.exit(2)

    process_file(args.infile, args.threshold, args.stitch_gaps,
                 args.max_stitch_ratio, args.stitch_flank_factor)

if __name__ == "__main__":
    main()
