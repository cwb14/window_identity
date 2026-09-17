"""Block orientation labels in the coord extractors (Step 7).

Orientation must come from ALL anchors of a block, not from its first/last anchor pair:
on real Poa data a single minority anchor at a block end mislabelled whole-chromosome,
perfectly collinear blocks (e.g. 1,747 anchors, 98% opposite-strand, labelled '+').
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import gene_coords_extractor_all4 as blk  # noqa: E402
import gene_coords_extractor_all4_pairs as gp  # noqa: E402


def make_block(strand_pairs, ascending_b=True, step=10_000):
    """Synthetic collinear block on A_chr1 / B_chr1.

    strand_pairs: list of (strand_in_A, strand_in_B), one per anchor, in A order.
    Returns (gene_coords, cluster) in the shapes the extractors consume.
    """
    coords, cluster = {}, []
    n = len(strand_pairs)
    for i, (sa, sb) in enumerate(strand_pairs):
        a, b = f"a{i}", f"b{i}"
        pos_b = i if ascending_b else n - 1 - i
        coords[a] = ("A_chr1", i * step, i * step + 1000, sa)
        coords[b] = ("B_chr1", pos_b * step, pos_b * step + 1000, sb)
        cluster.append([a, b])
    return coords, cluster


def strands_of(lines):
    return [ln.split("\t")[2] for ln in lines]


# ---------------------------------------------------------------- block partition

def test_block_minority_first_anchor_does_not_decide_label():
    # First anchor pair is opposite-strand, the other 9 agree -> collinear '+' block.
    pairs = [("+", "-")] + [("+", "+")] * 9
    coords, cluster = make_block(pairs)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["+"]


def test_block_majority_inverted_even_when_first_and_last_agree_plus():
    # First and last pairs say '+', but 8 of 10 anchors are opposite-strand.
    pairs = [("+", "+")] + [("+", "-")] * 8 + [("-", "-")]
    coords, cluster = make_block(pairs, ascending_b=False)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["-"]


def test_block_uniform_orientations_unchanged():
    coords, cluster = make_block([("+", "+"), ("-", "-"), ("+", "+")])
    assert strands_of(blk.process_clusters([cluster], coords)) == ["+"]
    coords, cluster = make_block([("+", "-"), ("-", "+"), ("+", "-")], ascending_b=False)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["-"]


def test_block_strand_tie_broken_by_anchor_order():
    tie = [("+", "+"), ("+", "-")]
    coords, cluster = make_block(tie, ascending_b=True)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["+"]
    tie = [("+", "-"), ("+", "+")]          # first anchor says '-', order says '+'
    coords, cluster = make_block(tie, ascending_b=True)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["+"]
    coords, cluster = make_block(tie, ascending_b=False)
    assert strands_of(blk.process_clusters([cluster], coords)) == ["-"]


def test_block_orientation_helper_reports_vote():
    coords, cluster = make_block([("+", "-")] + [("+", "+")] * 3)
    anchors = [(coords[a], coords[b]) for a, b in cluster]
    assert blk.block_orientation(anchors) == "+"


# ---------------------------------------------------------------- genepair partition

def test_genepair_isolated_flipped_anchor_uses_block_majority():
    # Anchor 2 of 5 is opposite-strand inside a '+' block. Both segments touching it
    # have disagreeing ends; the old first-anchor rule labelled the one starting at
    # anchor 2 as '-'. Every segment of this collinear block should be '+'.
    pairs = [("+", "+"), ("+", "+"), ("+", "-"), ("+", "+"), ("+", "+")]
    coords, cluster = make_block(pairs)
    block = [tuple(p) for p in cluster]
    assert strands_of(gp.process_pairs([block], coords)) == ["+"] * 4


def test_genepair_agreeing_segments_keep_their_own_orientation():
    # A locally consistent inverted run inside a '+'-majority block keeps '-'.
    pairs = [("+", "+")] * 4 + [("+", "-")] * 2
    coords, cluster = make_block(pairs)
    block = [tuple(p) for p in cluster]
    assert strands_of(gp.process_pairs([block], coords)) == ["+", "+", "+", "+", "-"]
