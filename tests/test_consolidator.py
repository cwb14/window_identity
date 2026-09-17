"""anchor_coord_consolidator.py: merging is confined to a block; gap stitching is not.

Commit 829e382 put the block id into the bin key so a merge could never run across blocks.
Stitching iterated over those same bins, so it only ever saw records of one block -- but gaps
exist only BETWEEN blocks, and stitching silently stopped doing anything. Fixtures below are
the PaA/PaB chr1 blocks where that was noticed.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import anchor_coord_consolidator as cons  # noqa: E402

M = 1_000_000


def rec(a0, a1, b0, b1, strand, block, ca="A_chr1", cb="B_chr1"):
    return f"{ca}:{int(a0)}..{int(a1)}\t{cb}:{int(b0)}..{int(b1)}\t{strand}\t{block}"


def run(tmp_path, capsys, lines, **kw):
    p = tmp_path / "in.coords"
    p.write_text("\n".join(lines) + "\n")
    cons.process_file(str(p), kw.pop("threshold", 15000), kw.pop("stitch_gaps", True), **kw)
    return capsys.readouterr().out.splitlines()


def spans(out):
    """(A start, A end, B start, B end, strand) per output record, in Mb."""
    res = []
    for ln in out:
        f = ln.split("\t")
        a = f[0].split(":")[1].split("..")
        b = f[1].split(":")[1].split("..")
        res.append((int(a[0]) / M, int(a[1]) / M, int(b[0]) / M, int(b[1]) / M, f[2]))
    return res


# PaA_PaB chr1, blocks 2 and 3: same strand, gap 2.00 Mb / 1.80 Mb.
BLOCK2 = rec(175.39 * M, 212.55 * M, 46.42 * M, 51.54 * M, "+", 2)
BLOCK3 = rec(214.55 * M, 284.37 * M, 53.34 * M, 71.36 * M, "+", 3)


def test_gap_between_adjacent_blocks_is_stitched(tmp_path, capsys):
    out = spans(run(tmp_path, capsys, [BLOCK2, BLOCK3]))
    assert (212.55, 214.55, 51.54, 53.34, "+") in out
    assert len(out) == 3


def test_stitching_off_adds_nothing(tmp_path, capsys):
    out = run(tmp_path, capsys, [BLOCK2, BLOCK3], stitch_gaps=False)
    assert len(out) == 2


def test_default_ratio_guard_is_5x(tmp_path, capsys):
    # PaA_PaB chr1 blocks 1 -> 2: gap 2.18 Mb / 0.70 Mb = 3.1x, rejected by the old 3x default.
    block1 = rec(0.07 * M, 173.21 * M, 0.35 * M, 45.71 * M, "+", 1)
    out = spans(run(tmp_path, capsys, [block1, BLOCK2]))
    assert (173.21, 175.39, 45.71, 46.42, "+") in out
    # 6x mismatch stays unstitched.
    far = rec(212.55 * M + 6 * M, 250 * M, 51.54 * M + M, 60 * M, "+", 3)
    assert len(run(tmp_path, capsys, [BLOCK2, far])) == 2


def test_opposite_strand_neighbours_are_not_stitched(tmp_path, capsys):
    # Blocks 3 (+), 4 (-), 5 (+): no stitch across the strand change, and the + blocks
    # either side of the inversion are not bridged over it.
    block4 = rec(284.45 * M, 301.65 * M, 71.41 * M, 82.17 * M, "-", 4)
    block5 = rec(301.82 * M, 319.55 * M, 82.20 * M, 97.69 * M, "+", 5)
    assert len(run(tmp_path, capsys, [BLOCK3, block4, block5])) == 3


def test_blocks_out_of_order_in_second_genome_are_not_stitched(tmp_path, capsys):
    # A rearrangement: consecutive in genome 1, not in genome 2.
    x = rec(10 * M, 20 * M, 50 * M, 60 * M, "+", 1)
    y = rec(22 * M, 30 * M, 70 * M, 80 * M, "+", 2)
    z = rec(32 * M, 40 * M, 62 * M, 68 * M, "+", 3)
    out = spans(run(tmp_path, capsys, [x, y, z]))
    assert (20.0, 22.0, 60.0, 70.0, "+") not in out


def test_merging_still_never_crosses_a_block(tmp_path, capsys):
    # Touching sub-threshold gene-pair segments: merged within a block, never across.
    same = [rec(0, 5000, 0, 5000, "+", 1), rec(5000, 9000, 5000, 9000, "+", 1)]
    assert len(run(tmp_path, capsys, same, stitch_gaps=False)) == 1
    across = [rec(0, 5000, 0, 5000, "+", 1), rec(5000, 9000, 5000, 9000, "+", 2)]
    assert len(run(tmp_path, capsys, across, stitch_gaps=False)) == 2


def test_record_order_unchanged_when_nothing_is_stitched(tmp_path, capsys):
    lines = [BLOCK3, rec(1 * M, 2 * M, 1 * M, 2 * M, "-", 9, "A_chr2", "B_chr2"), BLOCK2]
    assert run(tmp_path, capsys, lines, stitch_gaps=False) == \
        run(tmp_path, capsys, lines, stitch_gaps=True, max_stitch_ratio=1.0)


# ---------------------------------------------------------------- flank guard

def test_gap_longer_than_smaller_flank_is_not_stitched(tmp_path, capsys):
    # Two ~0.1 Mb spurious blocks 75 Mb apart (the PaA_chr6/PaB_chr1 case from the sweep):
    # same strand, consecutive, comparable gap sizes -- but no synteny to support the jump.
    left = rec(6.87 * M, 6.97 * M, 26.70 * M, 26.80 * M, "+", 1, "A_chr6", "B_chr1")
    right = rec(82.24 * M, 82.34 * M, 97.79 * M, 97.89 * M, "+", 2, "A_chr6", "B_chr1")
    assert len(run(tmp_path, capsys, [left, right])) == 2


def test_flank_guard_applies_on_both_genomes(tmp_path, capsys):
    # Gap fits the flanks on genome 1 (1 Mb vs 10 Mb blocks) but not on genome 2
    # (1 Mb gap vs a 0.5 Mb block).
    left = rec(0, 10 * M, 0, 0.5 * M, "+", 1)
    right = rec(11 * M, 21 * M, 1.5 * M, 11.5 * M, "+", 2)
    assert len(run(tmp_path, capsys, [left, right], max_stitch_ratio=0)) == 2


def test_gap_equal_to_smaller_flank_is_stitched(tmp_path, capsys):
    left = rec(0, 2 * M, 0, 2 * M, "+", 1)
    right = rec(4 * M, 10 * M, 4 * M, 10 * M, "+", 2)
    out = spans(run(tmp_path, capsys, [left, right]))
    assert (2.0, 4.0, 2.0, 4.0, "+") in out


def test_flank_factor_is_tunable_and_zero_disables(tmp_path, capsys):
    left = rec(0, 1 * M, 0, 1 * M, "+", 1)
    right = rec(3 * M, 4 * M, 3 * M, 4 * M, "+", 2)                  # gap 2 Mb = 2x flank
    assert len(run(tmp_path, capsys, [left, right])) == 2
    assert len(run(tmp_path, capsys, [left, right], stitch_flank_factor=2.0)) == 3
    assert len(run(tmp_path, capsys, [left, right], stitch_flank_factor=0)) == 3


# ---------------------------------------------------------------- stitching under genepair

def gene_pairs(block, a0, b0, n=3, step=M, strand="+", ca="A_chr1", cb="B_chr1"):
    """A block cut into n consecutive gene-pair segments, as -partition genepair emits."""
    return [rec(a0 + i * step, a0 + (i + 1) * step, b0 + i * step, b0 + (i + 1) * step,
                strand, block, ca, cb) for i in range(n)]


def test_gap_between_blocks_is_stitched_under_genepair(tmp_path, capsys):
    # Two blocks of 3 Mb each (3 segments apiece) with a 1 Mb gap between them. The segments
    # bordering the gap are only 1 Mb, but the BLOCKS they belong to are 3 Mb, so the gap is
    # backed by synteny and must be stitched.
    lines = gene_pairs(1, 0, 0) + gene_pairs(2, 4 * M, 4 * M)
    out = spans(run(tmp_path, capsys, lines, threshold=0))
    assert (3.0, 4.0, 3.0, 4.0, "+") in out
    assert len(out) == 7


def test_genepair_stitch_uses_block_extent_not_segment_length(tmp_path, capsys):
    # Gap (2 Mb) is longer than the bordering segments (1 Mb) but shorter than the blocks
    # (3 Mb): allowed. Bumping the gap past the block extent is not.
    ok = gene_pairs(1, 0, 0) + gene_pairs(2, 5 * M, 5 * M)
    assert len(spans(run(tmp_path, capsys, ok, threshold=0))) == 7
    far = gene_pairs(1, 0, 0) + gene_pairs(2, 10 * M, 10 * M)      # 7 Mb gap vs 3 Mb blocks
    assert len(spans(run(tmp_path, capsys, far, threshold=0))) == 6


def test_genepair_stitch_respects_block_majority_strand(tmp_path, capsys):
    # Block 1 is 2 '+' segments and 1 '-' segment: majority '+', so it pairs with the '+'
    # block 2 and the synthetic record is '+'.
    mixed = gene_pairs(1, 0, 0)[:2] + [rec(2 * M, 3 * M, 2 * M, 3 * M, "-", 1)]
    out = spans(run(tmp_path, capsys, mixed + gene_pairs(2, 4 * M, 4 * M), threshold=0))
    assert (3.0, 4.0, 3.0, 4.0, "+") in out


def test_opposite_strand_block_between_genepair_blocks_blocks_the_stitch(tmp_path, capsys):
    inverted = gene_pairs(9, 3 * M, 3 * M, n=1, strand="-")
    lines = gene_pairs(1, 0, 0) + inverted + gene_pairs(2, 4 * M, 4 * M)
    out = spans(run(tmp_path, capsys, lines, threshold=0))
    # only the inverted record itself may occupy that span -- no '+' stitch over it
    assert [s for s in out if s[:4] == (3.0, 4.0, 3.0, 4.0)] == [(3.0, 4.0, 3.0, 4.0, "-")]


def test_records_without_a_block_id_still_stitch(tmp_path, capsys):
    # Legacy 3-column coords: no block id, so each record must count as its own block.
    p = tmp_path / "legacy.coords"
    p.write_text(f"A_chr1:0..{2 * M}\tB_chr1:0..{2 * M}\t+\n"
                 f"A_chr1:{3 * M}..{6 * M}\tB_chr1:{3 * M}..{6 * M}\t+\n")
    cons.process_file(str(p), 0, True)
    out = spans(capsys.readouterr().out.splitlines())
    assert (2.0, 3.0, 2.0, 3.0, "+") in out
