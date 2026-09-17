"""riparian.py chromosome flipping: every homologous pair should draw untwisted.

The old flip rule judged orientation from the covariance of block positions, which needs at
least two blocks. A chromosome carried by one whole-chromosome inverted block (PaA/PiA chr2,
chr4, chr6) could not be judged, kept its native orientation, and drew as a twist, while
multi-block chromosomes (PiA chr5) were flipped: inconsistent display of the same situation.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import riparian as rp  # noqa: E402

M = 1_000_000


def track(genome, lens):
    t = rp.Track(genome, dict(lens))
    t.relayout(0)
    t.denom = t.span or 1.0
    return t


def block(ca, sa, ea, cb, sb, eb, strand):
    return {"ga": "U", "ca": ca, "sa": sa, "ea": ea,
            "gb": "D", "cb": cb, "sb": sb, "eb": eb, "strand": strand}


def twisted(r):
    return (r["a2"] - r["a1"]) * (r["b2"] - r["b1"]) < 0


def test_single_inverted_block_is_flipped():
    up, dn = track("U", {"U_chr2": 100 * M}), track("D", {"D_chr2": 100 * M})
    b = block("U_chr2", 0, 100 * M, "D_chr2", 0, 100 * M, "-")
    rp.optimise_track(dn, up, [b], do_order=False)
    assert dn.flip["D_chr2"] is True


def test_collinear_block_is_not_flipped():
    up, dn = track("U", {"U_chr1": 100 * M}), track("D", {"D_chr1": 100 * M})
    b = block("U_chr1", 0, 100 * M, "D_chr1", 0, 100 * M, "+")
    rp.optimise_track(dn, up, [b], do_order=False)
    assert dn.flip["D_chr1"] is False


def test_flip_votes_by_bp_not_block_count():
    # Two small '+' blocks, one large '-' block: bp majority is '-'.
    up, dn = track("U", {"U_chr3": 100 * M}), track("D", {"D_chr3": 100 * M})
    blocks = [block("U_chr3", 0, 60 * M, "D_chr3", 0, 60 * M, "-"),
              block("U_chr3", 60 * M, 80 * M, "D_chr3", 60 * M, 80 * M, "+"),
              block("U_chr3", 80 * M, 95 * M, "D_chr3", 80 * M, 95 * M, "+")]
    rp.optimise_track(dn, up, blocks, do_order=False)
    assert dn.flip["D_chr3"] is True


def test_flip_accounts_for_a_flipped_partner_above():
    # The partner is already displayed reversed, so a '+' block only reads straight
    # if this chromosome is reversed too (3+ genome stacks).
    up, dn = track("U", {"U_chr1": 100 * M}), track("D", {"D_chr1": 100 * M})
    up.flip["U_chr1"] = True
    b = block("U_chr1", 0, 100 * M, "D_chr1", 0, 100 * M, "+")
    rp.optimise_track(dn, up, [b], do_order=False)
    assert dn.flip["D_chr1"] is True


def test_flip_judged_against_dominant_partner_only():
    # 80 Mb '+' to U_chr1 dominates a 10 Mb '-' translocation to U_chr2.
    up = track("U", {"U_chr1": 100 * M, "U_chr2": 100 * M})
    dn = track("D", {"D_chr1": 100 * M})
    blocks = [block("U_chr1", 0, 80 * M, "D_chr1", 0, 80 * M, "+"),
              block("U_chr2", 0, 10 * M, "D_chr1", 85 * M, 95 * M, "-")]
    rp.optimise_track(dn, up, blocks, do_order=False)
    assert dn.flip["D_chr1"] is False


def test_no_flip_keeps_native_orientation():
    up, dn = track("U", {"U_chr2": 100 * M}), track("D", {"D_chr2": 100 * M})
    b = block("U_chr2", 0, 100 * M, "D_chr2", 0, 100 * M, "-")
    rp.optimise_track(dn, up, [b], do_order=False, do_flip=False)
    assert dn.flip["D_chr2"] is False


def test_every_homolog_pair_draws_untwisted():
    up = track("U", {"U_chr1": 100 * M, "U_chr2": 80 * M})
    dn = track("D", {"D_chr1": 90 * M, "D_chr2": 70 * M})
    blocks = [block("U_chr1", 0, 100 * M, "D_chr1", 0, 90 * M, "+"),
              block("U_chr2", 0, 80 * M, "D_chr2", 0, 70 * M, "-")]
    rp.optimise_track(dn, up, blocks, do_order=False)
    dn.relayout(0)
    scene = rp.build_scene(["U", "D"], {"U": up, "D": dn}, {("U", "D"): blocks},
                           rp.build_palette(list(up.lens)), "U", "bp", 1.0)
    assert not any(twisted(r) for r in scene["ribbons"])


def test_inversion_within_a_chromosome_still_twists():
    # Flipping is per chromosome: a minority inverted block inside a collinear
    # chromosome is real structure and must stay visible.
    up, dn = track("U", {"U_chr1": 100 * M}), track("D", {"D_chr1": 100 * M})
    blocks = [block("U_chr1", 0, 70 * M, "D_chr1", 0, 70 * M, "+"),
              block("U_chr1", 70 * M, 90 * M, "D_chr1", 70 * M, 90 * M, "-")]
    rp.optimise_track(dn, up, blocks, do_order=False)
    dn.relayout(0)
    scene = rp.build_scene(["U", "D"], {"U": up, "D": dn}, {("U", "D"): blocks},
                           rp.build_palette(list(up.lens)), "U", "bp", 1.0)
    by_strand = {r["strand"]: twisted(r) for r in scene["ribbons"]}
    assert by_strand == {"+": False, "-": True}


# ---------------------------------------------------------------- orientation report

def two_pair_scene():
    """U_chr1/D_chr1 collinear; U_chr2/D_chr2 mostly inverted (60 Mb '-', 20 Mb '+')."""
    up = track("U", {"U_chr1": 100 * M, "U_chr2": 80 * M})
    dn = track("D", {"D_chr1": 90 * M, "D_chr2": 80 * M})
    blocks = [block("U_chr1", 0, 100 * M, "D_chr1", 0, 90 * M, "+"),
              block("U_chr2", 0, 60 * M, "D_chr2", 0, 60 * M, "-"),
              block("U_chr2", 60 * M, 80 * M, "D_chr2", 60 * M, 80 * M, "+")]
    return up, dn, blocks


def test_relative_orientation_summarises_each_lower_chromosome():
    _, _, blocks = two_pair_scene()
    blocks.append(block("U_chr1", 0, 5 * M, "D_chr2", 80 * M, 85 * M, "-"))  # minor partner
    calls = rp.relative_orientation(blocks)
    assert calls["D_chr1"] == {"partner": "U_chr1", "bp": 90 * M, "opposite_bp": 0}
    assert calls["D_chr2"] == {"partner": "U_chr2", "bp": 80 * M, "opposite_bp": 60 * M}


def test_orientation_table_rows(tmp_path):
    up, dn, blocks = two_pair_scene()
    rp.optimise_track(dn, up, blocks, do_order=False)
    out = tmp_path / "riparian.orientation.tsv"
    rp.write_orientation(str(out), ["U", "D"], {"U": up, "D": dn}, {("U", "D"): blocks})
    rows = [ln.split("\t") for ln in out.read_text().splitlines() if not ln.startswith("#")]
    assert rows[0] == ["upper", "lower", "orientation", "opposite_frac", "flipped_in_plot"]
    assert rows[1:] == [["U_chr1", "D_chr1", "+", "0.00", "no"],
                        ["U_chr2", "D_chr2", "-", "0.75", "yes"]]


def test_plot_labels_carry_no_prime(tmp_path, monkeypatch):
    from matplotlib.axes import Axes
    up, dn, blocks = two_pair_scene()
    rp.optimise_track(dn, up, blocks, do_order=False)
    dn.relayout(0)
    assert dn.flip["D_chr2"]
    scene = rp.build_scene(["U", "D"], {"U": up, "D": dn}, {("U", "D"): blocks},
                           rp.build_palette(list(up.lens)), "U", "bp",
                           max(up.span, dn.span))    # as main() passes; sizes the scale bar

    rp.render_html(scene, str(tmp_path / "p.html"), 0.6, 0.6)
    page = (tmp_path / "p.html").read_text()
    assert ">chr2</text>" in page and "′" not in page

    # Labels are drawn with a stroke path effect, which renders text as outlines in every
    # vector format, so record the strings handed to Axes.text instead.
    drawn, real_text = [], Axes.text
    monkeypatch.setattr(Axes, "text",
                        lambda self, x, y, s, *a, **k: (drawn.append(s), real_text(self, x, y, s, *a, **k))[1])
    rp.render_mpl(scene, [str(tmp_path / "p.png")], 11, 3.8, 72, 0.6, 0.6)
    assert "chr2" in drawn and not any("′" in s for s in drawn)


def test_cli_writes_orientation_table(tmp_path):
    (tmp_path / "U_mod.fa.fai").write_text(f"U_chr1\t{100 * M}\t0\t60\t61\nU_chr2\t{80 * M}\t0\t60\t61\n")
    (tmp_path / "D_mod.fa.fai").write_text(f"D_chr1\t{90 * M}\t0\t60\t61\nD_chr2\t{80 * M}\t0\t60\t61\n")
    (tmp_path / "U.D.anchors.coords").write_text(
        f"U_chr1:0..{100 * M}\tD_chr1:0..{90 * M}\t+\n"
        f"U_chr2:0..{80 * M}\tD_chr2:0..{80 * M}\t-\n")
    prefix = tmp_path / "riparian"
    rp.main(["--coords", str(tmp_path / "U.D.anchors.coords"),
             "--fai", str(tmp_path / "U_mod.fa.fai"), str(tmp_path / "D_mod.fa.fai"),
             "--order", "U,D", "--min-chrom-len", "1", "--formats", "html",
             "-o", str(prefix)])
    table = (tmp_path / "riparian.orientation.tsv").read_text()
    assert "U_chr2\tD_chr2\t-\t1.00\tyes" in table


# ---------------------------------------------------------------- shared block filtering

def test_filter_pair_blocks_transposes_and_filters():
    chroms = {"U": {"U_chr1": 100 * M}, "D": {"D_chr1": 100 * M}}
    raw = [
        # written D-first: must come back transposed so 'a' is the upper genome
        {"ga": "D", "ca": "D_chr1", "sa": 0, "ea": 50 * M,
         "gb": "U", "cb": "U_chr1", "sb": 0, "eb": 50 * M, "strand": "-"},
        # on a sequence that is not a chromosome
        {"ga": "U", "ca": "U_scaf9", "sa": 0, "ea": M,
         "gb": "D", "cb": "D_chr1", "sb": 0, "eb": M, "strand": "+"},
        # shorter than min_block_len on both sides
        {"ga": "U", "ca": "U_chr1", "sa": 60 * M, "ea": 60 * M + 5000,
         "gb": "D", "cb": "D_chr1", "sb": 60 * M, "eb": 60 * M + 5000, "strand": "+"},
    ]
    kept, stats = rp.filter_pair_blocks(raw, "U", "D", chroms, 10000, 5.0, 0.5, "both")
    assert kept == [{"ga": "U", "ca": "U_chr1", "sa": 0, "ea": 50 * M,
                     "gb": "D", "cb": "D_chr1", "sb": 0, "eb": 50 * M, "strand": "-"}]
    assert stats == {"raw": 3, "off_chrom": 1, "small": 1, "clean": 1, "rescued": 0,
                     "dropped": 0}
