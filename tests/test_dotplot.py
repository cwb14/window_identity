"""dotplot.py: bp-scaled anchor dotplot with chromosome labels, scale bar, orientation flip
and optional Ks colouring."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import dotplot as dp  # noqa: E402
import ks_summary as ksum  # noqa: E402

M = 1_000_000


@pytest.fixture
def run(tmp_path):
    """U (chr1 100 Mb, chr2 50 Mb) vs D (chr1 90 Mb, chr2 60 Mb); D_chr2 is inverted."""
    (tmp_path / "U_mod.fa.fai").write_text(f"U_chr1\t{100 * M}\t0\t60\t61\nU_chr2\t{50 * M}\t0\t60\t61\n")
    (tmp_path / "D_mod.fa.fai").write_text(f"D_chr1\t{90 * M}\t0\t60\t61\nD_chr2\t{60 * M}\t0\t60\t61\n")
    (tmp_path / "U.bed").write_text(
        f"U_chr1\t{10 * M - 500}\t{10 * M + 500}\tu1\t0\t+\n"
        f"U_chr1\t{20 * M - 500}\t{20 * M + 500}\tu2\t0\t+\n"
        f"U_chr2\t{10 * M - 500}\t{10 * M + 500}\tu3\t0\t+\n")
    (tmp_path / "D.bed").write_text(
        f"D_chr1\t{10 * M - 500}\t{10 * M + 500}\td1\t0\t+\n"
        f"D_chr1\t{20 * M - 500}\t{20 * M + 500}\td2\t0\t+\n"
        f"D_chr2\t{50 * M - 500}\t{50 * M + 500}\td3\t0\t-\n")
    (tmp_path / "U.D.clean.anchors").write_text("###\nu1\td1\nu2\td2\n###\nu3\td3\n")
    (tmp_path / "U.D.anchors.coords").write_text(
        f"U_chr1:0..{100 * M}\tD_chr1:0..{90 * M}\t+\t1\n"
        f"U_chr2:0..{50 * M}\tD_chr2:0..{60 * M}\t-\t2\n")

    def kaks(name, ks, syn):
        r = dict.fromkeys(ksum.KAKS_COLUMNS, "NA")
        r.update({"Sequence": name, "Method": "YN", "Ks": ks, "Syn-Subs": syn})
        return "\t".join(r[c] for c in ksum.KAKS_COLUMNS)
    (tmp_path / "U.D.kaks.tsv").write_text("\n".join([
        "\t".join(ksum.KAKS_COLUMNS),
        kaks("u1-d1", "0.1", "20"),
        kaks("u2-d2", "NA", "NA"),      # no synonymous change -> Ks 0
        # u3-d3 absent: ParaAT never scored it (e.g. internal stop) -> grey
    ]) + "\n")
    return tmp_path


def base_args(run):
    return ["--anchors", str(run / "U.D.clean.anchors"),
            "--bed", str(run / "U.bed"), str(run / "D.bed"),
            "--fai", str(run / "U_mod.fa.fai"), str(run / "D_mod.fa.fai")]


def load(run, flip=True):
    genomes = dp.load_genomes([str(run / "U_mod.fa.fai"), str(run / "D_mod.fa.fai")])
    flips = (dp.flipped_chromosomes(str(run / "U.D.anchors.coords"), genomes) if flip
             else set())
    return genomes, flips


def test_inverted_query_chromosome_is_flipped(run):
    _, flips = load(run)
    assert flips == {"D_chr2"}


def test_dot_positions_use_chromosome_offsets_and_flips(run):
    genomes, flips = load(run)
    genes = dp.read_gene_positions([str(run / "U.bed"), str(run / "D.bed")])
    anchors = dp.read_anchors(str(run / "U.D.clean.anchors"))
    dots = dp.dot_positions(anchors, genes, genomes, flips)
    assert dots["pair"] == ["u1-d1", "u2-d2", "u3-d3"]
    assert dots["x"][2] == 100 * M + 10 * M              # U_chr2 sits after U_chr1
    assert dots["y"][2] == 90 * M + (60 * M - 50 * M)    # D_chr2 reversed for display


def test_ks_na_is_zero_and_unscored_pair_is_missing(run):
    ks = dp.read_ks_by_pair(str(run / "U.D.kaks.tsv"), max_ks=2.0)
    assert ks == {"u1-d1": 0.1, "u2-d2": 0.0}


def test_ks_scale_clips_outliers_at_95th_percentile():
    vals = np.append(np.linspace(0, 0.3, 100), 4.2)
    vmax, clipped = dp.ks_scale(vals)
    assert vmax == pytest.approx(np.percentile(vals, 95)) and clipped
    assert dp.ks_scale(np.zeros(10)) == (pytest.approx(1e-3), False)   # all identical


def test_figure_labels_scale_bar_and_ks_colours(run):
    genomes, flips = load(run)
    genes = dp.read_gene_positions([str(run / "U.bed"), str(run / "D.bed")])
    dots = dp.dot_positions(dp.read_anchors(str(run / "U.D.clean.anchors")), genes, genomes,
                            flips)
    ks = dp.read_ks_by_pair(str(run / "U.D.kaks.tsv"), max_ks=2.0)
    fig, ax = dp.draw(dots, genomes, ks, width=6)
    assert [t.get_text() for t in ax.get_xticklabels()] == ["chr1", "chr2"]
    assert [t.get_text() for t in ax.get_yticklabels()] == ["chr1", "chr2"]
    assert dp.scale_bar_label(ax) == "20 Mb"                  # ~1/10 of 150 Mb, rounded
    coloured, grey = dp.dot_layers(ax)
    assert sorted(coloured.get_array()) == [0.0, 0.1] and len(grey.get_offsets()) == 1
    assert len(fig.axes) == 2                                 # plot + colorbar

    fig, ax = dp.draw(dots, genomes, None, width=6)           # -kaks no: no colorbar
    assert len(fig.axes) == 1


def test_cli_writes_pdf_and_png(run):
    prefix = run / "U.D.dotplot"
    dp.main(base_args(run) + ["--coords", str(run / "U.D.anchors.coords"),
                              "--kaks", str(run / "U.D.kaks.tsv"), "-o", str(prefix)])
    assert (run / "U.D.dotplot.pdf").stat().st_size > 0
    assert (run / "U.D.dotplot.png").stat().st_size > 0


def test_cli_needs_coords_to_flip(run, capsys):
    with pytest.raises(SystemExit) as e:
        dp.main(base_args(run) + ["-o", str(run / "x")])
    assert e.value.code == 1 and "--coords" in capsys.readouterr().err
    dp.main(base_args(run) + ["--no-flip", "-o", str(run / "x"), "--formats", "png"])
    assert (run / "x.png").exists()


def test_colorbar_matches_plot_height(run):
    # Equal aspect shrinks the plot box after fig.colorbar() sizes the bar, so on a
    # square-ish genome pair the colorbar came out ~2/3 of the plot height.
    genomes, flips = load(run)
    genes = dp.read_gene_positions([str(run / "U.bed"), str(run / "D.bed")])
    dots = dp.dot_positions(dp.read_anchors(str(run / "U.D.clean.anchors")), genes, genomes,
                            flips)
    ks = dp.read_ks_by_pair(str(run / "U.D.kaks.tsv"), max_ks=2.0)
    fig, ax = dp.draw(dots, genomes, ks, width=6)
    fig.canvas.draw()
    cax = [a for a in fig.axes if a is not ax][0]
    plot, bar = ax.get_position(), cax.get_position()
    assert bar.y0 == pytest.approx(plot.y0)
    # Ks here is clipped, so matplotlib draws the extend arrow inside the allotted height:
    # bar + arrow span the plot exactly.
    assert bar.height * (1 + dp.KS_EXTEND_FRAC) == pytest.approx(plot.height, rel=1e-3)
