#!/usr/bin/env python3
"""Anchor dotplot on a shared bp scale, with Ks colouring.

The same anchors jcvi's ortholog dotplot draws (one dot per anchor gene pair), placed at gene
midpoints in bp rather than gene rank:

  * both axes share one bp-per-inch scale (equal aspect), so a single scale bar reads
    chromosome sizes on either genome, and local expansion shows as slope;
  * chromosomes are labelled chr1..chrN, genome prefix stripped;
  * query (y) chromosomes are drawn reverse-complemented exactly when riparian.py would
    flip them -- same blocks, same decision code -- so the two plots always agree;
  * with --kaks, dots are coloured by Ks (viridis, square-root scale from 0 to the pair's
    95th percentile, higher values clipped to the top colour). Ks 'NA' is a real 0 (see
    ks_summary.py). Grey marks anchors without a usable Ks: never scored by ParaAT (e.g. an
    internal stop codon), or at or above --ks-max.

python dotplot.py --anchors PaA.PaB.clean.anchors --bed PaA.bed PaB.bed \
    --fai PaA_mod.fa.fai PaB_mod.fa.fai --coords PaA.PaB.anchors.coords \
    --kaks PaA.PaB.kaks.tsv -o PaA.PaB.dotplot
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import ks_summary  # noqa: E402
import riparian as rp  # noqa: E402  (also imports numpy + matplotlib, Agg backend)
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import PowerNorm  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar  # noqa: E402

tlog, die = rp.tlog, rp.die

KS_CMAP = "viridis"
KS_GAMMA = 0.5            # square-root scale: spreads the dense low-Ks end
KS_CLIP_PERCENTILE = 95   # top of the colour scale; an outlier cannot flatten the gradient
KS_EXTEND_FRAC = 0.05    # clipped-values arrow, drawn inside the colorbar's allotted height
NO_KS_COLOUR = "#b3b3b3"
DOT_COLOUR = "#303030"    # when there is no Ks to colour by
DOT_SIZE = 2.0            # points^2
DRAW_SEED = 0


# -----------------------------
# Inputs
# -----------------------------

def load_genomes(fai_paths):
    """x and y genome from their FAIs: name, called chromosomes, and seq -> genome."""
    order, chroms, seq2genome = [], {}, {}
    for path in fai_paths:
        if not os.path.exists(path):
            die(f"missing FAI: {path}")
        g = rp.genome_label_from_fai(path)
        chrom2len, all2len, _ = rp.select_chromosomes(path)
        if not chrom2len:
            die(f"{g} is not chromosome-level ({path}); cannot draw a dotplot.")
        order.append(g)
        chroms[g] = chrom2len
        for name in all2len:
            seq2genome[name] = g
    if len(set(order)) != 2:
        die(f"need two distinct genomes, got {order}")
    return {"order": order, "chroms": chroms, "seq2genome": seq2genome}


def flipped_chromosomes(coords_path, genomes):
    """y chromosomes riparian.py would draw reverse-complemented against the x genome."""
    gx, gy = genomes["order"]
    blocks = rp.parse_coords(coords_path, genomes["seq2genome"])
    kept, _ = rp.filter_pair_blocks(blocks, gx, gy, genomes["chroms"], rp.MIN_BLOCK_LEN,
                                    rp.MAX_LEN_RATIO, rp.RESCUE_UNCOVERED, rp.RESCUE_RULE)
    tx, ty = rp.Track(gx, genomes["chroms"][gx]), rp.Track(gy, genomes["chroms"][gy])
    rp.optimise_track(ty, tx, kept, do_order=False)
    return {c for c in ty.order if ty.flip[c]}


def read_gene_positions(bed_paths):
    """gene id -> (sequence, midpoint bp)."""
    genes = {}
    for path in bed_paths:
        if not os.path.exists(path):
            die(f"missing BED: {path}")
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) < 4 or line.startswith(("#", "track", "browser")):
                    continue
                genes[f[3].strip()] = (f[0], (int(f[1]) + int(f[2])) / 2)
    return genes


def read_anchors(path):
    """(gene in x genome, gene in y genome) per anchor line; '###' block separators skipped."""
    if not os.path.exists(path):
        die(f"missing anchors: {path}")
    pairs = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and not f[0].startswith("#"):
                pairs.append((f[0], f[1]))
    return pairs


def read_ks_by_pair(kaks_path, max_ks):
    """'geneX-geneY' -> Ks for estimable pairs with 0 <= Ks < max_ks."""
    if not os.path.exists(kaks_path):
        die(f"missing Ks table: {kaks_path}")
    return {name: ks for name, ks in ks_summary.iter_ks(kaks_path)
            if ks is not None and 0.0 <= ks < max_ks}


# -----------------------------
# Geometry
# -----------------------------

def axis_offsets(chrom2len):
    """chromosome -> bp offset of its start along the axis, chromosomes end to end."""
    offsets, pos = {}, 0
    for c, L in chrom2len.items():
        offsets[c] = pos
        pos += L
    return offsets, pos


def dot_positions(anchors, genes, genomes, flips):
    """Axis positions (bp) of every anchor whose genes lie on called chromosomes."""
    gx, gy = genomes["order"]
    cx, cy = genomes["chroms"][gx], genomes["chroms"][gy]
    ox, _ = axis_offsets(cx)
    oy, _ = axis_offsets(cy)
    dots = {"pair": [], "x": [], "y": [], "skipped": 0}
    for a, b in anchors:
        if a not in genes or b not in genes or genes[a][0] not in cx or genes[b][0] not in cy:
            dots["skipped"] += 1
            continue
        (sa, pa), (sb, pb) = genes[a], genes[b]
        dots["pair"].append(f"{a}-{b}")
        dots["x"].append(ox[sa] + pa)
        dots["y"].append(oy[sb] + (cy[sb] - pb if sb in flips else pb))
    return dots


def ks_scale(values):
    """(vmax, clipped): the colour scale's top, and whether any value exceeds it."""
    vals = np.asarray(values, dtype=float)
    vmax = float(np.percentile(vals, KS_CLIP_PERCENTILE)) if vals.size else 0.0
    if vmax <= 0:
        pos = vals[vals > 0]
        vmax = float(pos.max()) if pos.size else 1e-3
    return vmax, bool((vals > vmax).any())


# -----------------------------
# Drawing
# -----------------------------

def _short(chrom):
    return re.sub(r"^.*?_", "", chrom)


def draw(dots, genomes, ks, width):
    """Figure and axes. ks is {'geneX-geneY': Ks} to colour by, or None."""
    gx, gy = genomes["order"]
    cx, cy = genomes["chroms"][gx], genomes["chroms"][gy]
    ox, span_x = axis_offsets(cx)
    oy, span_y = axis_offsets(cy)
    x = np.asarray(dots["x"]) / 1e6
    y = np.asarray(dots["y"]) / 1e6

    fig, ax = plt.subplots(figsize=(width, max(2.0, width * span_y / span_x)))
    ax.set_xlim(0, span_x / 1e6)
    ax.set_ylim(0, span_y / 1e6)
    ax.set_aspect("equal")

    for c in list(cx)[1:]:
        ax.axvline(ox[c] / 1e6, color="#d9d9d9", linewidth=0.5, zorder=0)
    for c in list(cy)[1:]:
        ax.axhline(oy[c] / 1e6, color="#d9d9d9", linewidth=0.5, zorder=0)
    ax.set_xticks([(ox[c] + L / 2) / 1e6 for c, L in cx.items()], [_short(c) for c in cx])
    ax.set_yticks([(oy[c] + L / 2) / 1e6 for c, L in cy.items()], [_short(c) for c in cy])
    ax.tick_params(length=0, labelsize=7)
    ax.set_xlabel(gx, fontsize=9, fontstyle="italic", fontweight="bold")
    ax.set_ylabel(gy, fontsize=9, fontstyle="italic", fontweight="bold")
    for s in ax.spines.values():
        s.set_linewidth(0.6)

    if ks is None:
        ax.scatter(x, y, s=DOT_SIZE, c=DOT_COLOUR, linewidths=0, rasterized=True, gid="dots")
    else:
        vals = np.array([ks.get(p, np.nan) for p in dots["pair"]])
        known = ~np.isnan(vals)
        ax.scatter(x[~known], y[~known], s=DOT_SIZE, c=NO_KS_COLOUR, linewidths=0,
                   rasterized=True, gid="no_ks")
        # Random draw order: dots overlap heavily at genome scale, so drawing sorted by Ks
        # would let whichever end is drawn last paint over the rest.
        idx = np.flatnonzero(known)
        idx = idx[np.random.default_rng(DRAW_SEED).permutation(idx.size)]
        vmax, clipped = ks_scale(vals[known])
        sc = ax.scatter(x[idx], y[idx], s=DOT_SIZE, c=vals[idx], cmap=KS_CMAP,
                        norm=PowerNorm(KS_GAMMA, vmin=0.0, vmax=vmax), linewidths=0,
                        rasterized=True, gid="ks")
        # Carved from the plot's own box so it tracks the equal-aspect height; fig.colorbar(ax=)
        # sizes the bar before the aspect shrinks the plot.
        cax = make_axes_locatable(ax).append_axes("right", size="2.5%", pad=0.08)
        cb = fig.colorbar(sc, cax=cax, extend="max" if clipped else "neither",
                          extendfrac=KS_EXTEND_FRAC)
        cb.set_label("Ks", fontsize=8)
        cb.ax.tick_params(labelsize=6)
        n_grey = int((~known).sum())
        if n_grey:
            cb.ax.text(0.5, -0.04, f"grey: Ks n/a\n(n = {n_grey:,})", transform=cb.ax.transAxes,
                       ha="center", va="top", fontsize=6)

    frac, label = rp.nice_scale_bar(span_x)
    # A fixed offset in points below the plot, level with the axis label, whatever the
    # plot's height (an axes-fraction offset drifts away on tall plots).
    bar = AnchoredSizeBar(ax.transData, frac * span_x / 1e6, label, loc="upper left",
                          bbox_to_anchor=(0.0, 0.0),
                          bbox_transform=offset_copy(ax.transAxes, fig=fig, y=-20,
                                                     units="points"),
                          frameon=False, borderpad=0, sep=2, size_vertical=0,
                          fontproperties={"size": 7})
    bar.set_gid("scale_bar")
    ax.add_artist(bar)
    return fig, ax


def scale_bar_label(ax):
    for artist in ax.artists:
        if artist.get_gid() == "scale_bar":
            return artist.txt_label.get_text()
    return None


def dot_layers(ax):
    """(Ks-coloured dots, grey dots) collections, as drawn by draw() with Ks."""
    by_gid = {c.get_gid(): c for c in ax.collections}
    return by_gid.get("ks"), by_gid.get("no_ks")


# -----------------------------
# Main
# -----------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Anchor dotplot on a shared bp scale: chromosome labels, scale bar, "
                    "riparian-consistent orientation flips, optional Ks colouring.")
    p.add_argument("--anchors", required=True,
                   help="anchor gene pairs, column 1 in the x genome "
                        "({ID1}.{ID2}.clean.anchors)")
    p.add_argument("--bed", nargs=2, required=True, metavar=("X_BED", "Y_BED"),
                   help="gene BEDs for the x and y genomes (gene id in column 4)")
    p.add_argument("--fai", nargs=2, required=True, metavar=("X_FAI", "Y_FAI"),
                   help="FAI per genome, x first; genome names come from the file names "
                        "('PaA_mod.fa.fai' -> PaA)")
    p.add_argument("--coords", help="syntenic block coords ({ID1}.{ID2}.anchors.coords); "
                                    "required unless --no-flip")
    p.add_argument("--kaks",
                   help="KaKs_Calculator table ({ID1}.{ID2}.kaks.tsv); colours dots by Ks")
    p.add_argument("--ks-max", type=float, default=2.0,
                   help="saturation cutoff: Ks >= this is drawn grey, like an anchor ParaAT "
                        "never scored (default: 2.0)")
    p.add_argument("--no-flip", action="store_true",
                   help="draw y chromosomes in their native orientation")
    p.add_argument("-o", "--out-prefix", default="dotplot", help="output path prefix")
    p.add_argument("--formats", default="pdf,png",
                   help="comma-separated: pdf,png (default: pdf,png)")
    p.add_argument("--width", type=float, default=8.0, help="figure width, inches (default: 8)")
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution (default: 300)")
    p.add_argument("-v", "--verbose", action="store_true", help="per-step counts")
    a = p.parse_args(argv)

    if not a.no_flip and not a.coords:
        die("--coords is needed to decide orientation flips (or pass --no-flip)")
    fmts = [f.strip().lower() for f in a.formats.split(",") if f.strip()]
    if not fmts or set(fmts) - {"pdf", "png"}:
        die(f"--formats must list pdf and/or png (got '{a.formats}')")

    genomes = load_genomes(a.fai)
    gx, gy = genomes["order"]
    flips = set() if a.no_flip else flipped_chromosomes(a.coords, genomes)
    genes = read_gene_positions(a.bed)
    dots = dot_positions(read_anchors(a.anchors), genes, genomes, flips)
    if not dots["pair"]:
        die(f"no anchors of {a.anchors} fall on chromosomes of {gx} and {gy}")
    ks = read_ks_by_pair(a.kaks, a.ks_max) if a.kaks else None

    tlog(f"Dotplot {gx} vs {gy}: {len(dots['pair']):,} anchors"
         + (f", {gy} flipped: {', '.join(_short(c) for c in rp.sort_chrom_names(flips))}"
            if flips else ""))
    tlog(f"  skipped {dots['skipped']:,} anchors off the called chromosomes", a.verbose)
    if ks is not None:
        vals = [ks[q] for q in dots["pair"] if q in ks]
        vmax, clipped = ks_scale(vals)
        tlog(f"  Ks: {len(vals):,} coloured ({sum(v == 0 for v in vals):,} with Ks = 0), "
             f"{len(dots['pair']) - len(vals):,} grey; colour scale 0..{vmax:.4g}"
             + (" (clipped)" if clipped else ""), a.verbose)

    fig, _ = draw(dots, genomes, ks, a.width)
    for f in fmts:
        path = f"{a.out_prefix}.{f}"
        fig.savefig(path, dpi=a.dpi, bbox_inches="tight", facecolor="white")
        tlog(f"  wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
