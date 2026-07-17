#!/usr/bin/env python3
"""Determine the true relative orientation of two syntenic segments from sequence.

Why this exists: the strand column in *.anchors.coords cannot be the last word. The old
gene-order rule mislabelled ~21% of spans (fixed in gene_coords_extractor_all4.py, now
~97-100% correct against k-mer truth), but even the fixed rule cannot resolve a span that
straddles an inversion boundary, and 2-3% remain wrong on the easy pair.

That residue matters because minimap2 and a strand-sensitive aligner fail differently.
minimap2 searches both strands, so a mislabelled segment costs it nothing. WFA aligns the
literal strings: given the wrong orientation it returns an all-gap path, or -- with a
heuristic enabled -- manufactures a plausible ~70% identity out of pure noise. That is what
sank the previous wavefront attempt. So for strand-sensitive backends we settle orientation
from the sequence and treat the strand column as advisory.

Method: a cheap minimap2 pass, reading only the best record's strand.

Not k-mer sketching, which was the obvious idea and does not survive contact with the data:
exact 15-mer survival is p^15, so at Atha/Chis divergence (~75% identity) only ~1.3% of
15-mers are shared. Sub-sampling to bound memory on 8Mb segments then leaves single-digit
k-mer counts and the test stops separating orientations exactly where it is needed most.
A smaller k trades that for background noise. minimap2 already solves this problem well
across low, medium and high divergence, and is a dependency either way.

LAST does not need this module -- lastal searches both strands natively, so orientation
comes out of the alignment for free.
"""

import os
import subprocess
import tempfile
from typing import Optional

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def orientation(seq1: str, seq2: str, minimap2_bin: str = "minimap2",
                preset: str = "asm20", threads: int = 1) -> Optional[bool]:
    """True if seq2 must be reverse-complemented to align with seq1.

    Returns None when minimap2 finds no alignment -- the segment is not confidently
    homologous (often a consolidator-stitched synthetic span) and should be skipped
    rather than aligned on a guess.
    """
    d = tempfile.mkdtemp(prefix="orient_")
    try:
        t_fa = os.path.join(d, "t.fa")
        q_fa = os.path.join(d, "q.fa")
        with open(t_fa, "w") as fh:
            fh.write(f">t\n{seq1}\n")
        with open(q_fa, "w") as fh:
            fh.write(f">q\n{seq2}\n")
        proc = subprocess.run(
            [minimap2_bin, "-t", str(threads), "--secondary=no", "-x", preset,
             t_fa, q_fa],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            return None
        best = None
        best_matches = -1
        for line in proc.stdout.strip().split("\n"):
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            matches = int(cols[9])
            if matches > best_matches:
                best_matches = matches
                best = cols[4]
        if best is None:
            return None
        return best == "-"
    finally:
        for f in ("t.fa", "q.fa"):
            try:
                os.remove(os.path.join(d, f))
            except OSError:
                pass
        try:
            os.rmdir(d)
        except OSError:
            pass
