#!/usr/bin/env python3
"""Map window-local alignments back to genome coordinates and emit contract-compliant PAF.

Every aligner backend returns AlnRecord in the SUBMITTED sequence frame. This module is the
single place that maps those to parent coordinates and builds the tags, so the coordinate
arithmetic is written and tested once rather than once per backend.

The PAF contract is not negotiable and every way of breaking it is silent -- Steps 11-17
parse positionally. See step10_dev/memo_step10_diagnosis.md for the full audit. In short:

  * cg:Z: must use M, not =/X. The awk in Step 11 sums ([0-9]+)M; with =/X it appends 0,
    k2p becomes NA on every line, and weighted_paf.py leaves the .tsv outputs zero bytes
    without failing the run (the shell has no `set -e`).
  * cs:Z: must be SHORT form. Long form makes every substitution unparseable, k2p comes out
    -0.000000, and a sed in Step 11 rewrites that to 0.000000 -- i.e. perfect identity.
  * de:f: must land at 1-based column 21. paf_to_bed.py reads columns[20] positionally with
    no tag-name check, so the whole NM..s2 tag prefix has to be present and in order.

Semantics below were reverse-engineered from minimap2 2.30 output on this dataset and hold
across all sampled records:
    de       = (mismatch + gap_opens) / (identical + mismatch + gap_opens)   [%.4f]
    invariant: sum(M in cg) == sum(':' runs in cs) + count('*' in cs)
    invariant: count(I in cg) == count('+' in cs);  count(D) == count('-')
"""

import hashlib
from dataclasses import dataclass
from typing import List, Tuple

# Op vocabulary used by every backend. I = query has extra bases relative to target
# (minimap2 'I' / cs '+'); D = target has extra bases (minimap2 'D' / cs '-').
OPS = ("=", "X", "I", "D")

# Tag carrying the segment identity, appended as the FINAL field of every emitted row.
# It MUST stay last: this module's contract (see the docstring above) is that de:f: sits at
# 1-based column 21, because paf_to_bed.py reads columns[20] positionally with no tag-name
# check. Appending is safe; inserting anywhere earlier silently corrupts every divergence
# value in the run.
SEG_TAG = "sd:Z:"


def segment_id(coord1: str, coord2: str, strand: str) -> str:
    """Deterministic 12-hex-char ID for one coords line.

    blake2b rather than hash(): the builtin is salted per process, so a resume in a second
    invocation would compute different IDs and re-align everything.
    """
    key = f"{coord1}|{coord2}|{strand}".encode()
    return hashlib.blake2b(key, digest_size=6).hexdigest()


@dataclass
class Segment:
    """One resolved line of the coords file.

    Coordinates are 1-based inclusive, exactly as written in *.anchors.coords, matching
    `samtools faidx name:start-end`.
    """

    t_name: str
    t_start: int
    t_end: int
    q_name: str
    q_start: int
    q_end: int
    q_revcomp: bool  # whether the SUBMITTED query was reverse-complemented
    seg_id: str = ""  # paf_emit.segment_id() of the source coords line; "" when unknown

    @property
    def t_len(self) -> int:
        return self.t_end - self.t_start + 1

    @property
    def q_len(self) -> int:
        return self.q_end - self.q_start + 1


@dataclass
class AlnRecord:
    """One alignment in the SUBMITTED frame. Backends produce these and nothing else."""

    q_start: int  # 0-based half-open, submitted query frame
    q_end: int
    t_start: int  # 0-based half-open, submitted target frame
    t_end: int
    strand: str  # '+'/'-' of submitted query vs submitted target
    ops: List[Tuple[str, int]]


def map_to_parent(rec: AlnRecord, seg: Segment):
    """Map a submitted-frame alignment to parent genome coordinates.

    Returns (q_start, q_end, t_start, t_end, strand), all 0-based half-open.

    The window's parent offset is start-1, not start: extraction slices [start-1:end] for a
    1-based-inclusive coord. Adding `start` (as the original run_minimap did) is off by one.

    When the submitted query was reverse-complemented the query coordinates must be MIRRORED,
    not shifted. With S = revcomp(P) and P = parent[q_start-1 : q_end] of length L,
    S[i] = comp(P[L-1-i]), so [qs, qe) in S maps to [L-qe, L-qs) in P, i.e.
    [q_start-1 + L-qe, q_start-1 + L-qs) = [q_end - qe, q_end - qs) in parent 0-based.
    The reported strand flips too, since the submitted query was already flipped once.
    """
    t_off = seg.t_start - 1
    t_ps = rec.t_start + t_off
    t_pe = rec.t_end + t_off

    if seg.q_revcomp:
        q_ps = seg.q_end - rec.q_end
        q_pe = seg.q_end - rec.q_start
        strand = "-" if rec.strand == "+" else "+"
    else:
        q_off = seg.q_start - 1
        q_ps = rec.q_start + q_off
        q_pe = rec.q_end + q_off
        strand = rec.strand

    return q_ps, q_pe, t_ps, t_pe, strand


def build_cg(ops) -> str:
    """CIGAR in minimap2's vocabulary: = and X collapse into M (match-or-mismatch)."""
    out = []
    run = 0
    for op, n in ops:
        if op in "=X":
            run += n
        else:
            if run:
                out.append(f"{run}M")
                run = 0
            out.append(f"{n}{op}")
    if run:
        out.append(f"{run}M")
    return "".join(out)


def build_cs(ops, tseq: str, qseq: str, t_start: int, q_start: int) -> str:
    """Short-form cs:Z: from the ops and the SUBMITTED sequences.

    :N identical run, *ab substitution (target base then query base), +seq query insertion,
    -seq target deletion. All sequence lowercase, per minimap2.
    """
    ti, qi = t_start, q_start
    out = []
    run = 0
    for op, n in ops:
        if op == "=":
            run += n
            ti += n
            qi += n
            continue
        if run:
            out.append(f":{run}")
            run = 0
        if op == "X":
            for k in range(n):
                out.append(f"*{tseq[ti + k].lower()}{qseq[qi + k].lower()}")
            ti += n
            qi += n
        elif op == "I":
            out.append(f"+{qseq[qi:qi + n].lower()}")
            qi += n
        elif op == "D":
            out.append(f"-{tseq[ti:ti + n].lower()}")
            ti += n
        else:
            raise ValueError(f"unknown op {op!r}")
    if run:
        out.append(f":{run}")
    return "".join(out)


def gap_compressed_de(ops) -> float:
    """minimap2's de: each gap RUN counts as one event, not one per base."""
    ident = sum(n for op, n in ops if op == "=")
    mism = sum(n for op, n in ops if op == "X")
    gapo = sum(1 for op, n in ops if op in "ID")
    denom = ident + mism + gapo
    return (mism + gapo) / denom if denom else 0.0


def counts(ops):
    ident = sum(n for op, n in ops if op == "=")
    mism = sum(n for op, n in ops if op == "X")
    gapbp = sum(n for op, n in ops if op in "ID")
    return ident, mism, gapbp


def format_paf(rec: AlnRecord, seg: Segment, tseq: str, qseq: str) -> str:
    """Emit one PAF line with minimap2's exact column and tag layout.

    Cols 10-12 and the NM..s2 tag block are never read by Steps 11-17 (verified) but are
    positional padding that de:f: at col 21 depends on, so they are emitted faithfully.
    """
    q_ps, q_pe, t_ps, t_pe, strand = map_to_parent(rec, seg)
    ident, mism, gapbp = counts(rec.ops)
    blen = ident + mism + gapbp
    de = gap_compressed_de(rec.ops)
    nm = mism + gapbp

    cg = build_cg(rec.ops)
    cs = build_cs(rec.ops, tseq, qseq, rec.t_start, rec.q_start)

    cols = [
        seg.q_name, str(seg.q_len), str(q_ps), str(q_pe), strand,
        seg.t_name, str(seg.t_len), str(t_ps), str(t_pe),
        str(ident), str(blen), "60",
        f"NM:i:{nm}", f"ms:i:{ident}", f"AS:i:{ident}", "nn:i:0",
        "tp:A:P", "cm:i:0", f"s1:i:{ident}", "s2:i:0",
        f"de:f:{de:.4f}", "rl:i:0",
        f"cg:Z:{cg}", f"cs:Z:{cs}",
    ]
    if seg.seg_id:
        cols.append(f"{SEG_TAG}{seg.seg_id}")
    return "\t".join(cols)


def shift_native_paf(line: str, seg: Segment) -> str:
    """Fix coordinates on a native minimap2 PAF line, leaving every other field untouched.

    Used by the minimap2 backend so its output stays byte-identical apart from the
    coordinate correction -- which lets the regression gate assert that de/k2p are unmoved.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 12:
        raise ValueError(f"short PAF line: {line!r}")
    rec = AlnRecord(
        q_start=int(cols[2]), q_end=int(cols[3]),
        t_start=int(cols[7]), t_end=int(cols[8]),
        strand=cols[4], ops=[],
    )
    q_ps, q_pe, t_ps, t_pe, strand = map_to_parent(rec, seg)
    cols[2], cols[3], cols[4] = str(q_ps), str(q_pe), strand
    cols[7], cols[8] = str(t_ps), str(t_pe)
    if seg.seg_id:
        cols.append(f"{SEG_TAG}{seg.seg_id}")
    return "\t".join(cols)
