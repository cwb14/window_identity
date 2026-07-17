#!/usr/bin/env python3
"""Partition a syntenic segment into maximal colinear chains.

Why this exists: a syntenic block routinely contains an internal inversion -- in this data,
minimap2's records split roughly 50/50 by strand inside blocks labelled '+'. A colinear
end-to-end alignment path CANNOT represent an inversion. Forced across one, a global aligner
bridges it with a long run of compensating mismatches and micro-indels: a confident-looking
alignment of non-homologous sequence, with every indel boundary inside it meaningless. No
penalty setting fixes that; it is a structural mismatch between the model and the sequence.

So before any base-level alignment, the segment is cut at structural boundaries. A cheap
minimap2 pass supplies anchor records (it reports per-record strand and is already run for
orientation); those are grouped into maximal runs that a single colinear path CAN represent.
Each chain is then aligned end-to-end on its own and emitted as its own PAF record with its
own strand.

The result: end-to-end coverage WITHIN each colinear run -- so the intergenic, repeat-rich
sequence minimap2 skips is still recovered -- while inversions come out as separate records
with the correct strand instead of garbage. Fragmentation rises above 1.0, but for a
principled reason: it now counts structural segments, not aligner timidity.
"""

from dataclasses import dataclass
from typing import List, Optional


@dataclass
class Chain:
    """A maximal colinear run, in SUBMITTED-segment (window-local) coordinates."""

    t_start: int
    t_end: int
    q_start: int  # always forward-frame, q_start < q_end
    q_end: int
    strand: str
    n_anchors: int

    @property
    def t_len(self):
        return self.t_end - self.t_start

    @property
    def q_len(self):
        return self.q_end - self.q_start


@dataclass
class Anchor:
    t_start: int
    t_end: int
    q_start: int
    q_end: int
    strand: str


def parse_paf_anchors(paf_text: str, min_len: int = 200) -> List[Anchor]:
    """Window-local anchors from a raw minimap2 PAF (coordinates NOT yet shifted)."""
    out = []
    for line in paf_text.strip().split("\n"):
        if not line:
            continue
        c = line.split("\t")
        if len(c) < 12:
            continue
        t_start, t_end = int(c[7]), int(c[8])
        q_start, q_end = int(c[2]), int(c[3])
        if (t_end - t_start) < min_len or (q_end - q_start) < min_len:
            continue
        out.append(Anchor(t_start, t_end, q_start, q_end, c[4]))
    return out


def _colinear_with(chain_anchors: List[Anchor], nxt: Anchor, slack: int) -> bool:
    """Can one colinear path still cover chain_anchors + nxt?

    Anchors are walked in ascending target order. On '+' the query must also advance; on '-'
    the query must retreat (minimap2 reports minus-strand query coords in the forward frame,
    so a colinear inverted run has ascending target and DESCENDING query).
    """
    prev = chain_anchors[-1]
    if nxt.strand != prev.strand:
        return False
    if nxt.strand == "+":
        return nxt.q_start >= prev.q_start - slack
    return nxt.q_end <= prev.q_end + slack


def colinear_chains(anchors: List[Anchor], slack: int = 100, min_chain_len: int = 0,
                    t_len: Optional[int] = None,
                    q_len: Optional[int] = None) -> List[Chain]:
    """Greedy maximal colinear chaining over target-sorted anchors.

    A strand change, or a query coordinate that moves the wrong way by more than `slack`,
    ends the chain -- that is a structural boundary an alignment path must not cross.

    When t_len/q_len are supplied the chains are then extended to TILE the whole segment.
    That extension is not cosmetic: a chain's span is derived from anchors, and anchors
    rarely reach a segment's edges, so without it the flanks are never aligned at all.
    Measured on real Atha/Chis segments, anchor envelopes covered as little as 7.9% of the
    query, and every TE in a flank was silently lost.
    """
    if not anchors:
        return []
    ordered = sorted(anchors, key=lambda a: (a.t_start, a.t_end))

    groups: List[List[Anchor]] = []
    current: List[Anchor] = [ordered[0]]
    for a in ordered[1:]:
        if _colinear_with(current, a, slack):
            current.append(a)
        else:
            groups.append(current)
            current = [a]
    groups.append(current)

    chains = []
    for g in groups:
        c = Chain(
            t_start=min(a.t_start for a in g),
            t_end=max(a.t_end for a in g),
            q_start=min(a.q_start for a in g),
            q_end=max(a.q_end for a in g),
            strand=g[0].strand,
            n_anchors=len(g),
        )
        if c.t_len >= min_chain_len and c.q_len >= min_chain_len:
            chains.append(c)
    chains = _drop_overlaps(chains)
    if t_len is not None and q_len is not None:
        chains = _tile(chains, t_len, q_len)
    return chains


def _tile(chains: List[Chain], t_len: int, q_len: int) -> List[Chain]:
    """Extend chains so they partition the segment, leaving no base unaligned.

    Outer edges reach the segment boundaries. Between two chains the true structural
    breakpoint lies somewhere in the unanchored gap; absent better evidence the gap is split
    at its midpoint rather than handed arbitrarily to one side.

    With a single chain this makes the result identical to a plain global alignment, which is
    the point: chaining must cost nothing when there is no structure to resolve.
    """
    if not chains:
        return chains
    ordered = sorted(chains, key=lambda c: c.t_start)
    # Query intervals must be ordered too, or midpoints would cross. An inversion flips
    # orientation locally but does not move the block, so this normally holds; if it does
    # not, the segment is too rearranged to tile safely and is left as-is.
    if any(a.q_start > b.q_start for a, b in zip(ordered, ordered[1:])):
        return ordered

    for i, c in enumerate(ordered):
        if i == 0:
            c.t_start = 0
            c.q_start = 0
        if i == len(ordered) - 1:
            c.t_end = t_len
            c.q_end = q_len
        if i + 1 < len(ordered):
            nxt = ordered[i + 1]
            t_mid = (c.t_end + nxt.t_start) // 2
            q_mid = (c.q_end + nxt.q_start) // 2
            c.t_end = t_mid
            nxt.t_start = t_mid
            c.q_end = q_mid
            nxt.q_start = q_mid
    return ordered


def _drop_overlaps(chains: List[Chain]) -> List[Chain]:
    """Keep chains whose target AND query spans are mutually disjoint.

    Overlapping spans would align the same bases twice and double-count them in the
    M-sum that weights every downstream divergence average.
    """
    kept: List[Chain] = []
    for c in sorted(chains, key=lambda x: -(x.t_len + x.q_len)):
        clash = False
        for k in kept:
            if not (c.t_end <= k.t_start or c.t_start >= k.t_end):
                clash = True
                break
            if not (c.q_end <= k.q_start or c.q_start >= k.q_end):
                clash = True
                break
        if not clash:
            kept.append(c)
    return sorted(kept, key=lambda x: x.t_start)


def inversion_count(chains: List[Chain]) -> int:
    """Strand alternations -- how many structural boundaries were found."""
    if len(chains) < 2:
        return 0
    return sum(1 for a, b in zip(chains, chains[1:]) if a.strand != b.strand)
