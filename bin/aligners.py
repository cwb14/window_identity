#!/usr/bin/env python3
"""Alignment backends. Each returns AlnRecord objects in the SUBMITTED sequence frame.

Backends align and nothing else -- no coordinate arithmetic, no tag building. That all lives
in paf_emit.py so it is written and tested once. See step10_dev/memo_step10_diagnosis.md.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import List

import chain
from paf_emit import AlnRecord

_CIGAR_RE = re.compile(r"(\d+)([=XIDM])")

# WFA2-lib's I/D are INVERTED relative to minimap2/PAF.
#
# WFAligner.alignEnd2End(pattern, text) with pattern=target, text=query emits a CIGAR where
# 'I' consumes the PATTERN and 'D' consumes the TEXT. Verified directly: a query missing
# 500bp yields '1000=500I1500=', and a query with 2000 extra bases yields '1201=2000D1799='.
# PAF is the other way round -- 'I' means the query has extra bases, 'D' means the target
# does. Left unswapped this silently inverts every +/- in cs:Z: and mis-sizes both spans.
_WFA_OP = {"=": "=", "X": "X", "I": "D", "D": "I", "M": "="}


def parse_wfa_cigar(cigar: str) -> List[tuple]:
    """WFA CIGAR -> PAF-convention ops, with I/D swapped."""
    ops = []
    for n, op in _CIGAR_RE.findall(cigar):
        mapped = _WFA_OP[op]
        n = int(n)
        if ops and ops[-1][0] == mapped:
            ops[-1] = (mapped, ops[-1][1] + n)
        else:
            ops.append((mapped, n))
    return ops


def consumed(ops):
    """(target_bases, query_bases) consumed by a set of PAF-convention ops."""
    t = sum(n for op, n in ops if op in "=XD")
    q = sum(n for op, n in ops if op in "=XI")
    return t, q


def wfa_binary(explicit=None):
    if explicit:
        return explicit
    here = os.path.dirname(os.path.abspath(__file__))
    local = os.path.join(here, "..", "dev", "wavefront", "wfa_align")
    if os.path.exists(local):
        return os.path.abspath(local)
    found = shutil.which("wfa_align")
    if found:
        return found
    raise FileNotFoundError(
        "wfa_align not found. Build it with:\n"
        "  cd window_identity/dev/wavefront && g++ -O3 -std=c++11 wfa_align.cpp \\\n"
        "    ../../../WFA2-lib/bindings/cpp/WFAligner.cpp -I../../../WFA2-lib \\\n"
        "    -I../../../WFA2-lib/bindings/cpp ../../../WFA2-lib/build/libwfa2.a -o wfa_align"
    )


def _wfa_cmd(opts):
    return [
        wfa_binary(opts.get("wfa_bin")),
        "-x", str(opts.get("wfa_mismatch", 6)),
        "-O", str(opts.get("wfa_gap_open1", 4)),
        "-E", str(opts.get("wfa_gap_ext1", 2)),
        "-o", str(opts.get("wfa_gap_open2", 100)),
        "-e", str(opts.get("wfa_gap_ext2", 1)),
        "-m", opts.get("wfa_model", "affine2p"),
        "-s", opts.get("wfa_span", "end2end"),
        "-M", opts.get("wfa_memory", "ultralow"),
        "-H", opts.get("wfa_heuristic", "none"),
    ]


# WFA2-lib status codes (wavefront/wfa.h)
WF_OK = 0
WF_PARTIAL = 1
_WF_STATUS = {0: "completed", 1: "partial", -100: "max_steps_reached",
              -200: "out_of_memory", -300: "unattainable"}


def _wfa_one(tseq, qseq, opts, _allow_fallback=True):
    """Align one pair of already-oriented sequences; returns ops or None.

    A non-zero WFA status is NOT silently treated as "no alignment". Doing so made a failing
    aligner look like a bad parameter choice: BiWFA (memory_mode=ultralow) returns -300 on
    real diverged segments once gap_open2 reaches ~400 with 2-piece affine, while high/med/low
    all succeed and agree on the score. Swallowed, that reported 0% indel recall for O2>=400
    and would have condemned the penalty rather than the memory mode.

    In BiWFA, -300 is a catch-all (wavefront_bialign.c:1330) for any status that is not
    OK/MAX_STEPS/OOM -- its "unattainable under configured heuristics" message is misleading
    and appears even with heuristics disabled.

    So: fall back to a working memory mode once, and if it still fails, say so loudly.
    """
    if not tseq or not qseq:
        return None
    res = subprocess.run(_wfa_cmd(opts), input=f"job\t{tseq}\t{qseq}\n",
                         capture_output=True, text=True, timeout=opts.get("timer"))
    if res.returncode != 0 or not res.stdout.strip():
        print(f"[wfa] binary failed (rc={res.returncode}): "
              f"{res.stderr.strip()[:200]}", file=sys.stderr)
        return None
    fields = res.stdout.strip().split("\t")
    if len(fields) < 4:
        print(f"[wfa] malformed output: {res.stdout[:120]!r}", file=sys.stderr)
        return None
    try:
        status = int(fields[1])
    except ValueError:
        return None

    if status not in (WF_OK, WF_PARTIAL):
        name = _WF_STATUS.get(status, str(status))
        mem = opts.get("wfa_memory", "ultralow")
        if _allow_fallback and mem == "ultralow" and opts.get("wfa_fallback_memory", "high"):
            retry = dict(opts)
            retry["wfa_memory"] = opts.get("wfa_fallback_memory", "high")
            ops = _wfa_one(tseq, qseq, retry, _allow_fallback=False)
            if ops is not None:
                print(f"[wfa] BiWFA returned {name} on a {len(tseq)}x{len(qseq)}bp pair; "
                      f"succeeded with memory={retry['wfa_memory']}", file=sys.stderr)
                return ops
        print(f"[wfa] alignment FAILED (status={name}, memory={mem}) on a "
              f"{len(tseq)}x{len(qseq)}bp pair -- emitting no record for it",
              file=sys.stderr)
        return None
    if not fields[3]:
        return None
    return parse_wfa_cigar(fields[3]) or None


def align_wfa_global(tseq, qseq, seg, opts) -> List[AlnRecord]:
    """One global end-to-end alignment across the whole segment.

    Only correct when the segment is colinear throughout. Across an internal inversion this
    produces a confident-looking bridge of non-homologous sequence -- see align_wfa().
    """
    ops = _wfa_one(tseq, qseq, opts)
    if not ops:
        return []
    t_used, q_used = consumed(ops)
    return [AlnRecord(q_start=0, q_end=q_used, t_start=0, t_end=t_used,
                      strand="+", ops=ops)]


def align_wfa(tseq, qseq, seg, opts) -> List[AlnRecord]:
    """Chain-aware wavefront alignment (the default).

    Cuts the segment at structural boundaries first, then aligns each maximal colinear chain
    end-to-end. Inversions come out as their own records with strand '-' rather than being
    bridged with garbage. Set opts['wfa_chain']=False to force a single global path.
    """
    if not opts.get("wfa_chain", True):
        return align_wfa_global(tseq, qseq, seg, opts)

    # seeding=True: this is inter-anchoring too, so it must not be bounded by the
    # alignment timer. Without it, a large segment's top-level anchor pass could be
    # killed at --timer and drop through to align_colinear with no chains -- which
    # is the blind-halving path. Observed once on a 35.9 Mb x 37.0 Mb block-mode
    # segment, where a 600 s cap was nowhere near enough.
    anchors = _minimap2_anchors(tseq, qseq, opts, seeding=True)
    chains = chain.colinear_chains(anchors, min_chain_len=opts.get("min_chain_len", 500),
                                   t_len=len(tseq), q_len=len(qseq))
    if not chains:
        # No usable anchors at this level. Subdivide rather than betting the whole
        # segment on one global path that will time out if it is large.
        return [AlnRecord(q_start=sq, q_end=sq + consumed(ops)[1],
                          t_start=st, t_end=st + consumed(ops)[0],
                          strand="+", ops=ops)
                for st, sq, ops in align_colinear(tseq, qseq, opts)]

    out = []
    for c in chains:
        t_sub = tseq[c.t_start:c.t_end]
        q_sub = qseq[c.q_start:c.q_end]
        if c.strand == "-":
            q_sub = _revcomp(q_sub)
        # A chain that will not align in one piece gets subdivided, not dropped.
        # `continue` here was a silent coverage loss: the hardest chains -- the ones
        # carrying the large indels -- are exactly the ones that failed.
        for st, sq, ops in align_colinear(t_sub, q_sub, opts):
            t_used, q_used = consumed(ops)
            if c.strand == "-":
                # ops are in the revcomped chain frame; mirror back into the forward
                # frame, where paf_emit.map_to_parent expects them.
                q_lo = c.q_end - (sq + q_used)
                q_hi = c.q_end - sq
            else:
                q_lo = c.q_start + sq
                q_hi = c.q_start + sq + q_used
            out.append(AlnRecord(q_start=q_lo, q_end=q_hi,
                                 t_start=c.t_start + st,
                                 t_end=c.t_start + st + t_used,
                                 strand=c.strand, ops=ops))
    return out


def _revcomp(s):
    return s.translate(_COMP)[::-1]


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# Subdivision defaults. A pair is attempted directly when t_len * q_len is at or
# under MAX_ALIGN_AREA; larger pairs are seeded and split first. Area rather than
# length is the right gate because WFA cost scales with the product, which is also
# how AnchorWave gates its own novel-anchor pass (`-fa3`, on len_r*len_q > T**2).
DEFAULT_MAX_ALIGN_AREA = 100_000 ** 2
DEFAULT_MAX_SUBDIVIDE_DEPTH = 8
MIN_SUBDIVIDE_LEN = 200  # below this, just align it; seeding cannot help


def _tile_from_chains(chains, t_len, q_len):
    """Ordered sub-intervals that TILE the pair: each chain, plus the gaps around them.

    Returns [(t0, t1, q0, q1), ...] covering [0,t_len) x [0,q_len) with no holes, so
    subdividing can never drop sequence -- the gaps between seeds are exactly the hard
    regions, and they must be recursed into, not discarded.
    """
    tiles = []
    t_at = q_at = 0
    for c in sorted(chains, key=lambda c: (c.t_start, c.q_start)):
        if c.t_start < t_at or c.q_start < q_at:
            continue  # overlaps what we already covered; skip rather than double-count
        if c.t_start > t_at or c.q_start > q_at:
            tiles.append((t_at, c.t_start, q_at, c.q_start))
        tiles.append((c.t_start, c.t_end, c.q_start, c.q_end))
        t_at, q_at = c.t_end, c.q_end
    if t_at < t_len or q_at < q_len:
        tiles.append((t_at, t_len, q_at, q_len))
    return [t for t in tiles if t[1] > t[0] or t[3] > t[2]]


def align_colinear(tseq, qseq, opts, depth=0):
    """Align an oriented, colinear pair, subdividing rather than giving up.

    Returns [(t_start, q_start, ops), ...] in the local frame of (tseq, qseq).

    A timeout or a WFA failure is a signal to re-seed and split, never a reason to
    emit nothing. That is the one behaviour that separated this pipeline from
    AnchorWave: AnchorWave densifies anchors inside a block until every remaining
    interval is small enough to align, and only the interval sizes ever shrink, so
    it always terminates with a full alignment. Measured on this data, one minimap2
    seed pass over a whole chromosome leaves gaps with a median of ~2.5 kb.

    Recursion stops on depth, on a pair too small to seed, or when seeding fails to
    make the interval strictly smaller (no progress).
    """
    if not tseq or not qseq:
        return []

    max_area = opts.get("max_align_area", DEFAULT_MAX_ALIGN_AREA)
    at_floor = (depth >= opts.get("max_subdivide_depth", DEFAULT_MAX_SUBDIVIDE_DEPTH)
                or len(tseq) < MIN_SUBDIVIDE_LEN or len(qseq) < MIN_SUBDIVIDE_LEN)

    # Attempt directly when the pair is small enough to be worth it -- or when it
    # cannot be split any further, in which case the area budget is irrelevant:
    # refusing to try is just dropping the sequence. That was a real defect; a
    # piece below the split floor got no alignment attempt at all.
    if len(tseq) * len(qseq) <= max_area or at_floor:
        try:
            ops = _wfa_one(tseq, qseq, opts)
        except subprocess.TimeoutExpired:
            ops = None
        if ops:
            return [(0, 0, ops)]

    if at_floor:
        return []

    anchors = _minimap2_anchors(tseq, qseq, opts, seeding=True)
    chains = [c for c in chain.colinear_chains(
        anchors, min_chain_len=opts.get("min_chain_len", 500),
        t_len=len(tseq), q_len=len(qseq)) if c.strand == "+"]
    tiles = _tile_from_chains(chains, len(tseq), len(qseq)) if chains else []

    # Seeding can fail to make progress two ways: no usable chain at all, or a
    # single chain spanning the whole pair (the common case -- the pair really is
    # colinear, just too big to align in one go). Returning [] here would be the
    # give-up this function exists to remove, so fall back to halving.
    #
    # AnchorWave's equivalent is its SLIDING_WINDOW branch, which cuts at the
    # argmax of a score-only DP extension. A proportional midpoint is cruder, but
    # it guarantees both halves are strictly smaller, so recursion terminates, and
    # each half is re-seeded on the way down -- which usually finds the anchors the
    # whole-span pass could not resolve.
    if not tiles or (len(tiles) == 1 and tiles[0] == (0, len(tseq), 0, len(qseq))):
        # Cut at the seed alignment's own large indels first. A blind midpoint
        # mis-pairs the halves whenever a big indel sits between them: on a
        # 20 kb x 23 kb pair carrying one 3 kb insertion it tiled the reference
        # fully but matched only 4,909 of 20,000 bases, because every piece below
        # the first cut was offset by the insertion and WFA end2end dutifully
        # aligned the wrong sequences to each other. The CIGAR knows exactly where
        # the query jumps, so cutting there keeps both sides in register.
        pts = _seed_breakpoints(tseq, qseq, opts)
        if pts:
            tiles = []
            t_at = q_at = 0
            for t_cut, q_cut in pts:
                if t_cut > t_at or q_cut > q_at:
                    tiles.append((t_at, t_cut, q_at, q_cut))
                    t_at, q_at = t_cut, q_cut
            if t_at < len(tseq) or q_at < len(qseq):
                tiles.append((t_at, len(tseq), q_at, len(qseq)))
        if not tiles or (len(tiles) == 1
                         and tiles[0] == (0, len(tseq), 0, len(qseq))):
            t_mid = len(tseq) // 2
            q_mid = max(1, min(len(qseq) - 1,
                               len(qseq) * t_mid // max(1, len(tseq))))
            tiles = [(0, t_mid, 0, q_mid), (t_mid, len(tseq), q_mid, len(qseq))]

    out = []
    for t0, t1, q0, q1 in tiles:
        for st, sq, ops in align_colinear(tseq[t0:t1], qseq[q0:q1], opts, depth + 1):
            out.append((t0 + st, q0 + sq, ops))
    return out


_CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def _seed_breakpoints(tseq, qseq, opts, min_indel=100):
    """(t_pos, q_pos) cut points at the large indels of the best seed alignment.

    Positions are local to (tseq, qseq) and paired, so cutting there leaves both
    sides in register -- which a proportional midpoint does not once an indel sits
    between the halves. Only forward-strand records are used; align_colinear is
    given an already-oriented pair.
    """
    paf = _minimap2_paf(tseq, qseq, opts, seeding=True)
    best = None
    for line in paf.strip().split("\n"):
        c = line.split("\t")
        if len(c) < 12 or c[4] != "+":
            continue
        try:
            score = int(c[9])
        except ValueError:
            continue
        if best is None or score > best[0]:
            best = (score, c)
    if best is None:
        return []
    c = best[1]
    cg = next((f[5:] for f in c[12:] if f.startswith("cg:Z:")), None)
    if not cg:
        return []
    t, q = int(c[7]), int(c[2])
    pts = []
    for num, op in _CIGAR_OP.findall(cg):
        n = int(num)
        if op in ("M", "=", "X"):
            t += n
            q += n
        elif op in ("D", "N"):          # reference advances, query does not
            if n >= min_indel:
                pts.append((t, q))
            t += n
            if n >= min_indel:
                pts.append((t, q))
        elif op == "I":                 # query advances, reference does not
            if n >= min_indel:
                pts.append((t, q))
            q += n
            if n >= min_indel:
                pts.append((t, q))
    return pts


def _minimap2_paf(tseq, qseq, opts, seeding=False):
    """Raw seed PAF for this pair; '' on any failure."""
    d = tempfile.mkdtemp(dir=opts.get("temp_base") or None, prefix="seed_")
    try:
        t_fa, q_fa = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        with open(t_fa, "w") as fh:
            fh.write(f">t\n{tseq}\n")
        with open(q_fa, "w") as fh:
            fh.write(f">q\n{qseq}\n")
        # The landmark pass is UNBOUNDED by default. It is the thing that makes
        # subdivision safe: without landmarks the only fallback is a proportional
        # midpoint, which mis-pairs the halves across a large indel and lets an
        # end2end aligner confidently align the wrong sequences together. A 60 s
        # budget silently pushed ~200 of the largest pairs per run down that path.
        # Letting minimap2 take the minutes it needs is always cheaper than the
        # alignment it protects. --seed-timer can still impose one deliberately.
        budget = opts.get("seed_timer") if seeding else opts.get("timer")
        try:
            res = subprocess.run(
                [opts.get("minimap2_bin", "minimap2"), "-t", str(opts.get("threads", 1)),
                 "--secondary=no", "-x", opts.get("preset", "asm20"), "-c", t_fa, q_fa],
                capture_output=True, text=True, timeout=budget)
        except (subprocess.TimeoutExpired, OSError) as exc:
            print(f"[seed] minimap2 pass failed on a {len(tseq)}x{len(qseq)}bp "
                  f"pair: {type(exc).__name__} -- falling back to proportional "
                  f"halving, which is NOT indel-aware", file=sys.stderr)
            return ""
        return res.stdout if res.returncode == 0 else ""
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _minimap2_anchors(tseq, qseq, opts, seeding=False):
    """Anchor pass. tseq/qseq are as submitted, so anchors are window-local.

    `-c` is required, not an optimisation. Without base-level alignment minimap2 reports the
    approximate CHAIN SPAN, which brackets an internal inversion inside a single '+' record
    and hides the very boundary we are looking for. With `-c` the same input resolves into
    '+' / '-' / '+', which is what makes chaining work at all.

    `--secondary=no` is also load-bearing: secondary records duplicate the inverted span, and
    a duplicate chain would align the same bases twice and double-count them in the M-sum
    that weights every downstream divergence average.
    """
    d = tempfile.mkdtemp(dir=opts.get("temp_base") or None, prefix="chain_")
    try:
        t_fa, q_fa = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        with open(t_fa, "w") as fh:
            fh.write(f">t\n{tseq}\n")
        with open(q_fa, "w") as fh:
            fh.write(f">q\n{qseq}\n")
        # The seed pass gets its own, shorter budget. Sharing opts['timer'] with the
        # alignment meant the seeder could be killed by the very timeout it exists to
        # prevent -- and a dead seeder turns into a skipped segment.
        budget = opts.get("seed_timer") if seeding else opts.get("timer")
        try:
            res = subprocess.run(
                [opts.get("minimap2_bin", "minimap2"), "-t", str(opts.get("threads", 1)),
                 "--secondary=no", "-x", opts.get("preset", "asm20"), "-c", t_fa, q_fa],
                capture_output=True, text=True, timeout=budget)
        except (subprocess.TimeoutExpired, OSError) as exc:
            # Returning [] lets the caller fall back; raising would abort the segment.
            print(f"[seed] minimap2 anchor pass failed on a {len(tseq)}x{len(qseq)}bp "
                  f"pair: {type(exc).__name__}", file=sys.stderr)
            return []
        if res.returncode != 0:
            return []
        return chain.parse_paf_anchors(res.stdout, min_len=opts.get("min_anchor_len", 200))
    finally:
        shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------------------------------- LAST


def _parse_maf(text) -> List[AlnRecord]:
    """Parse lastal MAF into PAF-convention ops.

    MAF gives gapped alignment rows; the first 's' line is the target (the lastdb subject),
    the second is the query. lastal searches both strands natively, so orientation comes out
    of the alignment and no separate orientation pass is needed for this backend.
    """
    records = []
    block = []
    for line in text.split("\n"):
        if line.startswith("a "):
            block = []
        elif line.startswith("s "):
            block.append(line.split())
            if len(block) == 2:
                records.append(_maf_block(block))
                block = []
    return [r for r in records if r]


def _maf_block(block):
    # s name start alnSize strand srcSize sequence
    _, _, t_start, t_size, t_strand, t_srclen, t_aln = block[0]
    _, _, q_start, q_size, q_strand, q_srclen, q_aln = block[1]
    t_start, q_start = int(t_start), int(q_start)
    t_size, q_size = int(t_size), int(q_size)
    q_srclen = int(q_srclen)

    ops = []
    for tc, qc in zip(t_aln, q_aln):
        if tc == "-":
            op = "I"
        elif qc == "-":
            op = "D"
        else:
            op = "=" if tc.upper() == qc.upper() else "X"
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + 1)
        else:
            ops.append((op, 1))
    if not ops:
        return None

    # MAF minus-strand coordinates are given on the reverse strand; convert to forward.
    if q_strand == "-":
        q_fwd_start = q_srclen - (q_start + q_size)
    else:
        q_fwd_start = q_start
    return AlnRecord(q_start=q_fwd_start, q_end=q_fwd_start + q_size,
                     t_start=t_start, t_end=t_start + t_size,
                     strand=("-" if q_strand != t_strand else "+"), ops=ops)


def align_last(tseq, qseq, seg, opts) -> List[AlnRecord]:
    d = tempfile.mkdtemp(dir=opts["temp_base"], prefix="last_")
    try:
        t_fa, q_fa = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        with open(t_fa, "w") as fh:
            fh.write(f">t\n{tseq}\n")
        with open(q_fa, "w") as fh:
            fh.write(f">q\n{qseq}\n")
        db = os.path.join(d, "db")
        subprocess.run(["lastdb", "-uNEAR" if opts.get("last_near") else "-uYASS", db, t_fa],
                       capture_output=True, check=True, timeout=opts.get("timer"))
        cmd = ["lastal"]
        for flag, key in (("-p", "last_matrix"), ("-a", "last_gap_open"),
                          ("-b", "last_gap_ext"), ("-P", "last_threads")):
            if opts.get(key) is not None:
                cmd += [flag, str(opts[key])]
        cmd += [db, q_fa]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=opts.get("timer"))
        if res.returncode != 0:
            return []
        maf = res.stdout
        if opts.get("last_split"):
            sp = subprocess.run(["last-split"], input=maf, capture_output=True, text=True,
                                timeout=opts.get("timer"))
            if sp.returncode == 0:
                maf = sp.stdout
        return _parse_maf(maf)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        raise
    finally:
        if not opts.get("debug"):
            shutil.rmtree(d, ignore_errors=True)
