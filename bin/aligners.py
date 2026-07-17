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

    anchors = _minimap2_anchors(tseq, qseq, opts)
    chains = chain.colinear_chains(anchors, min_chain_len=opts.get("min_chain_len", 500),
                                   t_len=len(tseq), q_len=len(qseq))
    if not chains:
        # No usable anchors: fall back to a single global path rather than emitting nothing.
        return align_wfa_global(tseq, qseq, seg, opts)

    out = []
    for c in chains:
        t_sub = tseq[c.t_start:c.t_end]
        q_sub = qseq[c.q_start:c.q_end]
        if c.strand == "-":
            q_sub = _revcomp(q_sub)
        ops = _wfa_one(t_sub, q_sub, opts)
        if not ops:
            continue
        t_used, q_used = consumed(ops)
        if c.strand == "-":
            # ops were computed against the revcomped chain; the query interval is still
            # [c.q_start, c.q_end) in the forward frame, and paf_emit maps it from there.
            q_lo = c.q_end - q_used
            q_hi = c.q_end
        else:
            q_lo = c.q_start
            q_hi = c.q_start + q_used
        out.append(AlnRecord(q_start=q_lo, q_end=q_hi,
                             t_start=c.t_start, t_end=c.t_start + t_used,
                             strand=c.strand, ops=ops))
    return out


def _revcomp(s):
    return s.translate(_COMP)[::-1]


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _minimap2_anchors(tseq, qseq, opts):
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
        res = subprocess.run(
            [opts.get("minimap2_bin", "minimap2"), "-t", str(opts.get("threads", 1)),
             "--secondary=no", "-x", opts.get("preset", "asm20"), "-c", t_fa, q_fa],
            capture_output=True, text=True, timeout=opts.get("timer"))
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
