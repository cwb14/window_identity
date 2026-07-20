#!/usr/bin/env python3
"""Extract syntenic segment pairs and align them, emitting contract-compliant PAF.

Backends are selected with --aligner. Each one only produces alignments; all coordinate
mapping and tag construction live in paf_emit.py, so the arithmetic is written and tested
once instead of once per backend. See step10_dev/memo_step10_diagnosis.md for the audit
behind that split and for the PAF contract Steps 11-17 depend on.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from multiprocessing import Pool, Manager

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import aligners
import fastaio
import paf_emit
from orient import orientation, revcomp

VERBOSE = False
_FASTA_CACHE = {}


def log(msg):
    if VERBOSE:
        print(msg)


def warn(msg):
    print(msg, file=sys.stderr)


# ------------------------------------------------------------------ sequence access


def get_fasta(path):
    """Indexed, per-process handle. The old code re-parsed the whole genome with SeqIO for
    every coords line -- 3088 full passes over an 18-78MB FASTA."""
    fa = _FASTA_CACHE.get(path)
    if fa is None:
        fa = pysam.FastaFile(path)
        _FASTA_CACHE[path] = fa
    return fa


def parse_syncoord(syncoord, work_dir="."):
    """'AthaAt4_chr4:11535683..11542200' -> ('AthaAt4', 'AthaAt4_chr4', 11535683, 11542200).

    Coordinates are 1-based inclusive, as written by gene_coords_extractor_all4.py.
    The accession is resolved against the canonical genome ID list, so IDs containing
    underscores ('annuaA_chr_chr1') resolve correctly.
    """
    name, rng = syncoord.split(":")
    start, end = (int(x) for x in rng.split(".."))
    accession = fastaio.accession_of(name, fastaio.genome_ids(work_dir))
    return accession, name, start, end


def extract(genome_dir, syncoord, do_revcomp=False, upper=True):
    accession, name, start, end = parse_syncoord(syncoord, genome_dir)
    fa = get_fasta(os.path.join(genome_dir, f"{accession}_mod.fa"))
    seq = fa.fetch(name, start - 1, end)  # 1-based inclusive -> 0-based half-open
    if upper:
        # WFA compares raw bytes and these genomes are ~20% soft-masked, so 'a' vs 'A'
        # would score as a mismatch. minimap2 is case-insensitive, so this is a no-op there.
        seq = seq.upper()
    return revcomp(seq) if do_revcomp else seq


# ------------------------------------------------------------------ backends


def align_minimap2(tseq, qseq, seg, opts):
    """Native minimap2, coordinates corrected, every other field untouched.

    Deliberately does not go through the emitter: keeping minimap2's own tags byte-identical
    is what lets the regression gate prove the coordinate fix did not move de/k2p.
    """
    d = tempfile.mkdtemp(dir=opts["temp_base"])
    try:
        t_fa, q_fa = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        with open(t_fa, "w") as fh:
            fh.write(f">{seg.t_name}\n{tseq}\n")
        with open(q_fa, "w") as fh:
            fh.write(f">{seg.q_name}\n{qseq}\n")
        cmd = [opts["minimap2_bin"], "-t", str(opts["threads"]), "--secondary=no",
               "--cs=short", "-x", opts["preset"], "-c", t_fa, q_fa]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=opts["timer"])
        if res.returncode != 0:
            log(f"minimap2 error: {res.stderr}")
            return []
        out = []
        for line in res.stdout.strip().split("\n"):
            if line:
                out.append(paf_emit.shift_native_paf(line, seg))
        return out
    finally:
        if not opts["debug"]:
            for f in ("t.fa", "q.fa"):
                try:
                    os.remove(os.path.join(d, f))
                except OSError:
                    pass
            try:
                os.rmdir(d)
            except OSError:
                pass


def _emit(records, seg, tseq, qseq):
    return [paf_emit.format_paf(r, seg, tseq, qseq) for r in records]


def align_wfa_lines(tseq, qseq, seg, opts):
    return _emit(aligners.align_wfa(tseq, qseq, seg, opts), seg, tseq, qseq)


def align_last_lines(tseq, qseq, seg, opts):
    return _emit(aligners.align_last(tseq, qseq, seg, opts), seg, tseq, qseq)


BACKENDS = {
    "minimap2": align_minimap2,
    "last": align_last_lines,
    "wfa": align_wfa_lines,
}


# ------------------------------------------------------------------ per-segment driver


def build_segment(syncoord1, syncoord2, q_revcomp, strand):
    _, t_name, t_start, t_end = parse_syncoord(syncoord1)
    _, q_name, q_start, q_end = parse_syncoord(syncoord2)
    return paf_emit.Segment(
        t_name, t_start, t_end, q_name, q_start, q_end, q_revcomp,
        paf_emit.segment_id(syncoord1, syncoord2, strand),
    )


def length_ratio(seg):
    lo, hi = sorted((seg.t_len, seg.q_len))
    return float("inf") if lo <= 0 else hi / lo


def guard_legacy_untagged(paf_path, p):
    """Abort --resume if paf_path holds well-formed rows written before sd:Z: existed.

    A row with >=12 tab-separated columns but no trailing sd:Z: tag is ambiguous by
    itself: it can be legacy data (written by a version of this pipeline that predates
    the tag), or it can be the one row a SIGKILL tore mid-write. Those two cases must be
    told apart by WHERE the untagged rows sit, not merely whether any exist:

    * cg:Z: and cs:Z: are roughly 99% of a row's bytes, so a SIGKILL landing inside a
      worker's flush almost always tears inside one of those fields -- not at the exact
      tab before sd:Z:. The result is a partial row that still has 20-24 columns (past
      the >=12 threshold below) but is missing the tag, because the tag is appended
      after cs:Z: and the write never got that far. Because every writer opens the
      output file, appends, and the process tree dies together, at most one row can ever
      be torn this way, and it is necessarily the LAST line of the file -- nothing can be
      appended after the kill.
    * A legacy pre-tag file has no such constraint: every row it holds was written
      without the tag, scattered across the whole file, not just the last line (see
      /data2/chris/poa3/PannuaA/alignment.paf for a real example -- 13446 rows, none
      tagged, no .params sidecar).

    So: an untagged row at any position OTHER than the final line can only be legacy
    data, since a tear can never land there. A file whose rows are ALL untagged is also
    legacy even when it has exactly one row -- "only the final line is untagged" would
    otherwise be indistinguishable from a single-row legacy file, and silently pruning
    that row to an empty PAF is just as wrong as pruning a multi-row legacy file.

    done_segment_ids() can only recognize done segments by their sd:Z: tag, so it cannot
    tell which segments legacy rows already cover. Blindly keeping them would mean
    --resume re-aligns every segment and appends new tagged rows alongside the old
    untagged ones, silently mixing two runs' output in one file, which is the exact
    corruption this whole feature exists to prevent. So: refuse and let the user decide,
    the same idiom as the parameter-mismatch and corrupt-sidecar guards below.

    A short, genuinely incomplete line (<12 columns -- cut off mid-write before even the
    fixed leading columns finished) is NOT flagged here; that is prune_torn_tail's job
    below.
    """
    if not os.path.exists(paf_path):
        return
    lines = []
    with open(paf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                lines.append(line)
    if not lines:
        return
    total = len(lines)
    last_idx = total - 1
    # Rows with no sd:Z: tag at all -- includes a short torn tail (<12 columns), which
    # is real "no tag" data even though it is not legacy and must not affect the trigger
    # logic below. Kept separate from untagged_idx so the abort message's count reflects
    # every untagged row, not just the ones that can carry the trigger.
    no_tag_idx = [
        i for i, line in enumerate(lines)
        if not line.split("\t")[-1].startswith(paf_emit.SEG_TAG)
    ]
    untagged_idx = [i for i in no_tag_idx if len(lines[i].split("\t")) >= 12]
    if not untagged_idx:
        return
    non_final_untagged = any(i != last_idx for i in untagged_idx)
    all_untagged = len(untagged_idx) == total
    if non_final_untagged or all_untagged:
        no_tag = len(no_tag_idx)
        p.error(
            f"{paf_path} has {total} row(s), {no_tag} of which have no sd:Z: tag.\n"
            f"That tag is how --resume tells which segments are already done; rows "
            f"without it were written by a version of this pipeline that predates the "
            f"tag. --resume cannot tell which segments those rows cover, so it will not "
            f"guess and risk silently mixing old and new rows into one {paf_path}.\n"
            f"Either treat {paf_path} as already complete and do not pass --resume, or "
            f"move it aside and re-run without --resume to regenerate it with tags from "
            f"scratch."
        )


def prune_torn_tail(paf_path):
    """Rewrite paf_path without lines that lack a parseable sd:Z: tag. Returns the count dropped.

    A SIGKILL landing inside a worker's flush can leave a torn final line. Rows here reach
    257 KB, far past the 8 KiB buffer, so a partial raw write is possible in principle. The
    torn line must be removed rather than merely ignored: leaving it in place would feed a
    malformed record to Step 11. Cost is one segment, which is simply re-aligned.

    Callers must run guard_legacy_untagged() on paf_path first. By the time this runs, any
    remaining untagged line is guaranteed short (<12 columns) rather than legacy data, so
    the two conditions below never disagree on a real input -- kept here as a second check
    rather than relying solely on length so a future caller that skips the guard fails
    closed (drops the line) instead of open (keeps an untagged row un-flagged).
    """
    if not os.path.exists(paf_path):
        return 0
    kept, dropped = [], 0
    with open(paf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) >= 12 and cols[-1].startswith(paf_emit.SEG_TAG):
                kept.append(line)
            else:
                dropped += 1
    if dropped:
        tmp = paf_path + ".partial"
        with open(tmp, "w") as fh:
            for line in kept:
                fh.write(line + "\n")
        os.replace(tmp, paf_path)  # atomic
    return dropped


def done_segment_ids(paf_path, skipped_path):
    """Segment IDs already accounted for, from both the PAF and the skip sidecar.

    Every segment ends in exactly one of the two files, so their union is the complete
    done-set and nothing is retried forever.

    Note: a Segment built without a seg_id (paf_emit.Segment defaults seg_id="") produces an
    untagged format_paf row, which the `cols[-1].startswith(SEG_TAG)` check below silently
    excludes from the done-set. That segment would simply be re-aligned as a duplicate rather
    than corrupted. Not reachable today via build_segment, which always sets seg_id.
    """
    done = set()
    if os.path.exists(paf_path):
        with open(paf_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[-1].startswith(paf_emit.SEG_TAG):
                    done.add(cols[-1][len(paf_emit.SEG_TAG):])
    if os.path.exists(skipped_path):
        with open(skipped_path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3 or cols[0] == "coord1":
                    continue
                done.add(paf_emit.segment_id(cols[0], cols[1], cols[2]))
    return done


def record_skip(skipped_path, syncoord1, syncoord2, strand, reason, ratio, len1, len2):
    """Append one skip row immediately.

    The old code buffered these in a Manager().list() and serialized only after pool.map
    returned, so a killed run lost every skip decision and resume could not tell a skipped
    segment from an unattempted one. Concurrent O_APPEND writes from pool workers are atomic
    on Linux -- not because of PIPE_BUF (that guarantee is for pipes/FIFOs), but because the
    kernel holds the inode lock for the duration of each write(2) to a regular file.
    """
    with open(skipped_path, "a") as fh:
        fh.write(f"{syncoord1}\t{syncoord2}\t{strand}\t{reason}\t{ratio}\t{len1}\t{len2}\n")


# ------------------------------------------------------------------ resume parameter guard


def params_sidecar_path(output_path):
    """Sidecar recording the alignment-affecting parameters used to produce output_path."""
    return output_path + ".params"


def align_params(args, timer):
    """The subset of main()'s `opts` that can change alignment output or skip decisions.

    Covers the aligner choice, preset, the length-ratio and timeout skip thresholds, the
    minimap2 binary, and every backend tuning knob actually read out of `opts` in main()
    (wfa_*, last_*, min_chain_len, min_anchor_len). minimap2_bin is included even though
    minimap2 is only one of three --aligner choices: align_minimap2() invokes it directly,
    and orientation() (used by the last/wfa backends to settle q_revcomp) also shells out
    to it, so a different binary can change alignment output or trigger a no_orientation
    skip regardless of which --aligner is selected. Deliberately excludes threads,
    processes, debug, verbose, and temp_base -- none of those can change what gets
    written, only how fast or how loudly the run reports itself.
    """
    return {
        "aligner": args.aligner, "preset": args.preset,
        "max_len_ratio": args.max_len_ratio, "timer": timer,
        "minimap2_bin": args.minimap2_bin,
        "wfa_bin": args.wfa_bin, "wfa_model": args.wfa_model, "wfa_span": args.wfa_span,
        "wfa_memory": args.wfa_memory, "wfa_heuristic": args.wfa_heuristic,
        "wfa_chain": args.wfa_chain, "wfa_mismatch": args.wfa_mismatch,
        "wfa_gap_open1": args.wfa_gap_open1, "wfa_gap_ext1": args.wfa_gap_ext1,
        "wfa_gap_open2": args.wfa_gap_open2, "wfa_gap_ext2": args.wfa_gap_ext2,
        "min_chain_len": args.min_chain_len, "min_anchor_len": args.min_anchor_len,
        "last_matrix": args.last_matrix, "last_gap_open": args.last_gap_open,
        "last_gap_ext": args.last_gap_ext, "last_near": args.last_near,
        "last_split": args.last_split,
    }


def write_params_sidecar(output_path, params):
    """Record this run's alignment parameters so a later --resume can be checked against them.

    Atomic, same pattern as prune_torn_tail(): write to a temp path in the same
    directory, then os.replace(). Without this, a kill mid-write leaves a truncated,
    unparseable sidecar -- which is exactly the corrupt-file case load_params_sidecar()
    below has to guard against, so the write side should not be the thing causing it.
    """
    path = params_sidecar_path(output_path)
    tmp = path + ".partial"
    with open(tmp, "w") as fh:
        json.dump(params, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)  # atomic


def load_params_sidecar(output_path, p):
    """Previously recorded parameters, or None if the sidecar does not exist (older run dirs).

    `p` is the ArgumentParser, used to abort via p.error() -- the same idiom the
    parameter-mismatch check below uses -- if the sidecar exists but cannot be parsed
    (hand-edited, truncated by a kill outside of write_params_sidecar's atomic replace,
    disk corruption). Without this, json.load() raises JSONDecodeError uncaught and the
    user sees a raw traceback instead of a clean, actionable message.
    """
    path = params_sidecar_path(output_path)
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        p.error(
            f"cannot parse parameter sidecar {path}: {exc}\n"
            f"A corrupt sidecar cannot be trusted to guarantee --resume is using the "
            f"same alignment parameters as the original run.\n"
            f"Either delete {path} and re-run without --resume to start fresh, or "
            f"restore the original sidecar file."
        )


def diff_params(old, new):
    """{key: (old_value, new_value)} for every key present in both dicts whose value differs.

    A key only `new` has (an older sidecar predating some parameter) is silently skipped --
    there is no recorded baseline for it, so treating that as a mismatch would false-positive
    on every tool upgrade rather than on an actual parameter change.
    """
    return {k: (old[k], new[k]) for k in new if k in old and old[k] != new[k]}


def _read_skip_rows(path):
    """Data rows from the skip sidecar, with the header and any malformed lines excluded."""
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return [r for r in rows if len(r) >= 4 and r[0] != "coord1"]


def process_line(args):
    line, genome_dir, opts, counter, skipped = args
    try:
        syncoord1, syncoord2, strand = line.strip().split("\t")
    except ValueError:
        warn(f"malformed coords line, skipping: {line.strip()!r}")
        counter.value += 1
        return

    try:
        seg_probe = build_segment(syncoord1, syncoord2, False, strand)

        # Guard before doing any work. A 242:1 span is not alignable end-to-end; riparian's
        # rescue machinery deliberately keeps such blocks for plot completeness, which is
        # the opposite of what an aligner wants, so only its ratio test is borrowed.
        if seg_probe.t_len <= 0 or seg_probe.q_len <= 0:
            record_skip(skipped, syncoord1, syncoord2, strand, "empty_span", 0,
                        seg_probe.t_len, seg_probe.q_len)
            counter.value += 1
            return
        ratio = length_ratio(seg_probe)
        if opts["max_len_ratio"] > 0 and ratio > opts["max_len_ratio"]:
            record_skip(skipped, syncoord1, syncoord2, strand, "length_ratio", f"{ratio:.1f}",
                        seg_probe.t_len, seg_probe.q_len)
            counter.value += 1
            return

        tseq = extract(genome_dir, syncoord1)

        if opts["aligner"] == "minimap2":
            # Preserve historical behaviour exactly: trust the strand column, and let
            # minimap2's own both-strand search cover it being wrong.
            q_revcomp = (strand == "-")
        else:
            # Strand-sensitive backend: settle orientation from sequence.
            q_fwd = extract(genome_dir, syncoord2, do_revcomp=False)
            needs_rc = orientation(tseq, q_fwd, opts["minimap2_bin"], opts["preset"])
            if needs_rc is None:
                record_skip(skipped, syncoord1, syncoord2, strand, "no_orientation",
                            f"{ratio:.1f}", seg_probe.t_len, seg_probe.q_len)
                counter.value += 1
                return
            q_revcomp = needs_rc

        seg = build_segment(syncoord1, syncoord2, q_revcomp, strand)
        qseq = extract(genome_dir, syncoord2, do_revcomp=q_revcomp)

        try:
            lines = BACKENDS[opts["aligner"]](tseq, qseq, seg, opts)
        except subprocess.TimeoutExpired:
            record_skip(skipped, syncoord1, syncoord2, strand, "timeout", f"{ratio:.1f}",
                        seg_probe.t_len, seg_probe.q_len)
            counter.value += 1
            return

        if lines:
            with open(opts["output_file"], "a") as fh:
                fh.write("\n".join(lines) + "\n")
        else:
            # The aligner ran and produced nothing (empty stdout or non-zero rc, see
            # align_minimap2). Recording it is what makes "attempted, produced nothing"
            # distinguishable from "never attempted" -- otherwise resume retries it forever.
            record_skip(skipped, syncoord1, syncoord2, strand, "no_alignment",
                        f"{ratio:.1f}", seg_probe.t_len, seg_probe.q_len)
    except Exception as exc:  # noqa: BLE001 - one bad segment must not kill the pool
        warn(f"error on {line.strip()!r}: {exc}")
    finally:
        counter.value += 1


# ------------------------------------------------------------------ main


def main():
    p = argparse.ArgumentParser(
        description="Extract syntenic segments and align them, emitting PAF.")
    group = p.add_mutually_exclusive_group()
    group.add_argument("-pairedIDs", nargs="+", help="Accession IDs to include.")
    group.add_argument("-singleID", help="Single accession ID to include.")
    p.add_argument("-t", "--threads", type=int, default=3, help="Threads per aligner call.")
    p.add_argument("-p", "--processes", type=int, default=10, help="Parallel processes.")
    p.add_argument("-c", "--coords", required=True, help="Syntenic coordinates file.")
    p.add_argument("--timer", type=str, help="Per-alignment timeout, e.g. 10m.")
    p.add_argument("--preset", choices=["asm5", "asm10", "asm20"], default="asm10",
                   help="minimap2 preset (default: asm10).")
    p.add_argument("--aligner", choices=sorted(BACKENDS), default="minimap2",
                   help="Alignment backend (default: minimap2).")
    p.add_argument("--minimap2-bin", default="minimap2",
                   help="minimap2 binary (default: minimap2 on PATH).")
    p.add_argument("--max-len-ratio", type=float, default=5.0,
                   help="Skip segments whose length ratio exceeds this; 0 disables "
                        "(default: 5.0).")
    wfa = p.add_argument_group("wfa backend")
    wfa.add_argument("--wfa-bin", default=None, help="wfa_align binary (default: dev/wavefront/wfa_align).")
    wfa.add_argument("--wfa-model", choices=["affine", "affine2p"], default="affine2p")
    wfa.add_argument("--wfa-span", choices=["end2end", "endsfree"], default="end2end")
    wfa.add_argument("--wfa-memory", choices=["high", "med", "low", "ultralow"],
                     default="ultralow", help="ultralow = BiWFA (default); required for Mb-scale.")
    wfa.add_argument("--wfa-heuristic", choices=["none", "adaptive"], default="none",
                     help="none (default). 'adaptive' is WFA2-lib's own default and costs "
                          "~14 identity points on 30kb segments -- exposed for benchmarking only.")
    wfa.add_argument("--no-wfa-chain", dest="wfa_chain", action="store_false",
                     help="Force ONE global path per segment instead of chaining. An "
                          "internal inversion then gets reported as a large indel pair "
                          "(a false TE-indel call) -- benchmarking only.")
    wfa.add_argument("--min-chain-len", type=int, default=500,
                     help="Discard colinear chains shorter than this (default: 500).")
    wfa.add_argument("--min-anchor-len", type=int, default=200,
                     help="Discard anchor records shorter than this (default: 200).")
    wfa.add_argument("--wfa-mismatch", type=int, default=6)
    wfa.add_argument("--wfa-gap-open1", type=int, default=4)
    wfa.add_argument("--wfa-gap-ext1", type=int, default=2)
    wfa.add_argument("--wfa-gap-open2", type=int, default=100)
    wfa.add_argument("--wfa-gap-ext2", type=int, default=1)

    last = p.add_argument_group("last backend")
    last.add_argument("--last-matrix", default=None, help="lastal -p scoring matrix (e.g. from last-train).")
    last.add_argument("--last-gap-open", type=int, default=None, help="lastal -a.")
    last.add_argument("--last-gap-ext", type=int, default=None, help="lastal -b.")
    last.add_argument("--last-near", action="store_true", help="lastdb -uNEAR (similar sequences).")
    last.add_argument("--last-split", action="store_true", help="Pipe lastal through last-split.")

    p.add_argument("--output", default="alignment.paf", help="Output PAF.")
    p.add_argument("--skipped", default="alignment_skipped.tsv",
                   help="Sidecar listing skipped segments and why.")
    p.add_argument("--resume", action="store_true",
                   help="Reuse an existing alignment.paf: align only segments not already "
                        "recorded in it or in the skip sidecar. Without this, both files are "
                        "truncated at startup.")
    p.add_argument("--genome-dir", default=".", help="Directory holding <ID>_mod.fa.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("-d", "--debug", action="store_true", help="Retain temp files.")
    args = p.parse_args()

    global VERBOSE
    VERBOSE = args.verbose

    units = {"d": 86400, "h": 3600, "m": 60, "s": 1}
    timer = None
    if args.timer:
        if args.timer[-1] not in units:
            p.error(f"invalid time unit in '{args.timer}'; use d/h/m/s")
        timer = int(args.timer[:-1]) * units[args.timer[-1]]

    with open(args.coords) as fh:
        lines = [l for l in fh if l.strip()]
    filtered = []
    for line in lines:
        syn1, syn2, _ = line.strip().split("\t")
        acc1 = fastaio.accession_of(syn1.split(":")[0], fastaio.genome_ids(args.genome_dir))
        acc2 = fastaio.accession_of(syn2.split(":")[0], fastaio.genome_ids(args.genome_dir))
        if args.pairedIDs:
            if acc1 in args.pairedIDs and acc2 in args.pairedIDs:
                filtered.append(line)
        elif args.singleID:
            if args.singleID in (acc1, acc2):
                filtered.append(line)
        else:
            filtered.append(line)

    if args.resume:
        # guard_legacy_untagged() must run before ANYTHING below writes to the run
        # directory. It is a pure check (reads paf_path, otherwise inert) and every
        # branch that follows -- the parameter compare, prune_torn_tail, the sidecar
        # write -- either mutates a file or records state that a later --resume trusts.
        # Running the guard first means a legacy-PAF abort leaves the directory exactly
        # as it found it: no sidecar, no pruning, nothing for a subsequent --resume to
        # misread as anything other than the same legacy PAF.
        guard_legacy_untagged(args.output, p)

        current_params = align_params(args, timer)
        old_params = load_params_sidecar(args.output, p)
        sidecar = params_sidecar_path(args.output)
        if old_params is None:
            # The shell driver defaults to -resume yes, so --resume is the common case,
            # not the exception -- warning here every time would train the user to
            # ignore the one warning that matters. A missing sidecar is only actually
            # informative when args.output already exists: that means some prior run
            # produced it before this guard existed, and its parameters genuinely
            # cannot be verified. A --resume into a directory with no prior output has
            # nothing to have drifted from, so it stays silent. Either way nothing is
            # written here -- see the single write at the end of this branch.
            if os.path.exists(args.output):
                warn(f"WARNING: {sidecar} not found; cannot verify --resume is using the "
                     f"same alignment parameters as the original run (older run "
                     f"directories predate this check). Proceeding.")
        else:
            mismatches = diff_params(old_params, current_params)
            if mismatches:
                detail = "\n".join(
                    f"  {k}: was {old!r}, now {new!r}"
                    for k, (old, new) in sorted(mismatches.items())
                )
                p.error(
                    f"--resume parameters differ from the run recorded in {sidecar}:\n"
                    f"{detail}\n"
                    f"Mixing results computed under different settings into one "
                    f"{args.output} would be silently wrong.\n"
                    f"Either re-run without --resume to start fresh, or restore the "
                    f"original parameters."
                )
            # Parameters agree on every key the old sidecar has; the rewrite that used to
            # happen here now happens once, at the end of this branch (see below), along
            # with the old_params-is-None case above -- so a key added by a newer version
            # of this script (silently skipped by diff_params -- see its docstring) is
            # still captured for the next --resume.

        dropped = prune_torn_tail(args.output)
        if dropped:
            print(f"Dropped {dropped} unparseable line(s) from {args.output}")
        done = done_segment_ids(args.output, args.skipped)
        before = len(filtered)
        filtered = [
            l for l in filtered
            if paf_emit.segment_id(*l.strip().split("\t")[:3]) not in done
        ]
        print(f"Resume: {before - len(filtered)} of {before} segments already done, "
              f"{len(filtered)} remaining")
        if not os.path.exists(args.skipped):
            with open(args.skipped, "w") as fh:
                fh.write("coord1\tcoord2\tstrand\treason\tratio\tlen1\tlen2\n")

        # Single write point, reached only once the legacy-tag guard and the parameter
        # mismatch check above have both had their chance to abort. Writing any earlier
        # would leave a sidecar behind even on a run that aborts, which a later --resume
        # would misread as "parameters differ" instead of reporting the true cause.
        write_params_sidecar(args.output, current_params)
    else:
        if os.path.exists(args.output) and os.path.getsize(args.output) > 0:
            with open(args.output) as fh:
                n_lines = sum(1 for _ in fh)
            warn(f"WARNING: {args.output} exists and holds {n_lines} line"
                 f"{'s' if n_lines != 1 else ''}; truncating it because --resume was not "
                 f"given. Pass --resume to keep existing results instead.")
        # Truncate. synmap_split appends per segment, so without this a re-run duplicates
        # every already-aligned row.
        open(args.output, "w").close()
        with open(args.skipped, "w") as fh:
            fh.write("coord1\tcoord2\tstrand\treason\tratio\tlen1\tlen2\n")
        write_params_sidecar(args.output, align_params(args, timer))

    total = len(filtered)
    if total == 0:
        print("Nothing to align.")
        return

    print(f"Aligning {total} segments with {args.aligner} (preset {args.preset}, "
          f"max_len_ratio {args.max_len_ratio})")

    temp_base = os.path.abspath("./align_temp")
    os.makedirs(temp_base, exist_ok=True)

    opts = {
        "aligner": args.aligner, "preset": args.preset, "threads": args.threads,
        "timer": timer, "debug": args.debug, "temp_base": temp_base,
        "output_file": args.output, "minimap2_bin": args.minimap2_bin,
        "max_len_ratio": args.max_len_ratio,
        "wfa_bin": args.wfa_bin, "wfa_model": args.wfa_model, "wfa_span": args.wfa_span,
        "wfa_memory": args.wfa_memory, "wfa_heuristic": args.wfa_heuristic,
        "wfa_mismatch": args.wfa_mismatch, "wfa_gap_open1": args.wfa_gap_open1,
        "wfa_gap_ext1": args.wfa_gap_ext1, "wfa_gap_open2": args.wfa_gap_open2,
        "wfa_gap_ext2": args.wfa_gap_ext2, "wfa_chain": args.wfa_chain,
        "min_chain_len": args.min_chain_len, "min_anchor_len": args.min_anchor_len,
        "last_matrix": args.last_matrix, "last_gap_open": args.last_gap_open,
        "last_gap_ext": args.last_gap_ext, "last_near": args.last_near,
        "last_split": args.last_split, "last_threads": args.threads,
    }

    manager = Manager()
    counter = manager.Value("i", 0)

    stop = threading.Event()

    def reporter():
        i = 1
        while not stop.is_set() and counter.value < total:
            time.sleep(60)
            if counter.value >= total:
                break
            print(f"{i} minute: {int(100 * counter.value / total)}% complete")
            i += 1

    th = threading.Thread(target=reporter, daemon=True)
    th.start()

    baseline_skip_rows = len(_read_skip_rows(args.skipped))

    tasks = [(line, args.genome_dir, opts, counter, args.skipped) for line in filtered]
    with Pool(processes=args.processes) as pool:
        pool.map(process_line, tasks)
    stop.set()

    # Re-read rather than count in-process: skips are recorded straight to disk by pool
    # workers (see record_skip), so the sidecar is authoritative. On --resume it also holds
    # rows from earlier invocations, so only the tail added since baseline_skip_rows belongs
    # to this run -- reporting the raw file total there would misattribute historical skips.
    all_skip_rows = _read_skip_rows(args.skipped)
    this_run_rows = all_skip_rows[baseline_skip_rows:]
    if this_run_rows:
        from collections import Counter as _C
        reasons = _C(r[3] for r in this_run_rows)
        print(f"Skipped {len(this_run_rows)} segments this run -> {args.skipped}")
        for reason, n in reasons.most_common():
            print(f"    {reason}: {n}")
        if len(all_skip_rows) != len(this_run_rows):
            print(f"    ({len(all_skip_rows)} total skipped across all runs recorded in "
                  f"{args.skipped})")

    if not args.debug:
        try:
            os.rmdir(temp_base)
        except OSError:
            pass


if __name__ == "__main__":
    main()
