#!/usr/bin/env python3
"""Extract syntenic segment pairs and align them, emitting contract-compliant PAF.

Backends are selected with --aligner. Each one only produces alignments; all coordinate
mapping and tag construction live in paf_emit.py, so the arithmetic is written and tested
once instead of once per backend. See step10_dev/memo_step10_diagnosis.md for the audit
behind that split and for the PAF contract Steps 11-17 depend on.
"""

import argparse
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


def parse_syncoord(syncoord):
    """'AthaAt4_chr4:11535683..11542200' -> ('AthaAt4', 'AthaAt4_chr4', 11535683, 11542200).

    Coordinates are 1-based inclusive, as written by gene_coords_extractor_all4.py.
    """
    name, rng = syncoord.split(":")
    start, end = (int(x) for x in rng.split(".."))
    accession = name.split("_")[0]
    return accession, name, start, end


def extract(genome_dir, syncoord, do_revcomp=False, upper=True):
    accession, name, start, end = parse_syncoord(syncoord)
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


def build_segment(syncoord1, syncoord2, q_revcomp):
    _, t_name, t_start, t_end = parse_syncoord(syncoord1)
    _, q_name, q_start, q_end = parse_syncoord(syncoord2)
    return paf_emit.Segment(t_name, t_start, t_end, q_name, q_start, q_end, q_revcomp)


def length_ratio(seg):
    lo, hi = sorted((seg.t_len, seg.q_len))
    return float("inf") if lo <= 0 else hi / lo


def process_line(args):
    line, genome_dir, opts, counter, skipped = args
    try:
        syncoord1, syncoord2, strand = line.strip().split("\t")
    except ValueError:
        warn(f"malformed coords line, skipping: {line.strip()!r}")
        counter.value += 1
        return

    try:
        seg_probe = build_segment(syncoord1, syncoord2, False)

        # Guard before doing any work. A 242:1 span is not alignable end-to-end; riparian's
        # rescue machinery deliberately keeps such blocks for plot completeness, which is
        # the opposite of what an aligner wants, so only its ratio test is borrowed.
        if seg_probe.t_len <= 0 or seg_probe.q_len <= 0:
            skipped.append(f"{syncoord1}\t{syncoord2}\t{strand}\tempty_span\t0\t"
                           f"{seg_probe.t_len}\t{seg_probe.q_len}")
            counter.value += 1
            return
        ratio = length_ratio(seg_probe)
        if opts["max_len_ratio"] > 0 and ratio > opts["max_len_ratio"]:
            skipped.append(f"{syncoord1}\t{syncoord2}\t{strand}\tlength_ratio\t{ratio:.1f}\t"
                           f"{seg_probe.t_len}\t{seg_probe.q_len}")
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
                skipped.append(f"{syncoord1}\t{syncoord2}\t{strand}\tno_orientation\t"
                               f"{ratio:.1f}\t{seg_probe.t_len}\t{seg_probe.q_len}")
                counter.value += 1
                return
            q_revcomp = needs_rc

        seg = build_segment(syncoord1, syncoord2, q_revcomp)
        qseq = extract(genome_dir, syncoord2, do_revcomp=q_revcomp)

        try:
            lines = BACKENDS[opts["aligner"]](tseq, qseq, seg, opts)
        except subprocess.TimeoutExpired:
            skipped.append(f"{syncoord1}\t{syncoord2}\t{strand}\ttimeout\t{ratio:.1f}\t"
                           f"{seg_probe.t_len}\t{seg_probe.q_len}")
            counter.value += 1
            return

        if lines:
            with open(opts["output_file"], "a") as fh:
                fh.write("\n".join(lines) + "\n")
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
        acc1 = syn1.split(":")[0].split("_")[0]
        acc2 = syn2.split(":")[0].split("_")[0]
        if args.pairedIDs:
            if acc1 in args.pairedIDs and acc2 in args.pairedIDs:
                filtered.append(line)
        elif args.singleID:
            if args.singleID in (acc1, acc2):
                filtered.append(line)
        else:
            filtered.append(line)
    total = len(filtered)
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
    skipped = manager.list()

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

    tasks = [(line, args.genome_dir, opts, counter, skipped) for line in filtered]
    with Pool(processes=args.processes) as pool:
        pool.map(process_line, tasks)
    stop.set()

    if skipped:
        with open(args.skipped, "w") as fh:
            fh.write("coord1\tcoord2\tstrand\treason\tratio\tlen1\tlen2\n")
            for row in skipped:
                fh.write(row + "\n")
        from collections import Counter as _C
        reasons = _C(r.split("\t")[3] for r in skipped)
        print(f"Skipped {len(skipped)}/{total} segments -> {args.skipped}")
        for reason, n in reasons.most_common():
            print(f"    {reason}: {n}")

    if not args.debug:
        try:
            os.rmdir(temp_base)
        except OSError:
            pass


if __name__ == "__main__":
    main()
