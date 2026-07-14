#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import re
import tempfile
import shutil
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from fastaio import fasta_stem, is_compressed, materialize_plain

# python liftover.py --genome Bdhap.fa Cdact.fa Eindi.fa --reference Sbico.pep --outn 1 --outs 0.99 --outc 0.9 --threads 20
#
# Two approaches to generate CDS (1: concatenating 'CDS' feature from gff; 2: Using the '##ATN' sequence).
# The only difference between the two is that concatenating 'CDS' seqs SOMETIMES omits a stop codon whereas using '##ATN' seq retains all stop codons.
# Here, I opted to retain the '##ATN' approach since it seems more consistent with the peptide sequence with "##STA".
# The inframe output enables running KaKs_Calculator (but first run ParaAT to generate AXT files). Then can probably feed kaks, synteny (jcvi), and orthofinder results to SOI dotplot.

def ensure_miniprot():
    """
    Return a miniprot binary, preferring one already installed.

    Order: PATH (e.g. the conda/mamba env) -> a local ./miniprot clone -> clone and build.
    PATH comes first because it is pinned and reproducible, and because cloning needs
    network access, git and a compiler at runtime -- none of which a compute node is
    guaranteed to have.
    """
    found = shutil.which("miniprot")
    if found:
        return found

    binary = os.path.join("miniprot", "miniprot")
    if os.path.isfile(binary):
        return binary

    print("miniprot not found on PATH, cloning and building...", file=sys.stderr)
    subprocess.check_call(["git", "clone", "https://github.com/lh3/miniprot"])
    subprocess.check_call(["make"], cwd="miniprot")
    return binary

def run_miniprot(genome, prot, threads, outn, outs, outc, write_gff_path):
    """
    Run miniprot on a single genome with the given protein file and parameters.
    write_gff_path: path to write GFF (may be a temp file if user didn't request gff)
    """
    cmd = [
        ensure_miniprot(),
        "--trans",
        "--gff",
        "--aln",
        "-t", str(threads),
        genome,
        prot,
        "-P", os.path.splitext(os.path.basename(genome))[0]
    ]
    if outn is not None:
        cmd += ["--outn", str(outn)]
    if outs is not None:
        cmd += ["--outs", str(outs)]
    if outc is not None:
        cmd += ["--outc", str(outc)]

    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    with open(write_gff_path, "w") as out:
        subprocess.check_call(cmd, stdout=out)

def parse_gff(gff_file):
    """
    Parse the miniprot GFF and return:
      • peptide translations           → pep_data
      • mRNA coords (BED6)             → mrna_bed
      • CDS coords (BED6)              → cds_bed
      • wobble-frame CDS sequences     → wobbleframe_data

    New behavior:
      - If multiple mRNAs share the same Target, they are renamed to share the
        first-seen mRNA ID as a common prefix with numbered suffixes:
        <first_tid>_1, <first_tid>_2, ...
      - Singletons keep their original ID.
      - All outputs (pep, mRNA bed/seq, cds bed) use the renamed IDs.
    """
    import re
    # First pass: collect mRNA entries and CDS entries
    # We need the sequences that precede each mRNA (##STA / ##ATN).
    current_sta = ""
    current_atn = ""

    # mRNAs: tid -> dict with target, coords, strand, sequences
    mrnas = {}
    # Target -> list of tids in encounter order
    target_to_tids = defaultdict(list)
    # Keep encounter order of tids globally (if needed)
    tid_order = []

    # Collect CDS entries as (chrom, start0, end0, strand, parent_tid)
    cds_entries = []

    target_pat = re.compile(r"Target=([^;\s]+)")
    id_pat = re.compile(r"ID=([^;\s]+)")
    parent_pat = re.compile(r"Parent=([^;\s]+)")

    with open(gff_file) as fh:
        for raw in fh:
            line = raw.rstrip()

            if line.startswith("#"):
                if line.startswith("##STA"):
                    # Amino-acid translation (for >pep)
                    # After tab, everything is the sequence
                    current_sta = line.split("\t", 1)[1]
                elif line.startswith("##ATN"):
                    # Nucleotide CDS with codon frames encoded; clean as in your code
                    raw_atn = line.split("\t", 1)[1]
                    clean = re.sub(r"~\d+~", "", raw_atn)
                    clean = clean.replace("-", "")
                    clean = re.sub(r"[a-z]", "", clean)
                    current_atn = clean
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            feat = cols[2]
            attrs = cols[8]

            if feat == "mRNA":
                m = id_pat.search(attrs)
                if not m:
                    continue
                tid = m.group(1)

                tmatch = target_pat.search(attrs)
                target_id = tmatch.group(1) if tmatch else "__NO_TARGET__"

                s0 = int(start) - 1
                e0 = int(end)

                # Store info for this mRNA (sequence state at this time)
                mrnas[tid] = {
                    "target": target_id,
                    "chrom": chrom,
                    "start0": s0,
                    "end0": e0,
                    "strand": strand,
                    "sta": current_sta,
                    "atn": current_atn,
                }
                target_to_tids[target_id].append(tid)
                tid_order.append(tid)

            elif feat == "CDS":
                pm = parent_pat.search(attrs)
                if not pm:
                    continue
                parent_tid = pm.group(1)
                s0 = int(start) - 1
                e0 = int(end)
                cds_entries.append((chrom, s0, e0, strand, parent_tid))

    # Build rename map:
    # - singleton target -> original id
    # - multi-target -> first seen tid as prefix, numbered by encounter order
    rename_map = {}
    for target, tids in target_to_tids.items():
        if len(tids) == 1:
            rename_map[tids[0]] = tids[0]
        else:
            base = tids[0]  # first-seen for this Target
            for i, old_tid in enumerate(tids, start=1):
                rename_map[old_tid] = f"{base}_{i}"

    # Construct outputs using renamed IDs
    pep_data = []
    mrna_bed = []
    cds_bed = []
    wobbleframe_data = []

    for tid in tid_order:
        if tid not in mrnas:
            continue
        rec = mrnas[tid]
        new_id = rename_map.get(tid, tid)

        # Peptide FASTA from ##STA
        if rec["sta"]:
            pep_data.append(f">{new_id}\n{rec['sta']}")

        # Wobble-frame (ATN) CDS sequence
        if rec["atn"]:
            wobbleframe_data.append(f">{new_id}\n{rec['atn']}")

        # mRNA BED
        mrna_bed.append(
            f"{rec['chrom']}\t{rec['start0']}\t{rec['end0']}\t{new_id}\t0\t{rec['strand']}"
        )

    # CDS BED with Parent renamed
    for chrom, s0, e0, strand, parent_tid in cds_entries:
        new_parent = rename_map.get(parent_tid, parent_tid)
        cds_bed.append(f"{chrom}\t{s0}\t{e0}\t{new_parent}\t0\t{strand}")

    return pep_data, mrna_bed, cds_bed, wobbleframe_data

def write_output(pep_data, mrna_bed, cds_bed, wobbleframe_data, prefix, genome,
                 walker_path, threads, outputs):
    """
    Write only requested outputs. If an intermediate is needed to build a requested output,
    create it temporarily and remove it if it wasn't explicitly requested.
    """
    want = set(outputs) if outputs else {"gff", "pep", "bed", "fai", "mRNA", "cds", "inframe"}  # default excludes TEsorter artifacts
    created_to_cleanup = []

    # 1) Peptide FASTA
    pep_file = f"{prefix}.pep"
    if pep_data and ("pep" in want or "inframe" in want):
        if "pep" in want:
            if not os.path.exists(pep_file):
                with open(pep_file, "w") as fh:
                    fh.write("\n".join(pep_data) + "\n")
        else:
            pep_tmp = tempfile.NamedTemporaryFile(delete=False, prefix=f"{prefix}.", suffix=".pep.tmp")
            pep_tmp.write(("\n".join(pep_data) + "\n").encode())
            pep_tmp.close()
            pep_file = pep_tmp.name
            created_to_cleanup.append(pep_file)

    # 2) mRNA coordinates + sequences
    bed_file = f"{prefix}.bed"
    mfa = f"{prefix}.mRNA"
    bed_needed = ("bed" in want) or ("mRNA" in want)
    if mrna_bed and bed_needed:
        if "bed" in want:
            if not os.path.exists(bed_file):
                with open(bed_file, "w") as fh:
                    fh.write("\n".join(mrna_bed) + "\n")
        else:
            bed_tmp = tempfile.NamedTemporaryFile(delete=False, prefix=f"{prefix}.", suffix=".bed.tmp")
            bed_tmp.write(("\n".join(mrna_bed) + "\n").encode())
            bed_tmp.close()
            bed_file = bed_tmp.name
            created_to_cleanup.append(bed_file)

        if "mRNA" in want:
            if not os.path.exists(mfa):
                print(f"Extracting mRNA to {mfa}", file=sys.stderr)
                raw = subprocess.check_output([
                    "bedtools", "getfasta", "-s",
                    "-fi", genome,
                    "-bed", bed_file,
                    "-name+"
                ]).decode()
                with open(mfa, "w") as out:
                    name = None
                    seq_parts = []
                    for L in raw.splitlines():
                        if L.startswith(">"):
                            if name:
                                out.write(f">{name}\n{''.join(seq_parts)}\n")
                            header = L[1:].split()[0]
                            name = header.split("::")[0]
                            seq_parts = []
                        else:
                            seq_parts.append(L.strip())
                    if name:
                        out.write(f">{name}\n{''.join(seq_parts)}\n")

    # 3) CDS FASTA (wobble-frame, retains stop codon)
    cds_file = f"{prefix}.cds"
    if wobbleframe_data and ("cds" in want or "inframe" in want):
        if "cds" in want:
            if not os.path.exists(cds_file):
                with open(cds_file, "w") as fh:
                    fh.write("\n".join(wobbleframe_data) + "\n")
        else:
            cds_tmp = tempfile.NamedTemporaryFile(delete=False, prefix=f"{prefix}.", suffix=".cds.tmp")
            cds_tmp.write(("\n".join(wobbleframe_data) + "\n").encode())
            cds_tmp.close()
            cds_file = cds_tmp.name
            created_to_cleanup.append(cds_file)

    # 4) CDS in-frame via cds_walker.py
    if "inframe" in want:
        inframe_file = f"{prefix}.cds.inframe"
        walker_script = os.path.join(walker_path, "cds_walker.py")
        if os.path.isfile(walker_script):
            print(f"Running in-frame CDS: python {walker_script} "
                  f"-c {cds_file} -p {pep_file} -o {inframe_file} -t {threads}", file=sys.stderr)
            subprocess.check_call([
                sys.executable, walker_script,
                "-c", cds_file,
                "-p", pep_file,
                "-o", inframe_file,
                "-t", str(threads)
            ])
        else:
            print(f"WARNING, cds_walker.py not found at {walker_script}. "
                  f"Skipping building {inframe_file}.", file=sys.stderr)

    # Clean up intermediates not requested explicitly
    for path in created_to_cleanup:
        try:
            os.remove(path)
        except OSError:
            pass

def maybe_faidx(genome, want):
    if "fai" in want:
        print(f"Indexing FASTA (samtools faidx): {genome}", file=sys.stderr)
        subprocess.check_call(["samtools", "faidx", genome])

def run_cdhit_on_proteins(reference_fa):
    """
    Run cd-hit on proteins in current directory using basename for outputs.
    Returns path to deduplicated FASTA.
    """
    base = os.path.splitext(os.path.basename(reference_fa))[0]
    cleaned = base + ".cdhit"
    print(f"Running cd-hit on proteins: cd-hit -i {reference_fa} -o {cleaned}", file=sys.stderr)
    subprocess.check_call(["cd-hit", "-i", reference_fa, "-o", cleaned])
    return cleaned  # without extension: cd-hit creates exactly this file

def count_fasta_records(path):
    """Number of '>' headers in a FASTA, or 0 if it is missing/empty."""
    if not path or not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def copy_fasta(src, dst):
    """Materialise dst as a copy of src (used when a filtering pass is a no-op)."""
    shutil.copyfile(src, dst)
    return dst


def seqkit_exclude(ids_file, in_fa, out_fa):
    """
    Drop the sequences named in ids_file from in_fa, writing out_fa.

    An empty ids_file means 'exclude nothing'. seqkit would technically do the right
    thing there (it warns '0 patterns loaded' and passes everything through), but
    calling it is pointless and the warning reads like a failure. Copy instead.
    """
    if os.path.getsize(ids_file) == 0:
        return copy_fasta(in_fa, out_fa)
    with open(out_fa, "w") as out:
        subprocess.check_call(["seqkit", "grep", "-v", "-f", ids_file, in_fa], stdout=out)
    return out_fa


def run_tesorter_two_pass(protein_fa, threads):
    """
    Two-pass TE screen of the reference proteome: TEsorter (HMM) then BLASTP against
    the TE peptides TEsorter found. Artifacts land in CWD, prefixed by the basename
    of protein_fa. Returns (path to the cleaned FASTA, artifact manifest).

    A proteome with no detectable TEs is a normal outcome, not an error. Curated
    reference annotations (Ensembl's Amborella pep.all, for one) are already
    TE-filtered upstream, so TEsorter classifies nothing and the '.rexdb.cls.pep'
    it writes is legitimately empty. Both passes below degrade to no-ops in that
    case: an empty TE set means there is nothing to build a BLAST DB from and
    nothing to remove, so the proteome passes through unchanged. Feeding the empty
    file to makeblastdb is what used to abort the run.
    """
    prefix = os.path.basename(protein_fa)  # e.g. Sbico.pep.cdhit -- keep the dots

    final_no_tes = f"{prefix}_no_TEs.faa"
    rnd1 = f"{prefix}_rnd1.TE_peps.faa"
    te_list = f"{prefix}.TE_peps.list"
    ids_for_seqkit = f"{prefix}.TE_peps.ids"
    dummy_faa = f"{prefix}.TE_peptides.dummy.faa"
    rnd2_list = f"{prefix}.rnd2.TE_peps.list"
    rnd2_ids = f"{prefix}.rnd2.TE_peps.ids"
    cls_tsv = f"{prefix}.rexdb.cls.tsv"
    cls_pep = f"{prefix}.rexdb.cls.pep"
    db_dir = "db"
    db_prefix = os.path.join(db_dir, "TEpep")

    artifacts = {
        "keep_if_outputs_TEsorter": [
            final_no_tes, rnd2_list, db_dir, dummy_faa, rnd1, te_list,
            cls_pep,
            f"{prefix}.rexdb.cls.lib",
            cls_tsv,
            f"{prefix}.rexdb.dom.tsv",
            f"{prefix}.rexdb.dom.faa",
            f"{prefix}.rexdb.dom.gff3",
            f"{prefix}.rexdb.domtbl",
        ],
        "extra_tmp": [ids_for_seqkit, rnd2_ids],
    }

    n_in = count_fasta_records(protein_fa)

    # Resume on the finished product only. Intermediates are not a safe resume point:
    # TEsorter reuses an existing '.domtbl' and skips hmmscan whenever the file is
    # merely non-empty, so a truncated domtbl from a killed run would silently
    # under-call TEs. '-fw' below forces the scan; this check is what keeps a
    # completed run from paying for it twice.
    if os.path.exists(final_no_tes) and os.path.getsize(final_no_tes) > 0:
        print(f"TE screen: reusing existing {final_no_tes} "
              f"({count_fasta_records(final_no_tes)} peptides).", file=sys.stderr)
        return final_no_tes, artifacts

    # 1) TEsorter (protein mode)
    print(f"Running: TEsorter {protein_fa} -st prot -p {threads} -fw", file=sys.stderr)
    subprocess.check_call(["TEsorter", protein_fa, "-st", "prot", "-p", str(threads), "-fw"])

    # 2) First pass list: first column of cls.tsv, minus comments and any leading '>'.
    #    TEsorter still writes cls.tsv when it classifies nothing, but do not depend on it.
    with open(te_list, "w") as out_list, open(ids_for_seqkit, "w") as out_ids:
        if os.path.exists(cls_tsv):
            with open(cls_tsv) as inp:
                for line in inp:
                    if line.startswith("#") or not line.strip():
                        continue
                    first = line.split("\t", 1)[0].strip()
                    out_list.write(first + "\n")
                    out_ids.write(first.lstrip(">") + "\n")

    with open(ids_for_seqkit) as fh:
        n_hmm = sum(1 for _ in fh)

    # 3) Remove HMM-labelled TE peptides (pass 1)
    print(f"Removing TE peptides (pass1, {n_hmm} hit) → {rnd1}", file=sys.stderr)
    seqkit_exclude(ids_for_seqkit, protein_fa, rnd1)

    # 4) Build the TE peptide BLAST DB from TEsorter's classified domains.
    #    With nothing classified there is no DB to build and no second pass to run.
    n_te_pep = count_fasta_records(cls_pep)
    if n_te_pep == 0:
        print(
            "TE screen: TEsorter classified 0 peptides as TE-derived, so the BLASTP "
            "second pass has no database to search and is skipped. The proteome is "
            "passed through unchanged. This is expected for curated reference "
            "proteomes that were already TE-filtered by their annotation pipeline; "
            "if you expected TEs here, check that the TEsorter database matches the "
            "clade (-db rexdb-plant for plants).",
            file=sys.stderr,
        )
        open(dummy_faa, "w").close()
        open(rnd2_list, "w").close()
        open(rnd2_ids, "w").close()
        copy_fasta(rnd1, final_no_tes)
        print(f"TE screen: {n_in} in → {n_in} out (0 removed).", file=sys.stderr)
        return final_no_tes, artifacts

    # Every peptide being a TE means the screen ate the proteome. miniprot would then
    # align nothing and the pipeline would fail several expensive steps later, so stop here.
    if count_fasta_records(rnd1) == 0:
        raise SystemExit(
            f"[FATAL] The TE screen removed all {n_in} reference peptides "
            f"({protein_fa}). Nothing is left to lift over. Check that the reference "
            f"is a proteome (not a TE library) and that the TEsorter database is right."
        )

    print(f"Building dummy TE peptide DB input ({n_te_pep} peptides) → {dummy_faa}", file=sys.stderr)
    with open(dummy_faa, "w") as out, open(cls_pep) as inp:
        i = 0
        for line in inp:
            if line.startswith(">"):
                i += 1
                out.write(f">TE{i}\n")
            else:
                out.write(line)

    os.makedirs(db_dir, exist_ok=True)
    print(f"makeblastdb -in {dummy_faa} -dbtype prot -parse_seqids -title TEpep -out {db_prefix}", file=sys.stderr)
    subprocess.check_call([
        "makeblastdb", "-in", dummy_faa, "-dbtype", "prot",
        "-parse_seqids", "-title", "TEpep", "-out", db_prefix
    ])

    # 5) blastp second pass: catch TE peptides the HMMs missed but that resemble the ones they caught.
    print(f"BLASTP second pass → {rnd2_list}", file=sys.stderr)
    with open(rnd2_list, "w") as out:
        subprocess.check_call([
            "blastp",
            "-query", rnd1,
            "-db", db_prefix,
            "-evalue", "1e-5",
            "-max_target_seqs", "5",
            "-max_hsps", "1",
            "-qcov_hsp_perc", "50",
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs scovhsp",
            "-num_threads", str(threads)
        ], stdout=out)

    # Unique qseqids -> the peptides to drop. No hits is fine: nothing to drop.
    seen = set()
    with open(rnd2_ids, "w") as out, open(rnd2_list) as inp:
        for line in inp:
            if not line.strip():
                continue
            qid = line.split("\t", 1)[0].strip()
            if qid not in seen:
                seen.add(qid)
                out.write(qid + "\n")

    print(f"Removing BLAST-matched TE peptides (pass2, {len(seen)} hit) → {final_no_tes}", file=sys.stderr)
    seqkit_exclude(rnd2_ids, rnd1, final_no_tes)

    n_out = count_fasta_records(final_no_tes)
    if n_out == 0:
        raise SystemExit(
            f"[FATAL] The TE screen removed all {n_in} reference peptides "
            f"({protein_fa}). Nothing is left to lift over."
        )
    print(f"TE screen: {n_in} in → {n_out} out "
          f"({n_hmm} removed by HMM, {len(seen)} by BLASTP).", file=sys.stderr)

    return final_no_tes, artifacts


def main():
    parser = argparse.ArgumentParser(
        description="Protein liftover using miniprot + generate selected outputs (.gff, .pep, .bed, .fai, .mRNA, .cds, .inframe, .cdhit, .TEsorter)"
    )
    parser.add_argument("--genome", "-g", nargs='+', required=True,
                        help="One or more genome FASTAs")
    parser.add_argument("--walker-path", default=None,
                        help="Directory containing cds_walker.py (default: script dir)")
    parser.add_argument("--reference", "-r", required=True,
                        help="Reference protein FASTA")
    parser.add_argument("--threads", "-t", type=int, default=20,
                        help="Threads for miniprot and TEsorter/BLAST")
    parser.add_argument("--outn", type=int, default=None,
                        help="--outn for miniprot")
    parser.add_argument("--outs", type=float, default=None,
                        help="--outs for miniprot")
    parser.add_argument("--outc", type=float, default=None,
                        help="--outc for miniprot")

    parser.add_argument("--cdhit", action="store_true",
                        help="Run 'cd-hit -i <prot> -o <prot>.cdhit' and use the deduplicated proteins")
    parser.add_argument("--TEsorter", action="store_true",
                        help="Run TEsorter two-pass cleaning on the (optionally cd-hit deduplicated) protein FASTA")

    parser.add_argument("--outputs", nargs='+',
                        choices=["gff", "pep", "bed", "fai", "mRNA", "cds", "inframe", "cdhit", "TEsorter"],
                        help="Which outputs to produce (choose one or more). Default: gff pep bed fai mRNA cds inframe cdhit")

    args = parser.parse_args()

    # Determine walker script path
    base_walker = args.walker_path or os.path.dirname(os.path.realpath(__file__))
    # Default outputs exclude TEsorter artifacts unless explicitly asked
    default_outputs = {"gff", "pep", "bed", "fai", "mRNA", "cds", "inframe", "cdhit"}
    want = set(args.outputs) if args.outputs else default_outputs

    # Compressed inputs: cd-hit, TEsorter, makeblastdb, samtools faidx and bedtools
    # getfasta all read plain FASTA only, so decompress once up front rather than
    # teaching each call site. Uncompressed inputs are passed through untouched.
    def _log(msg):
        print(msg, file=sys.stderr)

    reference = materialize_plain(args.reference, log=_log)
    genomes = [materialize_plain(g, log=_log) for g in args.genome]

    # Prepare protein FASTA (optionally de-duplicated with cd-hit, then optionally TEsorter two-pass)
    # Use basename to ensure artifacts are created in CWD, not alongside the source reference
    prot_for_miniprot = reference
    created_cdhit_base = None
    tesorter_artifacts = None

    if args.cdhit:
        cleaned = run_cdhit_on_proteins(reference)
        prot_for_miniprot = cleaned
        created_cdhit_base = cleaned

    if args.TEsorter:
        # Run TEsorter on whichever file is the current "prot_for_miniprot"
        final_no_tes, artifacts = run_tesorter_two_pass(prot_for_miniprot, args.threads)
        prot_for_miniprot = final_no_tes
        tesorter_artifacts = artifacts

    for genome in genomes:
        prefix = fasta_stem(genome)

        if "fai" in want:
            print(f"Indexing FASTA (samtools faidx): {genome}", file=sys.stderr)
            subprocess.check_call(["samtools", "faidx", genome])

        need_miniprot = bool({"gff", "pep", "bed", "mRNA", "cds", "inframe"} & want)
        if not need_miniprot:
            continue

        gff_path = f"{prefix}.gff"
        created_gff_this_run = False

        if os.path.exists(gff_path):
            pep, mrna, cds, wobbleframe = parse_gff(gff_path)
            write_output(pep, mrna, cds, wobbleframe, prefix, genome, base_walker, args.threads, want)
        else:
            if "gff" in want:
                run_miniprot(genome, prot_for_miniprot, args.threads, args.outn, args.outs, args.outc, gff_path)
                created_gff_this_run = True
                pep, mrna, cds, wobbleframe = parse_gff(gff_path)
                write_output(pep, mrna, cds, wobbleframe, prefix, genome, base_walker, args.threads, want)
            else:
                with tempfile.NamedTemporaryFile(delete=False, prefix=f"{prefix}.", suffix=".gff.tmp") as tmpgff:
                    tmp_gff_path = tmpgff.name
                try:
                    run_miniprot(genome, prot_for_miniprot, args.threads, args.outn, args.outs, args.outc, tmp_gff_path)
                    pep, mrna, cds, wobbleframe = parse_gff(tmp_gff_path)
                    write_output(pep, mrna, cds, wobbleframe, prefix, genome, base_walker, args.threads, want)
                finally:
                    try:
                        os.remove(tmp_gff_path)
                    except OSError:
                        pass

        if created_gff_this_run and "gff" not in want:
            try:
                os.remove(gff_path)
            except OSError:
                pass

    # Post-run cleanup policy

    # If we ran cd-hit but 'cdhit' wasn't requested as an output, remove its artifacts
    if created_cdhit_base and "cdhit" not in want:
        for suffix in ["", ".clstr"]:
            try:
                os.remove(created_cdhit_base + suffix)
            except OSError:
                pass

    # If we ran TEsorter but 'TEsorter' isn't in outputs, remove TEsorter artifacts
    if tesorter_artifacts and "TEsorter" not in want:
        for f in tesorter_artifacts["keep_if_outputs_TEsorter"]:
            if f == "db":
                # remove db directory
                try:
                    shutil.rmtree("db")
                except OSError:
                    pass
            else:
                try:
                    os.remove(f)
                except OSError:
                    pass
        # always drop little tmp id lists
        for f in tesorter_artifacts["extra_tmp"]:
            try:
                os.remove(f)
            except OSError:
                pass

if __name__ == "__main__":
    main()
