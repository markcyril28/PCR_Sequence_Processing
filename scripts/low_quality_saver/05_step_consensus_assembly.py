#!/usr/bin/env python3
"""Method 2: Reference-anchored consensus assembly (no trace re-editing).

Maps the KB-basecalled forward + reverse reads to the assigned reference
transcript and writes per-strand reference-coordinate sequences (with quality)
that the downstream consensus step can vote on.

Outputs (per read) under <out-dir>/<sample>/:
    <stem>.aligned.fa     read sequence projected onto reference coordinates
    <stem>.aligned.qual   per-ref-position quality (0 where the read had a gap)
    <stem>.ref_align.txt  Tracy-style visual alignment of read vs reference
"""
from __future__ import annotations

import argparse
import importlib
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
format_alignment = importlib.import_module(
    "scripts.common.02_alignment_viz"
).format_alignment
mask_trim = importlib.import_module("scripts.common.03_trim_reads").mask_trim


def _read_assignments(path: Path) -> dict[tuple[str, str], str]:
    out = {}
    with path.open() as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2] not in ("NO_HIT", "FAIL_THRESHOLD", "-", ""):
                out[(p[0], p[1])] = p[2]
    return out


def _project(ref: str, read: str, quals: list[int],
             gap_open: float, gap_extend: float
             ) -> tuple[str, list[int], str, str, float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2; aligner.mismatch_score = -1
    aligner.open_gap_score = gap_open; aligner.extend_gap_score = gap_extend
    aln = aligner.align(ref, read)[0]
    aref, aread = str(aln[0]), str(aln[1])

    bases, qs = [], []
    ri = 0
    for ref_b, read_b in zip(aref, aread):
        if ref_b == "-":
            # Insertion in read relative to reference: consume read base, drop it.
            if read_b != "-":
                ri += 1
            continue
        if read_b == "-":
            bases.append("-")
            qs.append(0)
        else:
            bases.append(read_b)
            qs.append(quals[ri] if ri < len(quals) else 0)
            ri += 1
    return "".join(bases), qs, aref, aread, float(aln.score)


def _orient_to_reference(ref: str, seq: str, quals: list[int],
                         gap_open: float, gap_extend: float
                         ) -> tuple[str, list[int], str, str, str, float, float]:
    """Pick the read orientation (as-is vs reverse-complement) that aligns
    better to ``ref``. Returns (proj, pq, aref, aread, orientation,
    score_as_is, score_revcomp). Score-based; robust to mislabeled files.
    """
    proj_fw, pq_fw, aref_fw, aread_fw, score_fw = _project(
        ref, seq, quals, gap_open, gap_extend)
    rc_seq = str(Seq(seq).reverse_complement())
    rc_quals = list(reversed(quals))
    proj_rc, pq_rc, aref_rc, aread_rc, score_rc = _project(
        ref, rc_seq, rc_quals, gap_open, gap_extend)
    if score_rc > score_fw:
        return proj_rc, pq_rc, aref_rc, aread_rc, "revcomp", score_fw, score_rc
    return proj_fw, pq_fw, aref_fw, aread_fw, "as_is", score_fw, score_rc


def _process_one(sample: str, ab1_str: str, ref_id: str, ref_seq: str,
                 out_dir_str: str, gap_open: float, gap_extend: float,
                 trim_q: int = 0, trim_min_len: int = 50) -> str:
    ab1 = Path(ab1_str)
    stem = ab1.stem
    rec = next(SeqIO.parse(str(ab1), "abi"), None)
    if rec is None:
        return f"  WARN: no records in ABI file: {ab1}, skipping"
    seq = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality", [])) or [0] * len(seq)
    trim_tag = ""
    if trim_q > 0:
        seq, quals, t_lo, t_hi = mask_trim(seq, quals,
                                           q_cutoff=trim_q,
                                           min_len=trim_min_len)
        trim_tag = f" trim=[{t_lo},{t_hi})"
    proj, pq, aref_str, aread_str, orientation, score_fw, score_rc = \
        _orient_to_reference(ref_seq, seq, quals, gap_open, gap_extend)
    out_sample = Path(out_dir_str) / sample
    out_sample.mkdir(parents=True, exist_ok=True)
    (out_sample / f"{stem}.aligned.fa").write_text(
        f">{stem} ref={ref_id} projected=1 orientation={orientation}\n" +
        "\n".join(proj[i:i+60] for i in range(0, len(proj), 60)) + "\n")
    (out_sample / f"{stem}.aligned.qual").write_text(
        " ".join(str(q) for q in pq) + "\n")
    (out_sample / f"{stem}.ref_align.txt").write_text(
        f"# orientation={orientation}  "
        f"score_as_is={score_fw:.1f}  score_revcomp={score_rc:.1f}\n"
        + format_alignment("PCR_Seq", aread_str, "Reference", aref_str))
    return (f"  {sample}/{stem}  ref={ref_id} len={len(proj)} "
            f"orient={orientation}{trim_tag}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-root", required=True, type=Path)
    ap.add_argument("--reference",   required=True, type=Path)
    ap.add_argument("--assignments", required=True, type=Path)
    ap.add_argument("--out-dir",     required=True, type=Path)
    ap.add_argument("--aligner",      default="biopython")
    ap.add_argument("--min-quality",  type=int,   default=20)  # informational only here
    ap.add_argument("--disagreement", default="higher_q")       # used by merge step
    ap.add_argument("--min-coverage", type=int,   default=1)
    ap.add_argument("--gap-open",   type=float, default=-4.0)
    ap.add_argument("--gap-extend", type=float, default=-1.0)
    ap.add_argument("--threads",    type=int,   default=1,
                    help="Parallel worker processes (default: 1).")
    ap.add_argument("--trim-q",       type=int, default=0,
                    help="Mott trim quality cutoff (0 = disable). Trimmed "
                         "flanks are masked to N before projection.")
    ap.add_argument("--trim-min-len", type=int, default=50,
                    help="Skip trimming a read if it would leave fewer bases.")
    ap.add_argument("--include", action="append", default=None,
                    help="Restrict to these sample subdir names (repeatable).")
    ap.add_argument("--forward-suffix", default="_forward",
                    help="Stem suffix marking forward-strand reads (default: _forward).")
    ap.add_argument("--reverse-suffix", default="_reverse",
                    help="Stem suffix marking reverse-strand reads (default: _reverse).")
    args = ap.parse_args()
    _include = set(args.include) if args.include else None

    if args.aligner != "biopython":
        print(f"  WARN: aligner={args.aligner} not bundled; using biopython.", flush=True)

    refs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(args.reference), "fasta")}
    assignments = _read_assignments(args.assignments)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tasks: list[tuple[str, str, str, str]] = []
    for sample_dir in sorted(p for p in args.sample_root.iterdir() if p.is_dir()):
        if _include is not None and sample_dir.name not in _include:
            continue
        for ab1 in sorted(sample_dir.glob("*.ab1")):
            stem = ab1.stem
            key = (sample_dir.name, stem)
            ref_id = assignments.get(key)
            if ref_id is None:
                if args.forward_suffix in stem:
                    mate_stem = stem.replace(args.forward_suffix, args.reverse_suffix)
                elif args.reverse_suffix in stem:
                    mate_stem = stem.replace(args.reverse_suffix, args.forward_suffix)
                else:
                    continue
                ref_id = assignments.get((sample_dir.name, mate_stem))
                if ref_id is None:
                    continue
            if ref_id not in refs:
                continue
            tasks.append((sample_dir.name, str(ab1), ref_id, refs[ref_id]))

    threads = max(1, int(args.threads))
    out_dir_str = str(args.out_dir)
    gap_open, gap_extend = args.gap_open, args.gap_extend

    trim_q, trim_min_len = int(args.trim_q), int(args.trim_min_len)
    if trim_q > 0:
        print(f"[consensus_assembly] end-trimming reads (Mott Q={trim_q}, min_len={trim_min_len})",
              flush=True)

    if threads <= 1 or len(tasks) <= 1:
        print(f"[consensus_assembly] dispatching {len(tasks)} reads sequentially", flush=True)
        for sample, ab1_str, ref_id, ref_seq in tasks:
            print(_process_one(sample, ab1_str, ref_id, ref_seq,
                               out_dir_str, gap_open, gap_extend,
                               trim_q, trim_min_len), flush=True)
    else:
        print(f"[consensus_assembly] dispatching {len(tasks)} reads across {threads} workers", flush=True)
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futs = [pool.submit(_process_one, sample, ab1_str, ref_id, ref_seq,
                                out_dir_str, gap_open, gap_extend,
                                trim_q, trim_min_len)
                    for sample, ab1_str, ref_id, ref_seq in tasks]
            for fut in as_completed(futs):
                print(fut.result(), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
