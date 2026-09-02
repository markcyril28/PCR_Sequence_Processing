#!/usr/bin/env python3
"""Merge Tracy per-strand FASTAs (forward + reverse-complement-of-reverse) into
one consensus per amplicon. Designed to be called by pipeline.sh.

Inputs (positional):
    --forward-fa  forward strand corrected FASTA from `tracy align`
    --reverse-fa  reverse strand corrected FASTA from `tracy align`
                  (will be reverse-complemented before merging)
    --forward-q   matching per-base Phred qualities (whitespace-separated ints)
    --reverse-q   matching per-base Phred qualities for the reverse read
    --reference   reference transcript FASTA (single record)
    --out-fa      output consensus FASTA path
    --min-q       per-strand quality floor; below -> N on that strand
    --policy      iupac | n | higher_q   (handling of forward/reverse mismatch)
    --min-cov     1 = accept single-strand bases, 2 = require both strands

Strategy:
    Both strand reads are pairwise-aligned to the reference (Biopython
    PairwiseAligner, global, simple match scoring). At every reference column
    the better-quality strand call (above min_q, after reverse-complementing
    the reverse strand) is taken; conflicts are resolved by the chosen policy.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}


def _read_single(fa: Path) -> str:
    rec = next(SeqIO.parse(str(fa), "fasta"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in FASTA file: {fa}")
    return str(rec.seq).upper()


def _read_quals(qpath: Path, expected_len: int) -> list[int]:
    if not qpath.exists():
        return [0] * expected_len
    vals = qpath.read_text().split()
    qs = [int(v) for v in vals]
    if len(qs) != expected_len:
        # Pad/truncate defensively.
        if len(qs) < expected_len:
            qs += [0] * (expected_len - len(qs))
        else:
            qs = qs[:expected_len]
    return qs


def _align_to_ref(ref: str, read: str, quals: list[int]) -> tuple[list[str], list[int], float]:
    """Return per-reference-position (base, quality, score). Insertions in the
    read are dropped; deletions are recorded as '-' with quality 0."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    aln = aligner.align(ref, read)[0]
    aligned_ref, aligned_read = str(aln[0]), str(aln[1])

    out_bases: list[str] = []
    out_quals: list[int] = []
    read_idx = 0
    for r_base, q_base in zip(aligned_ref, aligned_read):
        if r_base == "-":
            # Insertion in read relative to ref: skip (consume read base)
            if q_base != "-":
                read_idx += 1
            continue
        if q_base == "-":
            out_bases.append("-")
            out_quals.append(0)
        else:
            out_bases.append(q_base)
            out_quals.append(quals[read_idx] if read_idx < len(quals) else 0)
            read_idx += 1
    return out_bases, out_quals, float(aln.score)


def _orient_to_ref(ref: str, seq: str, quals: list[int]
                   ) -> tuple[str, list[int], str, float, float]:
    """Pick the orientation of a strand (as-is vs reverse-complement) that
    aligns better to ``ref``. Returns (seq, quals, orientation, score_as_is,
    score_revcomp). Robust to mislabeled F/R inputs and to upstream tools
    that may already orient their output.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2; aligner.mismatch_score = -1
    aligner.open_gap_score = -4; aligner.extend_gap_score = -1
    score_asis = float(aligner.align(ref, seq)[0].score)
    rc_seq = str(Seq(seq).reverse_complement())
    score_rc = float(aligner.align(ref, rc_seq)[0].score)
    if score_rc > score_asis:
        return rc_seq, list(reversed(quals)), "revcomp", score_asis, score_rc
    return seq, list(quals), "as_is", score_asis, score_rc


def _merge(ref: str,
           f_bases: list[str], f_quals: list[int],
           r_bases: list[str], r_quals: list[int],
           min_q: int, policy: str, min_cov: int) -> str:
    consensus = []
    for i, ref_base in enumerate(ref):
        fb = f_bases[i] if i < len(f_bases) else "-"
        rb = r_bases[i] if i < len(r_bases) else "-"
        fq = f_quals[i] if i < len(f_quals) else 0
        rq = r_quals[i] if i < len(r_quals) else 0

        f_ok = fb not in ("-", "N") and fq >= min_q
        r_ok = rb not in ("-", "N") and rq >= min_q

        if not f_ok and not r_ok:
            consensus.append("N" if min_cov >= 1 else ref_base)
            continue
        if f_ok and not r_ok:
            consensus.append(fb if min_cov == 1 else "N")
            continue
        if r_ok and not f_ok:
            consensus.append(rb if min_cov == 1 else "N")
            continue
        # Both strands have a confident call.
        if fb == rb:
            consensus.append(fb)
        elif policy == "higher_q":
            consensus.append(fb if fq >= rq else rb)
        elif policy == "n":
            consensus.append("N")
        else:  # iupac
            consensus.append(IUPAC.get(frozenset(fb + rb), "N"))
    return "".join(consensus)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--forward-fa", required=True, type=Path)
    p.add_argument("--reverse-fa", required=True, type=Path)
    p.add_argument("--forward-q", required=True, type=Path)
    p.add_argument("--reverse-q", required=True, type=Path)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--ref-id", required=True)
    p.add_argument("--out-fa", required=True, type=Path)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--min-q", type=int, default=20)
    p.add_argument("--policy", default="iupac", choices=["iupac", "n", "higher_q"])
    p.add_argument("--min-cov", type=int, default=1, choices=[1, 2])
    args = p.parse_args()

    # Pull the chosen reference record.
    ref_seq = None
    for rec in SeqIO.parse(str(args.reference), "fasta"):
        if rec.id == args.ref_id:
            ref_seq = str(rec.seq).upper()
            break
    if ref_seq is None:
        raise SystemExit(f"Reference id {args.ref_id} not found in {args.reference}")

    f_seq_raw = _read_single(args.forward_fa)
    r_seq_raw = _read_single(args.reverse_fa)

    f_quals_raw = _read_quals(args.forward_q, len(f_seq_raw))
    r_quals_raw = _read_quals(args.reverse_q, len(r_seq_raw))

    # Auto-orient each strand to the reference (handles either Tracy's raw
    # output or pre-oriented inputs, and is robust to mislabeled files).
    f_seq, f_quals, f_orient, f_score_asis, f_score_rc = _orient_to_ref(
        ref_seq, f_seq_raw, f_quals_raw)
    r_seq, r_quals, r_orient, r_score_asis, r_score_rc = _orient_to_ref(
        ref_seq, r_seq_raw, r_quals_raw)

    f_bases, f_quals_aln, _ = _align_to_ref(ref_seq, f_seq, f_quals)
    r_bases, r_quals_aln, _ = _align_to_ref(ref_seq, r_seq, r_quals)

    cons = _merge(ref_seq, f_bases, f_quals_aln, r_bases, r_quals_aln,
                  args.min_q, args.policy, args.min_cov)

    args.out_fa.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fa.open("w") as fh:
        fh.write(f">{args.sample_id} ref={args.ref_id} len={len(cons)} "
                 f"min_q={args.min_q} policy={args.policy} min_cov={args.min_cov} "
                 f"forward_orientation={f_orient} reverse_orientation={r_orient}\n")
        for i in range(0, len(cons), 80):
            fh.write(cons[i:i + 80] + "\n")

    n_count = cons.count("N")
    print(f"[consensus] {args.sample_id} -> {args.out_fa.name}  "
          f"len={len(cons)}  N={n_count} ({100*n_count/len(cons):.1f}%)  "
          f"f_orient={f_orient} r_orient={r_orient}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
