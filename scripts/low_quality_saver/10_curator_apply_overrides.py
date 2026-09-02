#!/usr/bin/env python3
"""Apply curator single-strand overrides to the 04a trace-review consensus.

When the curator inspects 03a_TraceReview/*ref_align.txt and decides one of
the two reviewed strands is unreliable, they record that decision in
`05_Report/curator_decisions.tsv`. This script rebuilds the affected
04a consensus FASTA from the chosen strand's `.reviewed.fa`, projected onto
reference coordinates so the output length matches the reference (matching the
two-strand consensus produced by `07_step_merge_consensus.py`). This keeps the
downstream `d_tri_consensus_pipeline.sh` (which requires equal-length inputs)
compatible with curator-overridden amplicons.

Inputs (relative to the batch's `low_quality_saver_pipeline/` directory):
    05_Report/curator_decisions.tsv     curator's per-amplicon strand choice
    03a_TraceReview/<sample>/<amp>F.reviewed.fa (or R)
    03a_TraceReview/_ref/<ref>.fa       reference FASTA

Outputs:
    04a_Consensus_TraceReview/<sample>/<amp>_consensus.fa  (reference-length,
        single strand, header carries curator_override=1, strand=F|R)

Idempotent: rerun safely; previous curator-overridden files are overwritten.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

SUFFIX_F = "F"
SUFFIX_R = "R"
LEGACY_SUFFIX_F = "_forward"
LEGACY_SUFFIX_R = "_reverse"

ALIGN_MATCH = 2
ALIGN_MISMATCH = -1
ALIGN_GAP_OPEN = -4
ALIGN_GAP_EXTEND = -1


def _read_fasta(path: Path) -> tuple[str, str]:
    rec = next(SeqIO.parse(str(path), "fasta"), None)
    if rec is None:
        raise SystemExit(f"ERROR: empty FASTA: {path}")
    return rec.id, str(rec.seq).upper()


def _project_to_reference(ref: str, seq: str, min_q: int = 0,
                          quals: list[int] | None = None) -> str:
    """Align ``seq`` to ``ref`` and return a string of len(ref) where each
    position holds the read base aligned to that reference position, or 'N'
    where the read has a gap. Picks the better of as-is vs reverse-complement.
    Quality scores, if provided, mask low-Q calls to 'N'.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = ALIGN_MATCH
    aligner.mismatch_score = ALIGN_MISMATCH
    aligner.open_gap_score = ALIGN_GAP_OPEN
    aligner.extend_gap_score = ALIGN_GAP_EXTEND

    score_asis = float(aligner.align(ref, seq)[0].score)
    rc_seq = str(Seq(seq).reverse_complement())
    score_rc = float(aligner.align(ref, rc_seq)[0].score)

    if score_rc > score_asis:
        seq = rc_seq
        if quals is not None:
            quals = list(reversed(quals))

    aln = aligner.align(ref, seq)[0]
    aref, aread = str(aln[0]), str(aln[1])

    bases: list[str] = []
    ri = 0
    for r_b, c_b in zip(aref, aread):
        if r_b == "-":
            # insertion in read relative to reference: drop it
            if c_b != "-":
                ri += 1
            continue
        if c_b == "-":
            bases.append("N")
        else:
            q = (quals[ri] if quals is not None and ri < len(quals) else min_q)
            if c_b in "ACGTRYSWKMBDHV" and q >= min_q:
                bases.append(c_b)
            else:
                bases.append("N")
            ri += 1
    return "".join(bases)


def _resolve_strand_path(strand_dir: Path, base_stem: str,
                         strand: str) -> Path | None:
    """Try the new (F/R) and legacy (_forward/_reverse) suffix layouts."""
    candidates = [
        strand_dir / f"{base_stem}{strand}.reviewed.fa",
    ]
    if strand == SUFFIX_F:
        candidates.append(strand_dir / f"{base_stem}{LEGACY_SUFFIX_F}.reviewed.fa")
    elif strand == SUFFIX_R:
        candidates.append(strand_dir / f"{base_stem}{LEGACY_SUFFIX_R}.reviewed.fa")
    for c in candidates:
        if c.exists():
            return c
    return None


def _resolve_ref_path(ref_dir: Path, ref_acc: str) -> Path | None:
    p = ref_dir / f"{ref_acc}.fa"
    return p if p.exists() else None


def _update_consensus_summary(summary_path: Path,
                              updates: list[dict]) -> None:
    """Rewrite the curator rows in consensus_summary.csv with fresh
    length/N_count/N_pct values. Non-curator rows are preserved as-is.
    """
    if not summary_path.exists() or not updates:
        return
    by_key = {(u["sample"], u["amplicon"]): u for u in updates}
    rows: list[dict] = []
    fieldnames: list[str] = []
    with summary_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            key = (row.get("sample", ""), row.get("amplicon", ""))
            if key in by_key:
                u = by_key[key]
                row["length"]  = str(u["length"])
                row["N_count"] = str(u["N_count"])
                row["N_pct"]   = f"{u['N_pct']:.2f}"
                strand_label = "F" if u["strand"] == "F" else "R"
                row["note"]    = f"curator: use {strand_label} only"
            rows.append(row)

    # Ensure 'note' column exists.
    if "note" not in fieldnames:
        fieldnames.append("note")

    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def apply_overrides(lqs_root: Path, *, min_q: int = 5) -> int:
    decisions_path = lqs_root / "05_Report" / "curator_decisions.tsv"
    if not decisions_path.exists():
        print(f"[curator] no decisions file at {decisions_path}; nothing to do")
        return 0

    tr_root = lqs_root / "03a_TraceReview"
    ref_dir = tr_root / "_ref"
    out_root = lqs_root / "04a_Consensus_TraceReview"
    summary_path = lqs_root / "05_Report" / "consensus_summary.csv"
    if not ref_dir.is_dir():
        raise SystemExit(f"ERROR: reference directory missing: {ref_dir}")

    updates: list[dict] = []
    n_applied = 0
    with decisions_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            amplicon = row.get("amplicon", "").strip()
            strand   = row.get("chosen_strand", "").strip().upper()
            ref_acc  = row.get("ref", "").strip()
            src_rel  = row.get("source_reviewed_fa", "").strip()
            out_rel  = row.get("new_consensus_path", "").strip()
            if not (amplicon and strand and ref_acc and src_rel and out_rel):
                print(f"[curator] skip incomplete row: {row}")
                continue
            if strand not in {SUFFIX_F, SUFFIX_R}:
                print(f"[curator] skip {amplicon}: unrecognized chosen_strand={strand!r}")
                continue

            src_fa = lqs_root / src_rel
            if not src_fa.exists():
                # Fall back to resolving by sample folder + amplicon stem.
                sample_dir = src_fa.parent
                src_fa = (_resolve_strand_path(sample_dir, amplicon, strand)
                          or src_fa)
            if not src_fa.exists():
                print(f"[curator] ERROR {amplicon}: source reviewed FASTA missing: {src_fa}")
                continue

            ref_fa = _resolve_ref_path(ref_dir, ref_acc)
            if ref_fa is None:
                print(f"[curator] ERROR {amplicon}: reference FASTA missing for {ref_acc}")
                continue

            _, ref_seq = _read_fasta(ref_fa)
            src_hdr, src_seq = _read_fasta(src_fa)
            projected = _project_to_reference(ref_seq, src_seq, min_q=0)

            sample_dir = src_fa.parent.name              # e.g. I_HC01468262
            out_path = lqs_root / out_rel
            out_path.parent.mkdir(parents=True, exist_ok=True)

            cov_note = "forward_only" if strand == SUFFIX_F else "reverse_only"
            header = (
                f">{sample_dir}|{amplicon} ref={ref_acc} "
                f"len={len(projected)} min_q={min_q} policy=single_strand "
                f"strand={strand} curator_override=1 strands={cov_note}"
            )
            with out_path.open("w") as fh:
                fh.write(header + "\n")
                for i in range(0, len(projected), 60):
                    fh.write(projected[i:i + 60] + "\n")
            n_calls = projected.count("N")
            print(f"[curator] {amplicon} strand={strand} ref={ref_acc} "
                  f"len={len(projected)} N={n_calls} -> {out_path.relative_to(lqs_root)}")
            updates.append({
                "sample":    sample_dir,
                "amplicon":  amplicon,
                "strand":    strand,
                "length":    len(projected),
                "N_count":   n_calls,
                "N_pct":     (100.0 * n_calls / max(len(projected), 1)),
            })
            n_applied += 1

    _update_consensus_summary(summary_path, updates)
    print(f"[curator] applied {n_applied} override(s)")
    return n_applied


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--lqs-root", required=True, type=Path,
        help="Batch's low_quality_saver_pipeline directory "
             "(contains 03a_TraceReview/, 04a_Consensus_TraceReview/, 05_Report/)."
    )
    ap.add_argument("--min-q", type=int, default=5,
                    help="Reported in the output FASTA header; "
                         "projection itself is quality-agnostic.")
    args = ap.parse_args()

    if not args.lqs_root.is_dir():
        raise SystemExit(f"ERROR: not a directory: {args.lqs_root}")
    apply_overrides(args.lqs_root, min_q=args.min_q)
    return 0


if __name__ == "__main__":
    sys.exit(main())
