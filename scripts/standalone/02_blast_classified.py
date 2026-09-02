#!/usr/bin/env python3
"""
02_blast_classified.py: BLAST every .ab1 read inside a tier-classified folder
against the reference DMP transcripts and emit pairwise alignments.

Walks <classified_root>/<tier>/<sample>/*.ab1, converts each trace to FASTA,
runs blastn (pairwise output) against the reference, and writes:
  <out_root>/<tier>/<sample>/<read>.fa
  <out_root>/<tier>/<sample>/<read>.blast.txt        # pairwise alignment(s)
  <out_root>/<tier>/<sample>/<read>.blast.tsv        # tabular best hits
  <out_root>/blast_summary.tsv                        # one row per read (best hit)
"""
from __future__ import annotations

# =============================================================================
# Edit these defaults, or override on the CLI.
# =============================================================================
DEFAULT_CLASSIFIED_ROOT = "III_PipelineOutput/low_quality_saver_pipeline/01_QualityClassified"
DEFAULT_REFERENCE       = "II_QualityClassified/PF05078_DMP_transcripts_cdhit.fa"
DEFAULT_OUT_ROOT        = "III_PipelineOutput/low_quality_saver_pipeline/02_BlastAssignment/classified_alignments"
DEFAULT_TIERS           = ["2_MediumQuality", "3_LowQuality"]
DEFAULT_EVALUE          = "1e-10"
DEFAULT_WORD_SIZE       = 11
DEFAULT_MAX_TARGETS     = 5
DEFAULT_THREADS         = 4
# =============================================================================

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


TAB_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
TAB_HEADER = TAB_FMT.split(" ", 1)[1].split()


def ab1_to_fasta(ab1: Path, out_fa: Path, seq_id: str) -> None:
    rec = next(SeqIO.parse(str(ab1), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in ABI file: {ab1}")
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    s = str(rec.seq)
    with out_fa.open("w") as fh:
        fh.write(f">{seq_id}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")


def ensure_db(reference: Path, db_dir: Path) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "ref"
    if not (db_dir / "ref.nhr").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(reference), "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True, stdout=subprocess.DEVNULL,
        )
    return db_prefix


def run_blast(query_fa: Path, db_prefix: Path, out_pairwise: Path, out_tab: Path,
              evalue: str, word_size: int, max_targets: int, threads: int) -> None:
    # Pairwise alignment view
    with out_pairwise.open("w") as fh:
        subprocess.run(
            ["blastn", "-query", str(query_fa), "-db", str(db_prefix),
             "-evalue", evalue, "-word_size", str(word_size),
             "-max_target_seqs", str(max_targets),
             "-num_threads", str(threads), "-outfmt", "0"],
            check=True, stdout=fh,
        )
    # Tabular hits
    with out_tab.open("w") as fh:
        fh.write("\t".join(TAB_HEADER) + "\n")
        subprocess.run(
            ["blastn", "-query", str(query_fa), "-db", str(db_prefix),
             "-evalue", evalue, "-word_size", str(word_size),
             "-max_target_seqs", str(max_targets),
             "-num_threads", str(threads), "-outfmt", TAB_FMT],
            check=True, stdout=fh,
        )


def best_hit_row(tab_file: Path) -> list[str] | None:
    lines = [l for l in tab_file.read_text().splitlines()[1:] if l.strip()]
    return lines[0].split("\t") if lines else None


def main() -> int:
    here = Path(__file__).resolve().parents[2]

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--classified-root", type=Path, default=here / DEFAULT_CLASSIFIED_ROOT)
    ap.add_argument("--reference",       type=Path, default=here / DEFAULT_REFERENCE)
    ap.add_argument("--out-root",        type=Path, default=here / DEFAULT_OUT_ROOT)
    ap.add_argument("--tiers",           nargs="+", default=DEFAULT_TIERS)
    ap.add_argument("--evalue",          default=DEFAULT_EVALUE)
    ap.add_argument("--word-size",       type=int, default=DEFAULT_WORD_SIZE)
    ap.add_argument("--max-targets",     type=int, default=DEFAULT_MAX_TARGETS)
    ap.add_argument("--threads",         type=int, default=DEFAULT_THREADS)
    args = ap.parse_args()

    if shutil.which("blastn") is None or shutil.which("makeblastdb") is None:
        raise SystemExit("ERROR: blastn / makeblastdb not on PATH (activate env or run setup.sh)")
    if not args.reference.is_file():
        raise SystemExit(f"ERROR: reference not found: {args.reference}")
    if not args.classified_root.is_dir():
        raise SystemExit(f"ERROR: classified root not found: {args.classified_root}")

    args.out_root.mkdir(parents=True, exist_ok=True)
    db_prefix = ensure_db(args.reference, args.out_root / "db")

    summary_path = args.out_root / "blast_summary.tsv"
    summary_cols = ["tier", "sample", "read", *TAB_HEADER]
    with summary_path.open("w") as sf:
        sf.write("\t".join(summary_cols) + "\n")

        n_reads = n_hits = 0
        for tier in args.tiers:
            tier_dir = args.classified_root / tier
            if not tier_dir.is_dir():
                print(f"[skip] tier folder missing: {tier_dir}", flush=True)
                continue
            for sample_dir in sorted(p for p in tier_dir.iterdir() if p.is_dir()):
                for ab1 in sorted(sample_dir.glob("*.ab1")):
                    n_reads += 1
                    stem = ab1.stem
                    out_dir = args.out_root / tier / sample_dir.name
                    out_dir.mkdir(parents=True, exist_ok=True)
                    fa  = out_dir / f"{stem}.fa"
                    pw  = out_dir / f"{stem}.blast.txt"
                    tab = out_dir / f"{stem}.blast.tsv"

                    ab1_to_fasta(ab1, fa, stem)
                    run_blast(fa, db_prefix, pw, tab,
                              args.evalue, args.word_size, args.max_targets, args.threads)

                    row = best_hit_row(tab)
                    if row:
                        n_hits += 1
                        sf.write("\t".join([tier, sample_dir.name, stem, *row]) + "\n")
                    else:
                        sf.write("\t".join([tier, sample_dir.name, stem,
                                            stem, "NO_HIT", *(["-"] * (len(TAB_HEADER) - 2))]) + "\n")
                    print(f"[blast] {tier}/{sample_dir.name}/{stem}  "
                          f"-> {(row[1] if row else 'NO_HIT')}", flush=True)

    print(f"\n[done] {n_hits}/{n_reads} reads had a BLAST hit")
    print(f"[done] summary: {summary_path}")
    print(f"[done] alignments under: {args.out_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
