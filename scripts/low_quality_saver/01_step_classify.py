#!/usr/bin/env python3
"""Step 1: Phred-tier classification of every .ab1 read.

Reads are paired by gene (forward + reverse of the same amplicon share a
base stem). The pair is classified by the LOWER mean Phred score of its two
reads, so both forward and reverse always land in the same tier folder.

For each *.phd.1 in <sample-root>/<sample>/, compute mean Phred, group by
gene, and copy each read's companion files (.ab1 .phd.1 .txt) into:
    <out-dir>/1_HighQuality/<sample>/
    <out-dir>/2_MediumQuality/<sample>/
    <out-dir>/3_LowQuality/<sample>/
Writes <out-dir>/quality_report.csv.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import shutil
from pathlib import Path


HIGH_THRESHOLD  = 40.0
MEDIUM_THRESHOLD = 20.0


def parse_phd_qualities(phd_path: Path) -> list[int]:
    qs: list[int] = []
    in_dna = False
    with phd_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN_DNA":
                in_dna = True
                continue
            if line == "END_DNA":
                break
            if not in_dna or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    qs.append(int(parts[1]))
                except ValueError:
                    pass
    return qs


def gene_key(stem: str, forward_suffix: str, reverse_suffix: str) -> tuple[str, str | None]:
    """Split a read stem into (gene_id, direction).

    Forward/reverse direction is identified by the trailing strand suffix
    (default ``_forward`` / ``_reverse``; override per-batch when the lab uses
    e.g. ``F`` / ``R`` instead). The matched suffix is stripped to yield the
    gene-pair key. If neither suffix matches, direction is ``None``.
    """
    if forward_suffix and stem.endswith(forward_suffix):
        return stem[: -len(forward_suffix)], "forward"
    if reverse_suffix and stem.endswith(reverse_suffix):
        return stem[: -len(reverse_suffix)], "reverse"
    return stem, None


def classify(mean_q: float, high: float, medium: float) -> str:
    if mean_q >= high:
        return "1_HighQuality"
    if mean_q >= medium:
        return "2_MediumQuality"
    return "3_LowQuality"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-root", required=True, type=Path)
    ap.add_argument("--out-dir",     required=True, type=Path)
    ap.add_argument("--high",   type=float, default=HIGH_THRESHOLD)
    ap.add_argument("--medium", type=float, default=MEDIUM_THRESHOLD)
    # Per-batch overrides: empty means "use every direct subdir of sample-root".
    ap.add_argument("--include", action="append", default=None,
                    help="Restrict to these sample subdir names (repeatable). "
                         "If omitted, every direct subdir of --sample-root is processed.")
    ap.add_argument("--forward-suffix", default="_forward",
                    help="Stem suffix identifying forward-strand reads (default: _forward).")
    ap.add_argument("--reverse-suffix", default="_reverse",
                    help="Stem suffix identifying reverse-strand reads (default: _reverse).")
    ap.add_argument("--read-glob", action="append", default=None,
                    help="Restrict reads to these filename globs, relative to each sample dir. "
                         "Repeat for forward and reverse patterns.")
    args = ap.parse_args()

    include = set(args.include) if args.include else None

    rows = []
    for sample_dir in sorted(p for p in args.sample_root.iterdir() if p.is_dir()):
        if sample_dir == args.out_dir or sample_dir.name.startswith("III_PipelineOutput"):
            continue
        if include is not None and sample_dir.name not in include:
            continue

        # Collect per-read metrics grouped by gene
        groups: dict[str, list[dict]] = {}
        for phd in sorted(sample_dir.glob("*.phd.1")):
            stem = phd.name[: -len(".phd.1")]
            ab1_name = f"{stem}.ab1"
            if args.read_glob and not any(fnmatch.fnmatch(ab1_name, pat) for pat in args.read_glob):
                continue
            qs = parse_phd_qualities(phd)
            n = len(qs)
            mean_q = sum(qs) / n if n else 0.0
            gene, direction = gene_key(stem, args.forward_suffix, args.reverse_suffix)
            groups.setdefault(gene, []).append(
                {"stem": stem, "direction": direction, "mean_q": mean_q, "n": n}
            )

        for gene, reads in sorted(groups.items()):
            pair_min_q = min(r["mean_q"] for r in reads)
            tier = classify(pair_min_q, args.high, args.medium)
            dest = args.out_dir / tier / sample_dir.name
            dest.mkdir(parents=True, exist_ok=True)

            for r in reads:
                stem = r["stem"]
                for comp in sample_dir.glob(f"{stem}.*"):
                    shutil.copy2(comp, dest / comp.name)
                rows.append({
                    "sample":         sample_dir.name,
                    "gene":           gene,
                    "direction":      r["direction"] or "",
                    "read":           stem,
                    "bases":          r["n"],
                    "mean_phred":     f"{r['mean_q']:.2f}",
                    "pair_min_phred": f"{pair_min_q:.2f}",
                    "tier":           tier,
                })
                print(f"  {sample_dir.name}/{stem}  n={r['n']} "
                      f"meanQ={r['mean_q']:.2f} pairMinQ={pair_min_q:.2f} -> {tier}")

    rep = args.out_dir / "quality_report.csv"
    with rep.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "sample", "gene", "direction", "read",
            "bases", "mean_phred", "pair_min_phred", "tier"])
        w.writeheader(); w.writerows(rows)
    print(f"  -> {rep} ({len(rows)} reads)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

