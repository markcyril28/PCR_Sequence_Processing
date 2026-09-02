#!/usr/bin/env python3
"""Step 5: write a run manifest + roll-up CSV.

Captures: which method(s) were selected, which steps actually produced output,
how many reads per quality tier, BLAST hit counts, and per-sample N% in the
final consensus FASTAs.

Accepts one or more methods via --method (space- or comma-separated). The
manifest stores all of them and lists every method-specific call directory.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import shutil
from pathlib import Path


def _count_reads_per_tier(classify_dir: Path) -> dict[str, int]:
    out = {}
    if not classify_dir.exists(): return out
    for tier in sorted(classify_dir.iterdir()):
        if tier.is_dir():
            out[tier.name] = sum(1 for _ in tier.rglob("*.ab1"))
    return out


def _blast_summary(tsv: Path) -> dict:
    """Return counts of assigned / no-hit / fail rows in assignments.tsv."""
    if not tsv.exists():
        return {"total": 0, "assigned": 0, "no_hit": 0, "fail_threshold": 0}
    assigned = no_hit = fail = total = 0
    with tsv.open() as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            total += 1
            ref = parts[2]
            if ref == "NO_HIT":
                no_hit += 1
            elif ref == "FAIL_THRESHOLD":
                fail += 1
            elif ref not in ("-", ""):
                assigned += 1
    return {"total": total, "assigned": assigned,
            "no_hit": no_hit, "fail_threshold": fail}


def _rel(p: Path, root: Path) -> str:
    try:
        return str(p.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(p)


def _split_methods(arg: str) -> list[str]:
    return [m for m in arg.replace(",", " ").split() if m]


def _split_call_dirs(arg: str | None) -> list[Path]:
    if not arg:
        return []
    return [Path(p) for p in arg.replace(",", " ").split() if p]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-root",   required=True, type=Path)
    ap.add_argument("--classify-dir",  required=True, type=Path)
    ap.add_argument("--blast-dir",     required=True, type=Path)
    ap.add_argument("--call-dir",      required=False, default="",
                    help="Space- or comma-separated list of method call dirs.")
    ap.add_argument("--consensus-dir", required=False, default="",
                    help="Space- or comma-separated list of merged-consensus dirs (one per method).")
    ap.add_argument("--report-dir",    required=True, type=Path)
    ap.add_argument("--method",        required=True,
                    help="One or more methods (space- or comma-separated).")
    ap.add_argument("--config",        required=True, type=Path)
    args = ap.parse_args()

    args.report_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.config, args.report_dir / "config_used.toml")

    methods   = _split_methods(args.method)
    call_dirs = _split_call_dirs(args.call_dir)
    consensus_dirs = _split_call_dirs(args.consensus_dir)
    output_root = args.output_root.resolve()

    manifest = {
        "timestamp":         dt.datetime.now().isoformat(timespec="seconds"),
        "output_root":       _rel(args.output_root, output_root) or ".",
        "methods":           methods,
        "tiers":             _count_reads_per_tier(args.classify_dir),
        "blast_assignments": _rel(args.blast_dir / "assignments.tsv", output_root),
        "blast_summary":     _blast_summary(args.blast_dir / "assignments.tsv"),
        "call_dirs":         [_rel(d, output_root) for d in call_dirs],
        "consensus_dirs":    [_rel(d, output_root) for d in consensus_dirs],
    }

    # Per-consensus N% summary. 07_step_merge_consensus.py emits files named
    # "<base>_consensus.fa"; older runs may have used "<base>.consensus.fa".
    # Each consensus dir corresponds to one method, so tag every row with the
    # source dir's basename so downstream consumers can tell methods apart.
    rows = []
    seen: set[Path] = set()
    for cdir in consensus_dirs:
        if not cdir.exists():
            continue
        source = cdir.name
        for pattern in ("*_consensus.fa", "*.consensus.fa"):
            for fa in sorted(cdir.rglob(pattern)):
                if fa in seen:
                    continue
                seen.add(fa)
                text = fa.read_text().splitlines()
                seq = "".join(l for l in text if not l.startswith(">")).upper()
                if not seq:
                    continue
                n = seq.count("N")
                amplicon = fa.stem
                for suffix in ("_consensus", ".consensus"):
                    if amplicon.endswith(suffix):
                        amplicon = amplicon[: -len(suffix)]
                        break
                rows.append({
                    "source":   source,
                    "sample":   fa.parent.name,
                    "amplicon": amplicon,
                    "length":   len(seq),
                    "N_count":  n,
                    "N_pct":    f"{100*n/len(seq):.2f}",
                    "path":     _rel(fa, output_root),
                })

    csv_path = args.report_dir / "consensus_summary.csv"
    fieldnames = ["source", "sample", "amplicon", "length", "N_count", "N_pct", "path"]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    manifest["consensus_summary_csv"] = _rel(csv_path, output_root)
    manifest["consensus_records"] = len(rows)

    (args.report_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"  -> {args.report_dir/'manifest.json'}")
    print(f"  -> {csv_path} ({len(rows)} consensus records)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
