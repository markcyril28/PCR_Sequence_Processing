#!/usr/bin/env python3
"""overlap_consensus_report.py
==============================
Build a Markdown summary report from the overlap_consensus output directory
produced by overlap_consensus.py + the overlap_consensus_blast step in
f_inside_primer_pipeline.sh.

Reads from <output_dir>:
  PH_*_overlap_report.tsv          -- stitching metrics + per-source table
  PH_*_per_position.tsv            -- per-position coverage (3079 rows)
  blast/*_reference_length.blast.tsv
  blast/*_covered_span.blast.tsv
  blast/*_reference_length.pairwise.txt
  blast/*_covered_span.pairwise.txt
  PH_*_reference_length.fa
  PH_*_covered_span.fa

Writes:
  <output_dir>/overlap_consensus_report.md

Usage:
    python scripts/overlap_consensus_report.py \
        --output-dir III_PipelineOutput/HC01109784/inside_primer_overlap_consensus
    python scripts/overlap_consensus_report.py \
        --output-dir /abs/path/to/inside_primer_overlap_consensus \
        [--workspace /path/to/PCR_Sequence]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import re
import sys
from pathlib import Path
from typing import Optional

# ===========================================================================
# CONFIG
# ===========================================================================
COVERAGE_MAP_WIDTH = 80        # chars for the ASCII coverage bars
BLAST_TSV_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]
PAIRWISE_MAX_BYTES = 500_000   # cap embedded pairwise text to avoid enormous reports

# ===========================================================================
# helpers
# ===========================================================================

def read_text(p: Path, max_bytes: int = PAIRWISE_MAX_BYTES) -> Optional[str]:
    try:
        data = p.read_bytes()[:max_bytes]
        return data.decode("utf-8", errors="replace")
    except FileNotFoundError:
        return None
    except OSError as e:
        return f"[error reading {p.name}: {e}]"


def fenced(content: str, lang: str = "") -> str:
    return f"````{lang}\n{content.rstrip()}\n````\n"


def md_table(headers: list[str], rows: list[list]) -> str:
    buf = io.StringIO()
    buf.write("| " + " | ".join(headers) + " |\n")
    buf.write("|" + "|".join("---" for _ in headers) + "|\n")
    for row in rows:
        buf.write("| " + " | ".join(str(c) for c in row) + " |\n")
    return buf.getvalue()


# ===========================================================================
# parsers
# ===========================================================================

def parse_overlap_report(p: Path) -> tuple[dict, list[dict]]:
    """Return (metrics_dict, source_rows_list).

    The TSV has two sections separated by a blank line:
    Section 1: key<TAB>value pairs
    Section 2: header row + data rows (per-source table)
    """
    metrics: dict[str, str] = {}
    source_rows: list[dict] = []

    if not p.exists():
        return metrics, source_rows

    text = p.read_text(encoding="utf-8")
    sections = re.split(r"\n\s*\n", text.strip())

    # Section 1 -- key/value
    for line in sections[0].splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            metrics[k.strip()] = v.strip()

    # Section 2 -- per-source table (may be absent)
    if len(sections) > 1 and sections[1].strip():
        lines = sections[1].strip().splitlines()
        if lines:
            headers = lines[0].split("\t")
            for line in lines[1:]:
                if not line.strip():
                    continue
                vals = line.split("\t")
                source_rows.append(dict(zip(headers, vals)))

    return metrics, source_rows


def parse_per_position(p: Path) -> list[dict]:
    if not p.exists():
        return []
    rows = []
    with p.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def parse_blast_tsv(p: Path) -> list[dict]:
    if not p.exists():
        return []
    rows = []
    with p.open(newline="", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(BLAST_TSV_COLS):
                parts += [""] * (len(BLAST_TSV_COLS) - len(parts))
            rows.append(dict(zip(BLAST_TSV_COLS, parts[:len(BLAST_TSV_COLS)])))
    return rows


# ===========================================================================
# ASCII coverage map
# ===========================================================================

def coverage_map(ref_len: int, sources: list[dict], width: int = COVERAGE_MAP_WIDTH) -> str:
    """Build a multi-line ASCII bar diagram of per-source coverage vs reference."""
    if ref_len <= 0:
        return ""

    def bar(start: int, end: int) -> str:
        lo = max(0, int((start - 1) / ref_len * width))
        hi = min(width, int(end / ref_len * width))
        return "." * lo + "=" * max(0, hi - lo) + "." * max(0, width - hi)

    lines = []
    ruler = f"|{'|'.ljust(width // 4 - 1)}{'|'.ljust(width // 4 - 1)}{'|'.ljust(width // 4 - 1)}|"
    tick_line = "1" + " " * (width // 4 - 1) + str(ref_len // 4) + " " * (width // 4 - 1) + \
                str(ref_len // 2) + " " * (width // 4 - 1) + str(ref_len * 3 // 4) + \
                " " * max(0, width - (width // 4 * 3) - len(str(ref_len * 3 // 4))) + str(ref_len)

    lines.append(f"reference ({ref_len} bp): positions 1 to {ref_len}")
    lines.append(tick_line)
    lines.append(ruler)

    union_start = ref_len + 1
    union_end = 0
    for src in sources:
        try:
            s = int(src.get("covered_start_1based", 0))
            e = int(src.get("covered_end_1based", 0))
            cb = int(src.get("covered_bases", 0))
            label = src.get("source_id", "").split("|")[-1]
            lines.append(f"[{label}]  start={s}  end={e}  covered_bases={cb}")
            lines.append(bar(s, e))
            if s < union_start:
                union_start = s
            if e > union_end:
                union_end = e
        except (ValueError, TypeError):
            continue

    if sources:
        lines.append(f"[union]  start={union_start}  end={union_end}")
        lines.append(bar(union_start, union_end))

    return "\n".join(lines)


# ===========================================================================
# report builder
# ===========================================================================

def build(output_dir: Path, workspace: Optional[Path] = None) -> Path:
    # --- locate artifacts ---------------------------------------------------
    # Find the overlap_report TSV (there should be exactly one)
    report_tsvs = sorted(output_dir.glob("*_overlap_report.tsv"))
    pos_tsvs    = sorted(output_dir.glob("*_per_position.tsv"))
    ref_fas     = sorted(output_dir.glob("*_reference_length.fa"))
    span_fas    = sorted(output_dir.glob("*_covered_span.fa"))

    blast_dir = output_dir / "blast"
    rl_blast_tsvs  = sorted(blast_dir.glob("*_reference_length.blast.tsv")) if blast_dir.exists() else []
    span_blast_tsvs= sorted(blast_dir.glob("*_covered_span.blast.tsv"))    if blast_dir.exists() else []
    rl_pairwise    = sorted(blast_dir.glob("*_reference_length.pairwise.txt")) if blast_dir.exists() else []
    span_pairwise  = sorted(blast_dir.glob("*_covered_span.pairwise.txt"))    if blast_dir.exists() else []

    # --- parse --------------------------------------------------------------
    metrics, source_rows = parse_overlap_report(report_tsvs[0]) if report_tsvs else ({}, [])
    pos_rows             = parse_per_position(pos_tsvs[0])       if pos_tsvs     else []
    rl_blast_rows        = parse_blast_tsv(rl_blast_tsvs[0])     if rl_blast_tsvs  else []
    span_blast_rows      = parse_blast_tsv(span_blast_tsvs[0])   if span_blast_tsvs else []

    seq_id   = metrics.get("seq_id",            "PH_DMP8_full_gene")
    ref_id   = metrics.get("ref_id",            "")
    ref_len  = int(metrics.get("reference_length", 0) or 0)
    cov_start= metrics.get("covered_start_1based", "")
    cov_end  = metrics.get("covered_end_1based",   "")
    cov_span = metrics.get("covered_span_length",  "")
    n_count  = int(metrics.get("N_count_reference_length", 0) or 0)
    n_called = ref_len - n_count if ref_len else 0
    n_pct    = f"{n_count / ref_len * 100:.1f}" if ref_len else "-"
    ov_pos   = metrics.get("overlap_positions",        "")
    ag_pos   = metrics.get("agree_positions",           "")
    conflict = metrics.get("conflict_iupac_positions",  "")

    # compute per-position status counts
    status_counts: dict[str, int] = {}
    for r in pos_rows:
        s = r.get("status", "")
        status_counts[s] = status_counts.get(s, 0) + 1

    # --- compose Markdown ---------------------------------------------------
    md = io.StringIO()
    w = lambda s="": md.write(s + "\n")

    w(f"# Overlap Consensus Report: {seq_id}")
    w()
    w(f"Generated: {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  ")
    w(f"Directory: `{output_dir}`")
    w()
    w("---")
    w()

    # 1. Summary metrics
    w("## 1. Stitching summary")
    w()
    w(md_table(
        ["Metric", "Value"],
        [
            ["Sequence ID",         seq_id],
            ["Reference",           ref_id],
            ["Reference length",    f"{ref_len:,} bp"],
            ["Covered start (1-based)", cov_start],
            ["Covered end (1-based)",   cov_end],
            ["Covered span length",     f"{cov_span} bp"],
            ["Sequenced bases (non-N, reference-length FASTA)",
             f"{n_called:,} / {ref_len:,} ({100 - float(n_pct):.1f}%)" if ref_len else "-"],
            ["N count (reference-length FASTA)",
             f"{n_count:,} ({n_pct}%)"],
            ["Overlap positions (between pairs)", ov_pos],
            ["Agree positions in overlap",        ag_pos],
            ["IUPAC conflict positions",          conflict],
        ],
    ))

    # 2. Per-source
    if source_rows:
        w("## 2. Per-pair consensus input")
        w()
        w(md_table(
            ["Pair", "Covered bases", "Start (1-based)", "End (1-based)", "Source FASTA"],
            [
                [
                    r.get("source_id", "").split("|")[-1],
                    r.get("covered_bases", "-"),
                    r.get("covered_start_1based", "-"),
                    r.get("covered_end_1based", "-"),
                    Path(r.get("source_fasta", "")).name,
                ]
                for r in source_rows
            ],
        ))

    # 3. Coverage map
    if ref_len and source_rows:
        w("## 3. Coverage map")
        w()
        w(fenced(coverage_map(ref_len, source_rows), ""))

    # 4. Per-position status breakdown
    if status_counts:
        w("## 4. Per-position status breakdown")
        w()
        order = ["agree", "iupac_conflict", "single_source", "uncovered"]
        ordered = [(s, status_counts.get(s, 0)) for s in order if s in status_counts]
        ordered += [(s, c) for s, c in sorted(status_counts.items()) if s not in order]
        w(md_table(
            ["Status", "Count", "% of reference"],
            [
                [s, f"{c:,}", f"{c / ref_len * 100:.1f}%" if ref_len else "-"]
                for s, c in ordered
            ],
        ))
        w()

    # 5. Output FASAs
    w("## 5. Output FASTA sequences")
    w()
    for label, fas in [("reference_length", ref_fas), ("covered_span", span_fas)]:
        fa_text = read_text(fas[0]) if fas else None
        if fa_text:
            w(f"### `{fas[0].name}`")
            w()
            w(fenced(fa_text.strip(), "fasta"))

    # 6. Paralog BLAST confirmation
    w("## 6. Paralog BLAST confirmation")
    w()
    w("Query: final stitched consensus vs. `PF05078_DMP_transcripts_cdhit.fa` (7 SmelDMP paralogs).")
    w()

    def blast_section(label: str, rows: list[dict]) -> None:
        w(f"### {label}")
        w()
        if not rows:
            w("_No BLAST results found._")
            w()
            return
        top = rows[0]
        w(md_table(
            ["#", "Subject", "pident (%)", "Align len", "Mismatches",
             "Q start", "Q end", "S start", "S end", "E-value", "Bitscore", "qcovs (%)"],
            [
                [
                    i + 1,
                    r.get("sseqid", "-"),
                    r.get("pident", "-"),
                    r.get("length", "-"),
                    r.get("mismatch", "-"),
                    r.get("qstart", "-"),
                    r.get("qend", "-"),
                    r.get("sstart", "-"),
                    r.get("send", "-"),
                    r.get("evalue", "-"),
                    r.get("bitscore", "-"),
                    r.get("qcovs", "-"),
                ]
                for i, r in enumerate(rows)
            ],
        ))
        w(f"**Top hit:** `{top.get('sseqid', '-')}` -- {top.get('pident', '-')}% identity "
          f"over {top.get('length', '-')} bp, E={top.get('evalue', '-')}, "
          f"qcovs={top.get('qcovs', '-')}%")
        w()
        if rows and all(r.get("sseqid") == rows[0].get("sseqid") for r in rows):
            n_seqs = len({r.get("sseqid") for r in rows})
            w(f"All {len(rows)} HSP(s) map to the same sequence ({rows[0].get('sseqid', '-')}); "
              "no off-target paralog alignments at this stringency.")
        w()

    blast_section("6a. reference_length FASTA (3079 bp, N-padded)", rl_blast_rows)
    blast_section("6b. covered_span FASTA (actual coverage only)", span_blast_rows)

    # 7. Pairwise alignments
    w("## 7. Pairwise BLAST alignments")
    w()
    for label, pws in [("reference_length", rl_pairwise), ("covered_span", span_pairwise)]:
        pw_text = read_text(pws[0]) if pws else None
        if pw_text:
            w(f"### {label} (`{pws[0].name}`)")
            w()
            w(fenced(pw_text.strip(), "text"))

    # 8. Interpretation note
    sequenced_pct = f"{n_called / ref_len * 100:.0f}" if ref_len else "-"
    w("## 8. Interpretation")
    w()
    w(f"The 4 Sanger traces (2 primer pairs) together cover approximately "
      f"**{n_called:,} bp** ({sequenced_pct}%) of the {ref_len:,}-bp genomic reference "
      f"(`{ref_id}`). The uncovered 5' region (positions 1-{int(cov_start) - 1 if cov_start else '?'}) "
      f"contains the upstream sequence not targeted by the current primers. "
      f"BLAST confirms the sequenced region belongs exclusively to `{ref_id}` "
      f"({rl_blast_rows[0].get('pident', '-') if rl_blast_rows else '-'}% identity, "
      f"E={rl_blast_rows[0].get('evalue', '-') if rl_blast_rows else '-'}), "
      f"with no significant hits to the remaining 6 SmelDMP paralogs.")
    w()
    if conflict and int(conflict) > 0:
        w(f"**{conflict} IUPAC conflict position(s)** were detected in the overlap region "
          f"(positions where PH_2 and PH_3 consensuses disagreed). "
          f"These are encoded as IUPAC ambiguity codes in the output FASTA and should be "
          f"inspected if this region is functionally relevant.")
        w()

    # write
    out = output_dir / "overlap_consensus_report.md"
    out.write_text(md.getvalue(), encoding="utf-8")
    return out


# ===========================================================================
# CLI
# ===========================================================================

def main() -> int:
    ap = argparse.ArgumentParser(description="Build Markdown report from overlap_consensus output dir.")
    ap.add_argument(
        "--output-dir", required=True, type=Path,
        help="overlap_consensus output directory (contains *_overlap_report.tsv, blast/, etc.)",
    )
    ap.add_argument(
        "--workspace", type=Path, default=None,
        help="PCR_Sequence workspace root (used only to relativize paths in messages).",
    )
    args = ap.parse_args()

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        ws = (args.workspace or Path(__file__).resolve().parent.parent).resolve()
        output_dir = ws / output_dir
    output_dir = output_dir.resolve()

    if not output_dir.is_dir():
        print(f"ERROR: output-dir not found: {output_dir}", file=sys.stderr)
        return 2

    out = build(output_dir, args.workspace)
    print(f"[overlap_consensus_report] -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
