#!/usr/bin/env python3
"""01_build_full_report.py
==========================
One big rerunnable script that aggregates every artifact this PCR_Sequence
workspace produces (chromatograms, quality classification, merge_consensus
pipeline result, trace-review per-strand alignment, attempted merged
trace-review consensus) into a single Markdown report. The report is
PDF-ready (e.g. via `pandoc Report.md -o Report.pdf` or VS Code "Markdown
PDF").

Called by ``c_build_report.sh`` and also available as a standalone CLI.

Usage:
    python scripts/reporting/01_build_full_report.py
    python scripts/reporting/01_build_full_report.py --workspace /path/to/PCR_Sequence
    python scripts/reporting/01_build_full_report.py --no-chromatograms

Idempotent: rerunning overwrites the report and refreshes PNGs.
"""

# ===========================================================================
# === USER-EDITABLE CONFIG (everything important up top, per repo style) ====
# ===========================================================================
from __future__ import annotations
from pathlib import Path
import os

# Default workspace = the parent of this scripts/ folder. Override with --workspace.
DEFAULT_WORKSPACE = Path(__file__).resolve().parents[2]

# --- input locations (relative to workspace) -----------------------------
#
# Pipeline outputs are partitioned per Sanger sequencing batch under:
#   {PIPELINE_OUTPUT_ROOT}/{batch_id}/{LQS_PIPELINE_DIRNAME}/...
#   {PIPELINE_OUTPUT_ROOT}/{batch_id}/{MERGE_PIPELINE_DIRNAME}/...
#   {PIPELINE_OUTPUT_ROOT}/{batch_id}/{REPORT_PIPELINE_DIRNAME}/...
# Use --batch <id> [--sample-dir <I_HCXXX>] to select which batch this report
# covers. Sample dirs may carry any name (e.g. I_HC01463503, I_HC01468262)
# and are resolved below the ``[project].sample_root`` configured for the LQS
# pipeline.
PIPELINE_OUTPUT_ROOT      = "III_PipelineOutput"
LQS_PIPELINE_DIRNAME      = "low_quality_saver_pipeline"
MERGE_PIPELINE_DIRNAME    = "merge_consensus_pipeline"
REPORT_PIPELINE_DIRNAME   = "full_report"

SAMPLE_GLOB              = "I_*"                                # per-sample input dirs (chromatograms live here)
LQS_QUALITY_REPORT_CSV   = "01_QualityClassified/quality_report.csv"  # under LQS_OUTPUT_DIR
LQS_BLAST                = "02_BlastAssignment/assignments.tsv"
LQS_TRACE_REVIEW_DIR     = "03a_TraceReview"
LQS_TRACE_REVIEW_SUMMARY = "03a_TraceReview/trace_review_summary.tsv"
LQS_CONSENSUS_TR_DIR     = "04a_Consensus_TraceReview"
LQS_CONSENSUS_SUMMARY    = "05_Report/consensus_summary.csv"
LQS_CONFIG_TOML          = os.environ.get(
    "LQS_CONFIG_TOML",
    os.environ.get("LQS_CONFIG", "b_low_quality_saver_pipeline_config.toml"),
)
MERGE_CONFIG_TOML        = "a_merge_consensus_config.toml"

# --- output locations ----------------------------------------------------
REPORT_MD_NAME = "Report.md"
ASSETS_SUBDIR  = "assets"            # PNG chromatograms go here, referenced relatively

# --- gel image (optional) -----------------------------------------------
# Globs searched in order, evaluated RELATIVE TO THE REPORT DIR (i.e. inside
# {output_root}/{batch}/full_report/). The first match is embedded in the
# report under a "Gel image" section. Set to [] to disable.
GEL_IMAGE_GLOBS = [
    "assets/gel_image.*",
    "gel_image.*",
    "*.jpg",
    "*.jpeg",
    "*.png",
]
GEL_IMAGE_CAPTION = "PCR gel image (Sanger sequencing submission)."
GEL_IMAGE_MAX_WIDTH = "75%"          # CSS width when rendered to PDF/HTML

# --- PCR DNA quality TSV (Spectrophotometer) ----------------------------
# Workspace-relative path to a TSV with spectrophotometer purity/concentration
# results for the sent samples. Rendered as an HTML table at the end of
# the report (after the gel image) so column widths can be controlled.
# Set to None to disable.
PCR_DNA_QUALITY_TSV = "Sent/PCR_DNA_Quality.tsv"
PCR_DNA_QUALITY_TITLE = "PCR DNA Quality (Spectrophotometer)"
# Per-column widths (CSS %). Last entry should be the Comment column.
# List length must match the number of columns in the TSV header.
PCR_DNA_QUALITY_COL_WIDTHS = ["10%", "9%", "9%", "12%", "60%"]

# --- chromatogram rendering ---------------------------------------------
# Readability is controlled by CHROMATOGRAM_WINDOW (bases per panel) and
# CHROMATOGRAM_WIDTH (inches per panel). Smaller window + wider panel = bigger
# per-base spacing. Approx. inches/base = CHROMATOGRAM_WIDTH / CHROMATOGRAM_WINDOW.
CHROMATOGRAM_DPI    = 130
CHROMATOGRAM_WIDTH  = 18             # inches per panel
CHROMATOGRAM_HEIGHT = 3.6            # inches per panel
CHROMATOGRAM_WINDOW = 80             # bases per panel (was 250 - way too compressed)
MAX_CHROMATOGRAM_PANELS = 24         # cap; covers reads up to ~1920 bases at 80/panel
CHROMATOGRAM_TRACE_LW   = 1.1        # trace line width
CHROMATOGRAM_BASE_FS    = 9          # base-call label font size
CHROMATOGRAM_TITLE_FS   = 10
CHROMATOGRAM_LEGEND_FS  = 9
CHANNEL_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#000000", "T": "#d62728"}

# --- alignment excerpt rendering ----------------------------------------
ALIGNMENT_HEAD_LINES = 80            # show the first N lines of each .ref_align.txt / .fwd_rev_align.txt
ALIGNMENT_TAIL_LINES = 20            # plus the last N lines
EDITS_TABLE_HEAD     = 15            # rows of edits.tsv to embed

# --- verdict thresholds (one of: GOOD / REVIEW / RESEQUENCE) -------------
GOOD_MAX_N_PCT          = 1.0        # consensus N% <= this AND tier == 2_MediumQuality => GOOD
REVIEW_MAX_N_PCT        = 10.0       # else if N% <= this => REVIEW MANUALLY
                                      # else => RESEQUENCE
LOW_TIER_FORCES_RESEQ   = False      # if True, any 3_LowQuality amplicon => RESEQUENCE regardless of N%
NO_BLAST_HIT_IS_RESEQ   = True       # both strands NO_HIT in BLAST => RESEQUENCE
MISSING_CONSENSUS_RESEQ = True       # no 04a consensus produced => RESEQUENCE

# --- manual verdict overrides -------------------------------------------
# (sample, amplicon) pairs the user has manually inspected and approved.
# Forces verdict = GOOD regardless of N% / tier / BLAST.
MANUAL_GOOD_OVERRIDES = {
    ("I_HC01463503", "SmelDMP4_Primer4"),
}

# (sample, amplicon) pairs where only the forward read is salvageable.
# The reverse strand is unreliable, so the merged consensus should NOT be used;
# only the forward .reviewed.fa (or original forward read) is usable downstream.
MANUAL_FORWARD_ONLY_OVERRIDES = {
    ("I_HC01463503", "SmelDMP1_990_Primer1_990"),
    # HC01468262: F score_as_is=564 vs R revcomp=-196; R is junk - use F only.
    ("I_HC01468262", "1_990_1_990"),
}

# (sample, amplicon) pairs where only the reverse read is salvageable.
# The forward strand is unreliable, so only the reverse .reviewed.fa is usable.
MANUAL_REVERSE_ONLY_OVERRIDES = {
    # HC01468262: R revcomp=766 vs F as_is=497 with many mismatches and noisy tail.
    ("I_HC01468262", "2_2"),
}

# --- amplicon display order ---------------------------------------------
# Ordered list of (gene_number, optional_size) tuples. Use None for size when
# the amplicon name has no size suffix (e.g. SmelDMP4_Primer4 vs SmelDMP1_730).
# Amplicons not in this list are appended afterwards in natural-sorted order.
GENE_ORDER = [
    (1, 730),
    (1, 990),
    (4, None),
    (10, 560),
    (12, None),
]

# ===========================================================================
# === implementation =======================================================
# ===========================================================================
import argparse
import csv
import datetime as dt
import importlib
import io
import re
import shutil
import sys
import textwrap
import tomllib
from collections import defaultdict
from typing import Optional

# Import the numbered shared module through its canonical package name.
sys.path.insert(0, str(DEFAULT_WORKSPACE))
try:
    mott_trim = importlib.import_module(
        "scripts.common.03_trim_reads"
    ).mott_trim
except Exception:  # pragma: no cover - degrade gracefully if module missing
    mott_trim = None  # type: ignore


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def read_text_safe(p: Path, max_bytes: int = 2_000_000) -> Optional[str]:
    try:
        data = p.read_bytes()[:max_bytes]
        return data.decode("utf-8", errors="replace")
    except FileNotFoundError:
        return None
    except OSError as e:
        return f"[error reading {p.name}: {e}]"


def head_tail(text: str, head: int, tail: int) -> str:
    lines = text.splitlines()
    if len(lines) <= head + tail:
        return text
    return "\n".join(lines[:head] + [f"... [{len(lines) - head - tail} lines omitted] ..."] + lines[-tail:])


def fenced(content: str, lang: str = "") -> str:
    # Use a 4-backtick fence so embedded triple-backticks (rare here) survive.
    return f"````{lang}\n{content.rstrip()}\n````\n"


def md_escape(s: str) -> str:
    return s.replace("|", r"\|")


def fmt_pct(x) -> str:
    try:
        return f"{float(x):.1f}%"
    except Exception:
        return str(x)


def amplicon_basename(read: str, forward_suffix: str = "_forward", reverse_suffix: str = "_reverse") -> str:
    """Strip the trailing strand suffix to get the amplicon base name.

    Default suffixes are ``_forward`` / ``_reverse`` (the older filename
    convention). Pass per-batch suffixes (e.g. ``F`` / ``R``) for labs that
    use a different naming scheme.
    """
    return re.sub(
        rf"(?:{re.escape(forward_suffix)}|{re.escape(reverse_suffix)})$",
        "",
        read,
    )


def amplicon_sort_key(amp_base: str):
    """Sort key honoring GENE_ORDER. Unknown amplicons sort after, naturally."""
    m = re.match(r"SmelDMP(\d+)(?:_(\d+))?", amp_base)
    if m:
        gene = int(m.group(1))
        size = int(m.group(2)) if m.group(2) else None
        key = (gene, size)
        if key in GENE_ORDER:
            return (0, GENE_ORDER.index(key), amp_base)
        return (1, gene, size if size is not None else -1, amp_base)
    return (2, amp_base)


# --------------------------------------------------------------------------
# data loaders (each returns a plain dict / list; file-missing -> empty)
# --------------------------------------------------------------------------
def load_quality_report(p: Path) -> dict:
    """Returns {(sample, amplicon, direction): row_dict}."""
    out = {}
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            key = (row["sample"], row["gene"], row["direction"])
            out[key] = row
    return out


def load_blast(p: Path) -> dict:
    """Returns {(sample, read): row_dict}."""
    out = {}
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["sample"], row["read"])] = row
    return out


def load_trace_review_summary(p: Path) -> dict:
    """Returns {(sample, read): row_dict}. sample here is the I_<name> form."""
    out = {}
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["sample"], row["read"])] = row
    return out


def load_consensus_summary(p: Path) -> dict:
    """Returns {(sample, amplicon): row_dict} for source==04a_Consensus_TraceReview."""
    out = {}
    if not p.exists():
        return out
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("source") and "TraceReview" not in row["source"]:
                continue
            out[(row["sample"], row["amplicon"])] = row
    return out


def _workspace_path(workspace: Path, fname: str) -> Path:
    p = Path(fname)
    return p if p.is_absolute() else workspace / p


def load_trim_settings(workspace: Path) -> dict:
    """Returns {'lqs': {...}, 'merge': {...}} with keys enabled/q_cutoff/min_len.
    Missing files / sections / keys are filled with defaults (enabled=False, q=0, min_len=50).
    """
    out: dict[str, dict] = {}
    for label, fname in (("lqs", LQS_CONFIG_TOML), ("merge", MERGE_CONFIG_TOML)):
        cfg_path = _workspace_path(workspace, fname)
        trim = {"enabled": False, "q_cutoff": 0, "min_len": 50, "_source": fname}
        if cfg_path.exists():
            try:
                data = tomllib.loads(cfg_path.read_text(encoding="utf-8"))
                t = data.get("trim", {}) or {}
                trim["enabled"]  = bool(t.get("enabled", False))
                trim["q_cutoff"] = int(t.get("q_cutoff", 0))
                trim["min_len"]  = int(t.get("min_len", 50))
            except Exception as e:
                trim["_error"] = str(e)
        out[label] = trim
    return out


def compute_trim_window(ab1: Path, q_cutoff: int, min_len: int) -> Optional[dict]:
    """Run Mott trim on an .ab1; return {'start','end','kept','total','mean_q_kept'} or None on failure."""
    if mott_trim is None or q_cutoff <= 0:
        return None
    try:
        from Bio import SeqIO
        rec = next(SeqIO.parse(str(ab1), "abi"), None)
    except Exception:
        return None
    if rec is None:
        return None
    quals = list(rec.letter_annotations.get("phred_quality", []))
    total = len(quals)
    if total == 0:
        return None
    start, end = mott_trim(quals, q_cutoff, min_len)
    kept = end - start
    mean_q_kept = (sum(quals[start:end]) / kept) if kept > 0 else 0.0
    return {
        "start": start, "end": end, "kept": kept, "total": total,
        "mean_q_kept": mean_q_kept,
    }


def load_merge_summaries(merge_root: Path) -> dict:
    """Walk merge_consensus_pipeline/<tier>/merge_summary.tsv.
    Returns {(sample_no_prefix, amplicon): {"tier": tier, "consensus_fa": rel_path, ...}}.
    """
    out = {}
    if not merge_root.exists():
        return out
    for tier_dir in sorted(merge_root.iterdir()):
        ts = tier_dir / "merge_summary.tsv"
        if not ts.exists():
            continue
        with ts.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["tier"] = tier_dir.name
                out[(row["sample"], row["amplicon"])] = row
    return out


# --------------------------------------------------------------------------
# chromatogram rendering (.ab1 -> PNG) - optional, requires biopython+matplotlib
# --------------------------------------------------------------------------
def render_chromatogram(ab1: Path, out_png: Path) -> Optional[str]:
    """Render a multi-panel chromatogram PNG. Returns None on success, error str on failure."""
    try:
        from Bio import SeqIO
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        return f"missing dependency: {e}"

    try:
        rec = next(SeqIO.parse(str(ab1), "abi"), None)
    except Exception as e:
        return f"could not parse {ab1.name}: {e}"
    if rec is None:
        return f"no records in {ab1.name}"

    raw = rec.annotations.get("abif_raw", {})
    # Channel order in ABIF: DATA9=G, DATA10=A, DATA11=T, DATA12=C
    channels = {
        "G": list(raw.get("DATA9", [])),
        "A": list(raw.get("DATA10", [])),
        "T": list(raw.get("DATA11", [])),
        "C": list(raw.get("DATA12", [])),
    }
    if not any(channels.values()):
        return f"no trace data in {ab1.name}"

    peaks = list(raw.get("PLOC2", raw.get("PLOC1", [])))
    bases = str(rec.seq)
    n_bases = min(len(bases), len(peaks))

    if n_bases == 0:
        return f"no called bases in {ab1.name}"

    n_panels = min(MAX_CHROMATOGRAM_PANELS,
                   max(1, (n_bases + CHROMATOGRAM_WINDOW - 1) // CHROMATOGRAM_WINDOW))

    fig, axes = plt.subplots(n_panels, 1,
                             figsize=(CHROMATOGRAM_WIDTH, CHROMATOGRAM_HEIGHT * n_panels),
                             squeeze=False)
    axes = axes[:, 0]

    for panel_idx, ax in enumerate(axes):
        b_lo = panel_idx * CHROMATOGRAM_WINDOW
        b_hi = min(n_bases, b_lo + CHROMATOGRAM_WINDOW)
        if b_lo >= b_hi:
            ax.set_visible(False)
            continue
        scan_lo = peaks[b_lo]
        scan_hi = peaks[b_hi - 1]
        # add a small margin
        scan_lo = max(0, scan_lo - 5)
        scan_hi = scan_hi + 5

        for base, trace in channels.items():
            ax.plot(range(scan_lo, min(scan_hi, len(trace))),
                    trace[scan_lo:scan_hi],
                    color=CHANNEL_COLORS[base], linewidth=CHROMATOGRAM_TRACE_LW,
                    label=base)

        # base call labels under the x-axis
        ymax = max((max(t[scan_lo:scan_hi]) if t[scan_lo:scan_hi] else 0)
                   for t in channels.values())
        for bi in range(b_lo, b_hi):
            sp = peaks[bi]
            ax.text(sp, -ymax * 0.06, bases[bi],
                    ha="center", va="top",
                    fontsize=CHROMATOGRAM_BASE_FS,
                    color=CHANNEL_COLORS.get(bases[bi].upper(), "#666"))

        ax.set_xlim(scan_lo, scan_hi)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{ab1.name}  bases {b_lo + 1}-{b_hi}",
                     fontsize=CHROMATOGRAM_TITLE_FS, loc="left")
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=8)
        if panel_idx == 0:
            ax.legend(loc="upper right", ncol=4,
                      fontsize=CHROMATOGRAM_LEGEND_FS, frameon=False)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=CHROMATOGRAM_DPI)
    plt.close(fig)
    return None


# --------------------------------------------------------------------------
# verdict logic
# --------------------------------------------------------------------------
def compute_verdict(*, tier: Optional[str], n_pct: Optional[float],
                    has_consensus: bool, both_no_blast_hit: bool,
                    sample: Optional[str] = None,
                    amplicon: Optional[str] = None) -> tuple[str, str]:
    """Returns (verdict_label, reason)."""
    if (sample, amplicon) in MANUAL_GOOD_OVERRIDES:
        n_str = f"{n_pct:.2f}" if n_pct is not None else "-"
        return ("GOOD", f"Manually inspected and approved (N%={n_str}).")
    if (sample, amplicon) in MANUAL_FORWARD_ONLY_OVERRIDES:
        return ("FORWARD_ONLY", "Use the forward .reviewed.fa only.")
    if (sample, amplicon) in MANUAL_REVERSE_ONLY_OVERRIDES:
        return ("REVERSE_ONLY", "Use the reverse .reviewed.fa only.")
    if both_no_blast_hit and NO_BLAST_HIT_IS_RESEQ:
        return ("RESEQUENCE", "Neither strand had a BLAST hit against the reference panel.")
    if not has_consensus and MISSING_CONSENSUS_RESEQ:
        return ("RESEQUENCE", "No merged trace-review consensus produced for this amplicon.")
    if tier == "3_LowQuality" and LOW_TIER_FORCES_RESEQ:
        return ("RESEQUENCE", "Read pair classified as 3_LowQuality.")
    if n_pct is None:
        return ("REVIEW", "Consensus exists but N% could not be determined.")
    if n_pct <= GOOD_MAX_N_PCT and tier == "2_MediumQuality":
        return ("GOOD", f"Consensus N%={n_pct:.2f} and pair classified as MediumQuality.")
    if n_pct <= GOOD_MAX_N_PCT:
        return ("GOOD", f"Consensus N%={n_pct:.2f} (tier={tier}).")
    if n_pct <= REVIEW_MAX_N_PCT:
        return ("REVIEW", f"Consensus N%={n_pct:.2f} - inspect manually before downstream use.")
    return ("RESEQUENCE", f"Consensus N%={n_pct:.2f}, exceeds review threshold of {REVIEW_MAX_N_PCT}%.")


VERDICT_BADGE = {
    "GOOD":         "**[GOOD]** ready for downstream use",
    "FORWARD_ONLY": "**[FORWARD_ONLY]** forward read salvageable; reverse unreliable - do not use merged consensus",
    "REVERSE_ONLY": "**[REVERSE_ONLY]** reverse read salvageable; forward unreliable - do not use merged consensus",
    "REVIEW":       "**[REVIEW]** inspect alignments before use",
    "RESEQUENCE":   "**[RESEQUENCE]** failed quality - resequence amplicon",
}


# --------------------------------------------------------------------------
# main report builder
# --------------------------------------------------------------------------
def build(workspace: Path, render_chromatograms: bool, *,
          batch_id: str, sample_dir_name: str,
          sample_root: str = ".",
          forward_suffix: str = "_forward",
          reverse_suffix: str = "_reverse") -> Path:
    """Build a per-batch Markdown report.

    Path layout (one output tree per Sanger batch):
        workspace/{PIPELINE_OUTPUT_ROOT}/{batch_id}/{LQS_PIPELINE_DIRNAME}/...
        workspace/{PIPELINE_OUTPUT_ROOT}/{batch_id}/{MERGE_PIPELINE_DIRNAME}/...
        workspace/{PIPELINE_OUTPUT_ROOT}/{batch_id}/{REPORT_PIPELINE_DIRNAME}/Report.md

    ``sample_dir_name`` is the on-disk folder under ``sample_root`` (e.g.
    ``I_HC01463503``) whose chromatograms back this batch; only its amplicons
    are summarized so distinct batches never appear in the same report.
    """
    batch_out   = workspace / PIPELINE_OUTPUT_ROOT / batch_id
    lqs_root    = batch_out / LQS_PIPELINE_DIRNAME
    merge_root  = batch_out / MERGE_PIPELINE_DIRNAME
    report_dir  = batch_out / REPORT_PIPELINE_DIRNAME

    qreport     = load_quality_report(lqs_root / LQS_QUALITY_REPORT_CSV)
    blast       = load_blast(lqs_root / LQS_BLAST)
    tr_summary  = load_trace_review_summary(lqs_root / LQS_TRACE_REVIEW_SUMMARY)
    cons_summ   = load_consensus_summary(lqs_root / LQS_CONSENSUS_SUMMARY)
    merge_summ  = load_merge_summaries(merge_root)
    trim_cfg    = load_trim_settings(workspace)

    raw_root = _workspace_path(workspace, sample_root)
    sample_path = Path(sample_dir_name)
    sdir = sample_path if sample_path.is_absolute() else raw_root / sample_path
    if not sdir.is_dir():
        raise SystemExit(f"ERROR: sample dir not found for batch {batch_id}: {sdir}")
    sample_dirs = [sdir]

    assets_dir = report_dir / ASSETS_SUBDIR
    report_dir.mkdir(parents=True, exist_ok=True)
    if render_chromatograms:
        assets_dir.mkdir(parents=True, exist_ok=True)

    md = io.StringIO()
    w = lambda s="": md.write(s + "\n")

    # --- top-level header --------------------------------------------------
    w(f"# PCR_Sequence Full Pipeline Report")
    w()
    w("## Pipeline overview")
    w()
    w("Raw Sanger `.ab1` traces are processed to produce one consensus per amplicon")
    w("in two passes:")
    w()
    w("**Pass 1 - direct merge (`merge_consensus_pipeline`).** Reads are scored by Phred")
    w("quality and classified into `1_HighQuality`, `2_MediumQuality`, or `3_LowQuality`.")
    w("For each amplicon, the forward and reverse reads are aligned and merged into one")
    w("consensus FASTA, with an accompanying alignment report.")
    w()
    w("**Pass 2 - rescue of weaker reads (`low_quality_saver_pipeline`).** Each strand is")
    w("BLASTed against the reference panel (`PF05078_DMP_transcripts_cdhit.fa`) to")
    w("identify the most likely reference. The `trace_review` step then re-calls noisy")
    w("bases against that reference, producing `*.reviewed.fa`, `*.ref_align.txt`, and")
    w("`*.edits.tsv` outputs. The two reviewed strands are merged to form the 04a")
    w("reviewed consensus, and ambiguous `N` calls are tallied.")
    w()
    w("**Verdict** (based on the reviewed consensus):")
    w()
    w(f"- `GOOD` - N% \u2264 {GOOD_MAX_N_PCT:g}%.")
    w(f"- `REVIEW` - N% \u2264 {REVIEW_MAX_N_PCT:g}%.")
    w("- `RESEQUENCE` - amplicon has too many `N`s, no consensus was produced, or")
    w("  neither strand had a BLAST hit against the reference panel.")
    w()

    # --- end-trimming settings (Mott's modified algorithm) ---------------
    w("## End-trimming settings")
    w()
    w("Both pipelines use Mott's modified trimming algorithm to remove low-quality")
    w("flanks before downstream steps. `merge_consensus_pipeline` slices/masks the reads")
    w("fed into `pairwise_consensus`; `low_quality_saver_pipeline` masks the reads fed")
    w("into BLAST, `trace_review`, and `consensus_assembly`. Note that `tracy_decompose`")
    w("retains its own trimming, and `step_merge_consensus` only operates on already-called")
    w("FASTAs.")
    w()
    w("| Pipeline | Config file | Enabled | Q cutoff | Min kept length |")
    w("|---|---|:---:|---:|---:|")
    for label, lbl_disp in (("merge", "merge_consensus_pipeline"),
                            ("lqs",   "low_quality_saver_pipeline")):
        t = trim_cfg.get(label, {})
        en = "yes" if t.get("enabled") else "no"
        q  = t.get("q_cutoff", "-")
        ml = t.get("min_len", "-")
        src = t.get("_source", "-")
        w(f"| {lbl_disp} | `{src}` | {en} | {q} | {ml} |")
    w()

    # --- summary rows collected per amplicon below; rendered after the loop
    summary_rows: list[dict] = []

    detail_chunks: list[str] = []

    for sdir in sample_dirs:
        sample_with_prefix = sdir.name              # e.g. I_HC01463503
        sample_clean = re.sub(r"^I_", "", sample_with_prefix)

        # collect amplicons by scanning forward+reverse .ab1 pairs in sample dir
        amplicons: dict[str, dict] = defaultdict(dict)
        for ab1 in sorted(sdir.glob("*.ab1")):
            stem = ab1.stem
            base = amplicon_basename(stem, forward_suffix, reverse_suffix)
            if stem.endswith(forward_suffix):
                amplicons[base]["forward_ab1"] = ab1
            elif stem.endswith(reverse_suffix):
                amplicons[base]["reverse_ab1"] = ab1

        sample_chunk = io.StringIO()
        sw = lambda s="": sample_chunk.write(s + "\n")
        sw(f"## Sample `{sample_with_prefix}`")
        sw()

        for amp_base, files in sorted(amplicons.items(), key=lambda kv: amplicon_sort_key(kv[0])):
            sw(f"### Amplicon `{amp_base}`")
            sw()

            # -- compute trim windows once per strand (used in tables below) --
            # Use the LQS trim setting for per-strand display; both pipelines now use
            # identical defaults (q_cutoff=20, min_len=50). If LQS is disabled, fall
            # back to merge settings so the report still shows something.
            _t_lqs   = trim_cfg.get("lqs",   {}) or {}
            _t_merge = trim_cfg.get("merge", {}) or {}
            _t_active = _t_lqs if _t_lqs.get("enabled") else _t_merge
            _q  = int(_t_active.get("q_cutoff", 0)) if _t_active.get("enabled") else 0
            _ml = int(_t_active.get("min_len", 50))
            trim_windows: dict[str, Optional[dict]] = {}
            for _label in ("forward", "reverse"):
                _ab1 = files.get(f"{_label}_ab1")
                trim_windows[_label] = compute_trim_window(_ab1, _q, _ml) if _ab1 else None

            # -- 1) Quality classification ----------------------------------
            f_q = qreport.get((sample_clean, amp_base, "forward"))
            r_q = qreport.get((sample_clean, amp_base, "reverse"))
            tier = (f_q or r_q or {}).get("tier")
            sw("**1. Quality classification**")
            sw()
            sw("| Strand | Bases | Mean Phred | %Q20 | %Q30 | Tier | Trim window | Kept | Mean Q kept |")
            sw("|---|---:|---:|---:|---:|---|---|---:|---:|")
            for label, row in [("forward", f_q), ("reverse", r_q)]:
                tw = trim_windows.get(label)
                if tw is not None:
                    tw_cell  = f"[{tw['start']},{tw['end']})/{tw['total']}"
                    kept_cell = f"{tw['kept']}"
                    mq_cell  = f"{tw['mean_q_kept']:.1f}"
                else:
                    tw_cell = kept_cell = mq_cell = "-"
                if row is None:
                    sw(f"| {label} | (missing) | - | - | - | - | {tw_cell} | {kept_cell} | {mq_cell} |")
                else:
                    sw(f"| {label} | {row['bases']} | {row['mean_phred']} "
                       f"| {fmt_pct(row['pct_q20'])} | {fmt_pct(row['pct_q30'])} "
                       f"| {md_escape(row['tier'])} | {tw_cell} | {kept_cell} | {mq_cell} |")
            sw()

            # -- BLAST assignment summary -----------------------------------
            f_b = blast.get((sample_with_prefix, f"{amp_base}{forward_suffix}"))
            r_b = blast.get((sample_with_prefix, f"{amp_base}{reverse_suffix}"))
            sw("**BLAST assignment**")
            sw()
            sw("| Strand | Best ref | %Identity | %QCov | E-value |")
            sw("|---|---|---:|---:|---:|")
            for label, row in [("forward", f_b), ("reverse", r_b)]:
                if row is None:
                    sw(f"| {label} | (no row) | - | - | - |")
                else:
                    sw(f"| {label} | {md_escape(row['best_ref'])} | {row['pident']} "
                       f"| {row['qcov']} | {row['evalue']} |")
            sw()
            both_no_hit = (f_b and r_b and f_b["best_ref"] == "NO_HIT" and r_b["best_ref"] == "NO_HIT")

            # -- 2) Original chromatograms ---------------------------------
            sw("**2. Original chromatograms**")
            sw()
            for label in ("forward", "reverse"):
                ab1 = files.get(f"{label}_ab1")
                if ab1 is None:
                    sw(f"_{label}: missing .ab1_  ")
                    continue
                if render_chromatograms:
                    png = assets_dir / f"{sample_with_prefix}__{ab1.stem}.png"
                    err = render_chromatogram(ab1, png)
                    if err:
                        sw(f"_{label} chromatogram render failed: {err}_  ")
                        sw(f"_source: `{ab1.relative_to(workspace)}`_  ")
                    else:
                        rel = png.relative_to(report_dir).as_posix()
                        sw(f"_{label}_ (`{ab1.name}`):")
                        sw()
                        sw(f"![{label} chromatogram of {amp_base}]({rel})")
                        sw()
                else:
                    sw(f"_{label}_ source: `{ab1.relative_to(workspace)}` (chromatogram rendering disabled)  ")
            sw()

            # -- 3) merge_consensus_pipeline result ------------------------
            ms = merge_summ.get((sample_clean, amp_base))
            sw("**3. `merge_consensus_pipeline` result**")
            sw()
            if ms is None:
                sw("_No `merge_consensus_pipeline` entry for this amplicon._")
                sw()
            else:
                tier_dir = merge_root / ms["tier"] / sample_with_prefix
                cons_fa  = tier_dir / f"{amp_base}_consensus.fa"
                rep_tsv  = tier_dir / f"{amp_base}_consensus.report.tsv"
                aln_txt  = tier_dir / f"{amp_base}_consensus.fwd_rev_align.txt"
                sw(f"- tier: `{ms['tier']}`")
                sw(f"- status: `{ms.get('status', '')}`")
                sw(f"- consensus FASTA: `{cons_fa.relative_to(workspace).as_posix()}`")
                sw()
                rep = read_text_safe(rep_tsv)
                if rep:
                    sw("Per-pair report:")
                    sw(fenced(rep, "tsv"))
                cons_text = read_text_safe(cons_fa)
                if cons_text:
                    sw("Consensus FASTA:")
                    sw(fenced(cons_text, "fasta"))
                aln = read_text_safe(aln_txt)
                if aln:
                    sw("Forward vs reverse alignment:")
                    sw(fenced(aln, "text"))

            # -- 4) trace_review per strand --------------------------------
            tr_dir = lqs_root / LQS_TRACE_REVIEW_DIR / sample_with_prefix
            sw("**4. `trace_review` per-strand review**")
            sw()
            for label, _sfx in (("forward", forward_suffix), ("reverse", reverse_suffix)):
                read_id = f"{amp_base}{_sfx}"
                tr_row = tr_summary.get((sample_with_prefix, read_id))
                reviewed = tr_dir / f"{read_id}.reviewed.fa"
                ref_aln  = tr_dir / f"{read_id}.ref_align.txt"
                edits    = tr_dir / f"{read_id}.edits.tsv"
                sw(f"_{label}_:  ")
                if tr_row:
                    sw(f"- ref: `{tr_row.get('ref', '')}`  edits: `{tr_row.get('edits', '')}`  "
                       f"lowQ_frac: `{tr_row.get('lowQ_frac', '')}`  status: `{tr_row.get('status', '')}`")
                else:
                    sw("- _no trace_review summary row_")
                rv = read_text_safe(reviewed)
                if rv:
                    sw("Reviewed FASTA:")
                    sw(fenced(rv, "fasta"))
                ra = read_text_safe(ref_aln)
                if ra:
                    ref_name = (tr_row or {}).get("ref", "") or "unknown"
                    sw(f"Reference alignment vs `{ref_name}`:")
                    sw(fenced(ra, "text"))
                sw()

            # -- 5) attempted merged trace_review consensus (04a) ----------
            sw("**5. Attempted merged `trace_review` consensus (04a)**")
            sw()
            cons_row = cons_summ.get((sample_with_prefix, amp_base))
            cons_path = lqs_root / LQS_CONSENSUS_TR_DIR / sample_with_prefix / f"{amp_base}_consensus.fa"
            n_pct: Optional[float] = None
            if cons_row:
                try:
                    n_pct = float(cons_row.get("N_pct"))
                except (TypeError, ValueError):
                    n_pct = None
                sw(f"- length: `{cons_row.get('length', '')}`  N_count: `{cons_row.get('N_count', '')}`  "
                   f"N%: `{cons_row.get('N_pct', '')}`")
            else:
                sw("- _no entry in consensus_summary.csv_")
            cons_txt = read_text_safe(cons_path)
            if cons_txt:
                sw("Consensus FASTA:")
                sw(fenced(cons_txt, "fasta"))
            else:
                sw(f"- _no consensus file at `{cons_path.relative_to(workspace).as_posix()}`_")
            sw()

            # -- 6) verdict -------------------------------------------------
            verdict, reason = compute_verdict(
                tier=tier,
                n_pct=n_pct,
                has_consensus=cons_txt is not None,
                both_no_blast_hit=bool(both_no_hit),
                sample=sample_with_prefix,
                amplicon=amp_base,
            )
            sw("**6. Verdict**")
            sw()
            sw(f"{VERDICT_BADGE[verdict]}  ")
            sw(f"_Reason: {reason}_")
            sw()
            sw("---")
            sw()

            summary_rows.append({
                "sample": sample_with_prefix,
                "amplicon": amp_base,
                "tier": tier or "-",
                "n_pct": "-" if n_pct is None else f"{n_pct:.2f}",
                "verdict": verdict,
                "reason": reason,
                # quality
                "f_bases":  (f_q or {}).get("bases", "-"),
                "r_bases":  (r_q or {}).get("bases", "-"),
                "f_meanq":  (f_q or {}).get("mean_phred", "-"),
                "r_meanq":  (r_q or {}).get("mean_phred", "-"),
                "f_q30":    fmt_pct((f_q or {}).get("pct_q30")) if f_q else "-",
                "r_q30":    fmt_pct((r_q or {}).get("pct_q30")) if r_q else "-",
                # trim windows (per strand) - compact "lo-hi"; total shown in F/R bases columns
                "f_trim":   (f"{trim_windows['forward']['start']}-{trim_windows['forward']['end']}")
                            if trim_windows.get("forward") else "-",
                "r_trim":   (f"{trim_windows['reverse']['start']}-{trim_windows['reverse']['end']}")
                            if trim_windows.get("reverse") else "-",
                # blast
                "f_ref":    (f_b or {}).get("best_ref", "-"),
                "r_ref":    (r_b or {}).get("best_ref", "-"),
                "f_pident": (f_b or {}).get("pident", "-"),
                "r_pident": (r_b or {}).get("pident", "-"),
                # merge_consensus_pipeline
                "merge_status":  (ms or {}).get("status", "-") if ms else "-",
                "merge_length":  None,
                "merge_overlap": None,
                "merge_agree":   None,
                "merge_disagree": None,
                "merge_n":       None,
                # trace_review
                "f_tr_edits":   (tr_summary.get((sample_with_prefix, f"{amp_base}{forward_suffix}")) or {}).get("edits", "-"),
                "r_tr_edits":   (tr_summary.get((sample_with_prefix, f"{amp_base}{reverse_suffix}")) or {}).get("edits", "-"),
                "f_tr_lowq":    (tr_summary.get((sample_with_prefix, f"{amp_base}{forward_suffix}")) or {}).get("lowQ_frac", "-"),
                "r_tr_lowq":    (tr_summary.get((sample_with_prefix, f"{amp_base}{reverse_suffix}")) or {}).get("lowQ_frac", "-"),
                # final consensus
                "cons_length":  (cons_row or {}).get("length", "-") if cons_row else "-",
                "cons_n_count": (cons_row or {}).get("N_count", "-") if cons_row else "-",
                "has_consensus": cons_txt is not None,
            })

            # back-fill merge stats from the per-pair report TSV if available
            if ms is not None:
                _rep_path = merge_root / ms["tier"] / sample_with_prefix / f"{amp_base}_consensus.report.tsv"
                _rep_text = read_text_safe(_rep_path)
                if _rep_text:
                    _kv: dict[str, str] = {}
                    for _line in _rep_text.splitlines():
                        if "\t" in _line:
                            _k, _v = _line.split("\t", 1)
                            _kv[_k.strip()] = _v.strip()
                    last = summary_rows[-1]
                    last["merge_length"]   = _kv.get("length")
                    last["merge_overlap"]  = _kv.get("overlap")
                    last["merge_agree"]    = _kv.get("agree")
                    last["merge_disagree"] = _kv.get("disagree")
                    last["merge_n"]        = _kv.get("N_calls")

        detail_chunks.append(sample_chunk.getvalue())

    # --- write top summary tables ----------------------------------------
    w("## Summary")
    w()
    if not summary_rows:
        w("_No amplicons found._")
        w()
    else:
        # Overall verdict counts
        counts: dict[str, int] = defaultdict(int)
        for r in summary_rows:
            counts[r["verdict"]] += 1
        n_total = len(summary_rows)
        n_samples = len({r["sample"] for r in summary_rows})
        w(f"**Totals:** {n_total} amplicon(s) across {n_samples} sample(s) - "
          + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
        w()

        # 1. Verdict overview
        w("### 1. Verdict overview")
        w()
        w("| Sample | Amplicon | Tier | Consensus N% | Verdict | Reason |")
        w("|---|---|---|---:|---|---|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} | {md_escape(r['tier'])} "
              f"| {r['n_pct']} | **{r['verdict']}** | {md_escape(r['reason'])} |")
        w()

        # 2. Per-strand quality classification (split F/R for readability)
        w("### 2. Per-strand quality (Phred)")
        w()
        w("**Forward strand**")
        w()
        w("| Sample | Amplicon | Tier | bases | meanQ | %Q30 | kept |")
        w("|---|---|---|---:|---:|---:|---|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} | {md_escape(r['tier'])} "
              f"| {r['f_bases']} | {r['f_meanq']} | {r['f_q30']} | {r['f_trim']} |")
        w()
        w("**Reverse strand**")
        w()
        w("| Sample | Amplicon | Tier | bases | meanQ | %Q30 | kept |")
        w("|---|---|---|---:|---:|---:|---|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} | {md_escape(r['tier'])} "
              f"| {r['r_bases']} | {r['r_meanq']} | {r['r_q30']} | {r['r_trim']} |")
        w()

        # 3. BLAST assignment
        w("### 3. BLAST assignment vs reference panel")
        w()
        w("| Sample | Amplicon | F best ref | F %ident | R best ref | R %ident |")
        w("|---|---|---|---:|---|---:|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} "
              f"| {md_escape(str(r['f_ref']))} | {r['f_pident']} "
              f"| {md_escape(str(r['r_ref']))} | {r['r_pident']} |")
        w()

        # 4. merge_consensus_pipeline result
        w("### 4. `merge_consensus_pipeline` (forward+reverse direct merge)")
        w()
        w("| Sample | Amplicon | Status | Length | Overlap | Agree | Disagree | N calls |")
        w("|---|---|---|---:|---:|---:|---:|---:|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} | {md_escape(str(r['merge_status']))} "
              f"| {r['merge_length'] or '-'} | {r['merge_overlap'] or '-'} "
              f"| {r['merge_agree'] or '-'} | {r['merge_disagree'] or '-'} "
              f"| {r['merge_n'] or '-'} |")
        w()

        # 5. trace_review per-strand edits
        w("### 5. `trace_review` per-strand edits")
        w()
        w("| Sample | Amplicon | F edits | F lowQ_frac | R edits | R lowQ_frac |")
        w("|---|---|---:|---:|---:|---:|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} "
              f"| {r['f_tr_edits']} | {r['f_tr_lowq']} "
              f"| {r['r_tr_edits']} | {r['r_tr_lowq']} |")
        w()

        # 6. Final merged trace_review consensus (04a)
        w("### 6. Final merged `trace_review` consensus (04a)")
        w()
        w("| Sample | Amplicon | Has FASTA | Length | N count | N% | Verdict |")
        w("|---|---|:---:|---:|---:|---:|---|")
        for r in summary_rows:
            w(f"| {r['sample']} | {r['amplicon']} "
              f"| {'yes' if r['has_consensus'] else 'no'} "
              f"| {r['cons_length']} | {r['cons_n_count']} | {r['n_pct']} "
              f"| **{r['verdict']}** |")
        w()

    w("---")
    w()
    w("## Per-amplicon details")
    w()

    # --- append detail sections ------------------------------------------
    for chunk in detail_chunks:
        md.write(chunk)

    # --- gel image (optional, at end) ------------------------------------
    # Globs are evaluated under this batch's report dir so each batch can ship
    # its own gel image (or none) without colliding with the other batches.
    gel_path = None
    for _glob in GEL_IMAGE_GLOBS:
        _matches = sorted(report_dir.glob(_glob))
        # exclude per-strand chromatogram PNGs that live in assets/
        _matches = [p for p in _matches if "__" not in p.name]
        if _matches:
            gel_path = _matches[0]
            break
    if gel_path is not None:
        assets_dir.mkdir(parents=True, exist_ok=True)
        gel_dst = assets_dir / f"gel_image{gel_path.suffix.lower()}"
        try:
            if not gel_dst.exists() or gel_dst.stat().st_mtime < gel_path.stat().st_mtime:
                shutil.copy2(gel_path, gel_dst)
        except Exception as _e:
            print(f"[warn] could not copy gel image {gel_path} -> {gel_dst}: {_e}")
        rel = f"{ASSETS_SUBDIR}/{gel_dst.name}"
        w("---")
        w()
        w("## Gel Image of Sent Samples (May 13, Wednesday)")
        w()
        w(f'<img src="{rel}" alt="Gel image" style="width:{GEL_IMAGE_MAX_WIDTH};max-width:100%;" />')
        w()
        w(f"*{GEL_IMAGE_CAPTION}*")
        w()
        w("**Recommendation:** Purification.")
        w()

    # --- PCR DNA quality (spectrophotometer) TSV -------------------------
    if PCR_DNA_QUALITY_TSV:
        _tsv_path = workspace / PCR_DNA_QUALITY_TSV
        if _tsv_path.is_file():
            try:
                _rows = [
                    line.rstrip("\n").split("\t")
                    for line in _tsv_path.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ]
            except Exception as _e:
                print(f"[warn] could not read PCR DNA quality TSV {_tsv_path}: {_e}")
                _rows = []
            if len(_rows) >= 2:
                w("---")
                w()
                w(f"## {PCR_DNA_QUALITY_TITLE}")
                w()
                _hdr = _rows[0]
                _widths = (PCR_DNA_QUALITY_COL_WIDTHS + [""] * len(_hdr))[:len(_hdr)]

                def _esc(_s: str) -> str:
                    return (_s.replace("&", "&amp;")
                              .replace("<", "&lt;")
                              .replace(">", "&gt;"))

                w('<table style="width:100%;border-collapse:collapse;" border="1">')
                w("  <colgroup>")
                for _wd in _widths:
                    w(f'    <col style="width:{_wd};">' if _wd else "    <col>")
                w("  </colgroup>")
                w("  <thead>")
                w("    <tr>" + "".join(f"<th>{_esc(c)}</th>" for c in _hdr) + "</tr>")
                w("  </thead>")
                w("  <tbody>")
                for _r in _rows[1:]:
                    _r = (_r + [""] * len(_hdr))[:len(_hdr)]
                    w("    <tr>" + "".join(f"<td>{_esc(c)}</td>" for c in _r) + "</tr>")
                w("  </tbody>")
                w("</table>")
                w()

    out_md = report_dir / REPORT_MD_NAME
    out_md.write_text(md.getvalue(), encoding="utf-8")
    return out_md


def _load_batches_from_toml(workspace: Path) -> list[dict]:
    """Return the [[batches]] table from low_quality_saver_pipeline_config.toml.

    Each entry is normalized into a dict with id, sample_dir, forward_suffix,
    reverse_suffix (defaults filled from [reads] when not overridden). Used by
    main() to drive one report per batch from a single CLI invocation.
    """
    cfg_path = _workspace_path(workspace, LQS_CONFIG_TOML)
    if not cfg_path.is_file():
        return []
    data = tomllib.loads(cfg_path.read_text(encoding="utf-8"))
    project = data.get("project", {}) or {}
    sample_root = str(project.get("sample_root", "."))
    reads = data.get("reads", {}) or {}
    default_fwd = str(reads.get("forward_suffix", "_forward"))
    default_rev = str(reads.get("reverse_suffix", "_reverse"))
    # Honor the same active_batches gate the bash drivers use, so a single
    # commented-out entry in the TOML silences the batch everywhere.
    _active = (
        (data.get("project") or {}).get("active_batches")
        or (data.get("io") or {}).get("active_batches")
        or data.get("active_batches")
        or None
    )
    active_set = {str(x) for x in _active} if _active else None
    out: list[dict] = []
    for b in data.get("batches", []) or []:
        if not isinstance(b, dict):
            continue
        if "id" not in b or "sample_dir" not in b:
            continue
        if active_set is not None and str(b["id"]) not in active_set:
            continue
        out.append({
            "id":             str(b["id"]),
            "sample_dir":     str(b["sample_dir"]),
            "sample_root":    sample_root,
            "forward_suffix": str(b.get("forward_suffix", default_fwd)),
            "reverse_suffix": str(b.get("reverse_suffix", default_rev)),
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the one-big PCR_Sequence pipeline Markdown report.")
    ap.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE,
                    help="PCR_Sequence workspace root (default: inferred from this script).")
    ap.add_argument("--no-chromatograms", action="store_true",
                    help="Skip rendering .ab1 -> PNG (faster, text-only report).")
    ap.add_argument("--batch", action="append", default=None,
                    help="Build only this batch id (repeatable). "
                         "Default: iterate every [[batches]] entry from the LQS config.")
    ap.add_argument("--sample-dir",
                    help="Override the sample dir for a single-batch ad-hoc run "
                         "(requires exactly one --batch).")
    ap.add_argument("--forward-suffix",
                    help="Override forward-strand stem suffix for an ad-hoc single-batch run.")
    ap.add_argument("--reverse-suffix",
                    help="Override reverse-strand stem suffix for an ad-hoc single-batch run.")
    args = ap.parse_args()

    workspace = args.workspace.resolve()
    if not workspace.is_dir():
        print(f"ERROR: workspace not a directory: {workspace}", file=sys.stderr)
        return 2

    batches = _load_batches_from_toml(workspace)
    if args.batch:
        wanted = set(args.batch)
        batches = [b for b in batches if b["id"] in wanted]
        # Ad-hoc single-batch overrides take precedence so the user can run
        # without editing the TOML (e.g. a one-off prefix change).
        if args.sample_dir:
            if len(batches) != 1:
                print("ERROR: --sample-dir requires exactly one matching --batch", file=sys.stderr)
                return 2
            batches[0]["sample_dir"] = args.sample_dir
        if args.forward_suffix:
            for b in batches:
                b["forward_suffix"] = args.forward_suffix
        if args.reverse_suffix:
            for b in batches:
                b["reverse_suffix"] = args.reverse_suffix
    if not batches:
        print("ERROR: no batches resolved -- check [[batches]] in "
              f"{LQS_CONFIG_TOML} or pass --batch/--sample-dir.", file=sys.stderr)
        return 2

    rc = 0
    for b in batches:
        print(f"[report] === batch {b['id']} (sample_dir={b['sample_dir']}, sample_root={b['sample_root']}) ===")
        try:
            out_md = build(
                workspace,
                render_chromatograms=not args.no_chromatograms,
                batch_id        = b["id"],
                sample_dir_name = b["sample_dir"],
                sample_root     = b["sample_root"],
                forward_suffix  = b["forward_suffix"],
                reverse_suffix  = b["reverse_suffix"],
            )
            print(f"[report] wrote {out_md}")
        except SystemExit as e:
            print(f"[report] batch {b['id']} FAILED: {e}", file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
