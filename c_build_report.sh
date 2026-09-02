#!/usr/bin/env bash
# =============================================================================
# c_build_report.sh - run the full Markdown report builder (WSL).
#
# Usage (from this directory):
#   ./c_build_report.sh                  # build report with chromatograms
#   ./c_build_report.sh --no-chromatograms
#
# Edit BUILD_PDF / PDF_ENGINE below to control PDF conversion.
# Override Python via:  PYTHON=/path/to/python3.11+
# =============================================================================
set -Eeuo pipefail

# -------------------- USER-EDITABLE CONFIG ----------------------------------
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPORT_PY="${REPORT_PY:-$ROOT/scripts/reporting/01_build_full_report.py}"
LOG_DIR="${LOG_DIR:-$ROOT/logs}"
# 01_build_full_report.py reads its batch list from this TOML's [[batches]] table
# unless --batch is passed on the command line. Each batch gets its own report.
LQS_CONFIG="${LQS_CONFIG:-$ROOT/b_low_quality_saver_pipeline_config.toml}"
READ_TOML="${READ_TOML:-$ROOT/scripts/common/01_read_toml.py}"

# Per-batch report dir layout (kept in sync with 01_build_full_report.py):
#   III_PipelineOutput/<batch_id>/full_report/Report.md
PIPELINE_OUTPUT_ROOT="III_PipelineOutput"
REPORT_SUBDIR="full_report"

# Convert Report.md -> Report.pdf after building? true / false.
BUILD_PDF=true
# Pandoc PDF engine: weasyprint | wkhtmltopdf | xelatex | pdflatex
# Leave empty ("") to let pandoc auto-pick.
PDF_ENGINE="weasyprint"

# PDF layout (used when PDF_ENGINE=weasyprint).
PDF_FONT_SIZE="9pt"
PDF_MARGIN="0.5in"

# -------------------- pick a Python 3.11+ interpreter -----------------------
_pick_python() {
  local cand
  for cand in "${PYTHON:-}" \
              "$HOME/miniconda3/bin/python3" "$HOME/anaconda3/bin/python3" \
              python3.13 python3.12 python3.11 python3 /usr/bin/python3; do
    [[ -z "$cand" ]] && continue
    command -v "$cand" >/dev/null 2>&1 || continue
    "$cand" -c 'import sys, Bio' >/dev/null 2>&1 || continue
    "$cand" -c 'import sys; sys.exit(0 if sys.version_info>=(3,11) else 1)' 2>/dev/null \
      && { echo "$cand"; return 0; }
  done
  return 1
}

PY="$(_pick_python)" || {
  echo "ERROR: need Python >= 3.11 with Biopython." >&2
  echo "  Set PYTHON=/path/to/python3.11+ and re-run." >&2
  exit 2
}

mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/$(date +%Y%m%d_%H%M%S)_full_report.log"

echo "[report] python:    $PY ($("$PY" --version 2>&1))"
echo "[report] build pdf: $BUILD_PDF${PDF_ENGINE:+ (engine=$PDF_ENGINE)}"
echo "[report] log:       $LOG"
echo

"$PY" "$REPORT_PY" --workspace "$ROOT" "$@" 2>&1 | tee "$LOG"

# -------------------- resolve which batches were built ----------------------
# Honor --batch on the CLI (mirrors 01_build_full_report.py); otherwise iterate
# every [[batches]] entry from the LQS config so PDF conversion covers them all.
CLI_BATCHES=()
for arg in "$@"; do
  case "$arg" in
    --batch=*) CLI_BATCHES+=("${arg#--batch=}") ;;
  esac
done
# Also support `--batch ID` (two-token form) by re-scanning the arg list.
_prev=""
for arg in "$@"; do
  if [[ "$_prev" == "--batch" ]]; then
    CLI_BATCHES+=("$arg")
  fi
  _prev="$arg"
done

if (( ${#CLI_BATCHES[@]} > 0 )); then
  BATCH_IDS=("${CLI_BATCHES[@]}")
else
  # Read [[batches]].id column from the TOML.
  mapfile -t BATCH_IDS < <("$PY" "$READ_TOML" "$LQS_CONFIG" --batches | awk 'NR>1 {print $1}')
fi

if (( ${#BATCH_IDS[@]} == 0 )); then
  echo "[report] no batches resolved -- nothing to PDF-convert" | tee -a "$LOG"
  exit 0
fi

# -------------------- optional Markdown -> PDF conversion (per batch) -------
if [[ "${BUILD_PDF,,}" != "true" ]]; then
  for bid in "${BATCH_IDS[@]}"; do
    md="$ROOT/$PIPELINE_OUTPUT_ROOT/$bid/$REPORT_SUBDIR/Report.md"
    echo "[report] batch $bid -> $md"
  done
  exit 0
fi

if ! command -v pandoc >/dev/null 2>&1; then
  echo "[pdf] skipped: pandoc not installed (apt install pandoc)" | tee -a "$LOG"
  exit 0
fi

# Resolve a working PDF engine ONCE (used for every batch). The bare
# `weasyprint` on PATH may belong to a stale env (e.g. python 3.5) and crash
# with SyntaxError. Prefer running weasyprint as a module from the same Python
# we used for the report builder.
ENGINE_NAME=""
if [[ "$PDF_ENGINE" == "weasyprint" ]]; then
  if "$PY" -c 'import weasyprint' >/dev/null 2>&1; then
    WRAPPER_DIR="$LOG_DIR/.pdf_engine_bin"
    mkdir -p "$WRAPPER_DIR"
    WRAPPER="$WRAPPER_DIR/weasyprint"
    cat > "$WRAPPER" <<EOF
#!/usr/bin/env bash
exec "$PY" -m weasyprint "\$@"
EOF
    chmod +x "$WRAPPER"
    export PATH="$WRAPPER_DIR:$PATH"
    ENGINE_NAME="weasyprint"
    echo "[pdf] using weasyprint via: $PY -m weasyprint" | tee -a "$LOG"
  elif command -v weasyprint >/dev/null 2>&1 && weasyprint --version >/dev/null 2>&1; then
    ENGINE_NAME="weasyprint"
  else
    echo "[pdf] weasyprint not usable. Install it in the report Python env:" | tee -a "$LOG"
    echo "[pdf]   $PY -m pip install weasyprint" | tee -a "$LOG"
    exit 1
  fi
elif [[ -n "$PDF_ENGINE" ]]; then
  ENGINE_NAME="$PDF_ENGINE"
fi

# Convert each batch's Report.md -> Report.pdf in its own directory so that
# assets paths (chromatograms, gel image) resolve against the report dir.
FAILED_BATCHES=()
for bid in "${BATCH_IDS[@]}"; do
  REPORT_MD="$ROOT/$PIPELINE_OUTPUT_ROOT/$bid/$REPORT_SUBDIR/Report.md"
  if [[ ! -f "$REPORT_MD" ]]; then
    echo "[pdf] batch $bid: $REPORT_MD not found -- skipping" | tee -a "$LOG"
    continue
  fi

  REPORT_PDF="${REPORT_MD%.md}.pdf"
  REPORT_DIR="$(dirname "$REPORT_MD")"
  REPORT_MD_NAME="$(basename "$REPORT_MD")"
  REPORT_PDF_NAME="$(basename "$REPORT_PDF")"

  REPORT_CSS="$REPORT_DIR/.report_style.css"
  cat > "$REPORT_CSS" <<EOF
@page { size: A4; margin: ${PDF_MARGIN}; }
html { font-size: ${PDF_FONT_SIZE}; }
body { font-family: "DejaVu Sans", "Helvetica", sans-serif; line-height: 1.3; }
h1 { font-size: 1.6em; }
h2 { font-size: 1.3em; }
h3 { font-size: 1.1em; }
h4, h5, h6 { font-size: 1em; }
pre, code { font-family: "DejaVu Sans Mono", monospace; font-size: 0.85em; }
pre { white-space: pre-wrap; overflow-wrap: break-word; }
table { font-size: 0.9em; border-collapse: collapse; margin: 0.5em 0; }
th, td { padding: 3px 6px; border: 1px solid #888; }
th { background: #e8e8e8; text-align: left; }
tr:nth-child(even) td { background: #f6f6f6; }
img { max-width: 100%; }
EOF

  REPORT_TPL="$REPORT_DIR/.report_template.html"
  cat > "$REPORT_TPL" <<'EOF'
<!DOCTYPE html>
<html lang="$lang$">
<head>
  <meta charset="utf-8" />
  <title>$if(title-prefix)$$title-prefix$ - $endif$$pagetitle$</title>
$for(css)$
  <link rel="stylesheet" href="$css$" />
$endfor$
</head>
<body>
$if(title)$<h1 class="title">$title$</h1>$endif$
$body$
</body>
</html>
EOF

  PANDOC_ARGS=(
    "$REPORT_MD_NAME" -o "$REPORT_PDF_NAME"
    --resource-path="$REPORT_DIR"
    --metadata pagetitle="PCR_Sequence Full Pipeline Report ($bid)"
    --template="$(basename "$REPORT_TPL")"
    --css="$(basename "$REPORT_CSS")"
  )
  [[ -n "$ENGINE_NAME" ]] && PANDOC_ARGS+=(--pdf-engine="$ENGINE_NAME")

  echo "[pdf] batch $bid: (cd $REPORT_DIR && pandoc ${PANDOC_ARGS[*]})" | tee -a "$LOG"
  if ( cd "$REPORT_DIR" && pandoc "${PANDOC_ARGS[@]}" ) 2>&1 | tee -a "$LOG"; then
    echo "[pdf] batch $bid: done -> $REPORT_PDF" | tee -a "$LOG"
  else
    echo "[pdf] batch $bid: FAILED" | tee -a "$LOG"
    FAILED_BATCHES+=("$bid")
  fi
done

if (( ${#FAILED_BATCHES[@]} > 0 )); then
  echo "[pdf] FAILED batches: ${FAILED_BATCHES[*]}" | tee -a "$LOG"
  echo "[pdf] Try: $PY -m pip install weasyprint   (or apt install texlive-xetex and set PDF_ENGINE=xelatex)" | tee -a "$LOG"
  exit 1
fi
