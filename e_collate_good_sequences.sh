#!/usr/bin/env bash
# =============================================================================
# IV_collate_good_sequences.sh
#
# Collates "good" amplicon consensus sequences across the three pipeline
# output folders under III_PipelineOutput/ and writes a single 4th folder
# (z_Final_Good_Sequences/) containing:
#
#   01_FASTA/
#     all_good_sequences.fa          - one combined multi-FASTA
#     per_amplicon/<ref>.fa          - one FASTA per chosen consensus
#     good_sequences_manifest.tsv    - provenance + quality metrics
#   02_BLAST/
#     blast_pairwise.aln             - blastn -outfmt 0 (human-readable)
#     blast_tabular.tsv              - blastn -outfmt 6
#   03_Report/
#     Report.md                      - markdown summary
#     Report.pdf                     - pandoc -> weasyprint PDF
#
# Selection policy (deterministic):
#   1. Every Tri_Consensus FASTA (status=merged, BLAST-confirmed) is kept.
#      These are the gold standard: both batches agree, reference verified.
#   2. For amplicons NOT covered by Tri_Consensus, the batch-A or batch-B
#      single-batch consensus is included only if its N percentage is
#      below MAX_N_PCT (default 20%).
#
# Usage:
#   ./IV_collate_good_sequences.sh
#
# Override Python:  PYTHON=/path/to/python3.11+
# =============================================================================
set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Resource configuration: auto-detect cores for BLAST and per-ref parallelism.
# ---------------------------------------------------------------------------
CPU=$(nproc 2>/dev/null || echo 4)
MAX_PARALLEL="$CPU"

# wait_for_slot: block until fewer than $1 background jobs are running.
wait_for_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

# -------------------- USER-EDITABLE CONFIG ----------------------------------
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_OUTPUT="$ROOT/III_PipelineOutput"

BATCH_A_ID="HC01463503"
BATCH_B_ID="HC01468262"
TRI_DIR_NAME="Tri_Consensus"
OUTPUT_DIR_NAME="z_Final_Good_Sequences"

# Per-batch consensus directories (relative to PIPELINE_OUTPUT).
BATCH_A_CONSENSUS="$BATCH_A_ID/low_quality_saver_pipeline/04a_Consensus_TraceReview/I_$BATCH_A_ID"
BATCH_B_CONSENSUS="$BATCH_B_ID/low_quality_saver_pipeline/04a_Consensus_TraceReview/I_$BATCH_B_ID"
BATCH_A_SUMMARY="$BATCH_A_ID/low_quality_saver_pipeline/05_Report/consensus_summary.csv"
BATCH_B_SUMMARY="$BATCH_B_ID/low_quality_saver_pipeline/05_Report/consensus_summary.csv"

# BLAST reference (same multi-gene DMP reference used by tri_consensus_pipeline).
BLAST_REFERENCE="$ROOT/II_QualityClassified/PF05078_DMP_transcripts_cdhit.fa"
BLAST_EVALUE="1e-5"
BLAST_WORD_SIZE=9
BLAST_THREADS="$CPU"
BLAST_MIN_PIDENT=70.0
BLAST_MIN_QCOV=20.0

# Quality threshold for accepting a single-batch consensus when no
# tri-consensus exists for that reference accession.
MAX_N_PCT=20.0

# Reference accessions to exclude from the final good sequences regardless of
# source or quality. Match the per-amplicon filename stem (e.g. the value used
# for OUT_PER_DIR/<ref>.fa).
EXCLUDE_REFS=(
  "SMEL4.1_12g004810.1.01"
)

# PDF rendering.
BUILD_PDF=true
PDF_ENGINE="weasyprint"
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
    "$cand" -c 'import sys; sys.exit(0 if sys.version_info>=(3,11) else 1)' 2>/dev/null \
      && { echo "$cand"; return 0; }
  done
  return 1
}
PY="$(_pick_python)" || { echo "ERROR: need Python >= 3.11." >&2; exit 2; }

command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found (run setup.sh)." >&2; exit 2; }
command -v blastn      >/dev/null 2>&1 || { echo "ERROR: blastn not found (run setup.sh)." >&2; exit 2; }

# -------------------- logging -----------------------------------------------
LOG_DIR="${LOG_DIR:-$ROOT/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/IV_collate_good_sequences_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[log] writing to $LOG_FILE"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# -------------------- preflight ---------------------------------------------
[[ -d "$PIPELINE_OUTPUT" ]] || die "III_PipelineOutput not found: $PIPELINE_OUTPUT"
[[ -d "$PIPELINE_OUTPUT/$BATCH_A_ID" ]] || die "batch A dir missing: $PIPELINE_OUTPUT/$BATCH_A_ID"
[[ -d "$PIPELINE_OUTPUT/$BATCH_B_ID" ]] || die "batch B dir missing: $PIPELINE_OUTPUT/$BATCH_B_ID"
[[ -d "$PIPELINE_OUTPUT/$TRI_DIR_NAME" ]] || die "Tri_Consensus dir missing: $PIPELINE_OUTPUT/$TRI_DIR_NAME"
[[ -f "$BLAST_REFERENCE" ]] || die "BLAST reference not found: $BLAST_REFERENCE"

OUT_ROOT="$PIPELINE_OUTPUT/$OUTPUT_DIR_NAME"
OUT_FASTA_DIR="$OUT_ROOT/01_FASTA"
OUT_PER_DIR="$OUT_FASTA_DIR/per_amplicon"
OUT_BLAST_DIR="$OUT_ROOT/02_BLAST"
OUT_REPORT_DIR="$OUT_ROOT/03_Report"
mkdir -p "$OUT_PER_DIR" "$OUT_BLAST_DIR" "$OUT_REPORT_DIR"

COMBINED_FA="$OUT_FASTA_DIR/all_good_sequences.fa"
MANIFEST="$OUT_FASTA_DIR/good_sequences_manifest.tsv"
BLAST_PAIRWISE="$OUT_BLAST_DIR/blast_pairwise.aln"
BLAST_TABULAR="$OUT_BLAST_DIR/blast_tabular.tsv"
REPORT_MD="$OUT_REPORT_DIR/Report.md"
REPORT_PDF="$OUT_REPORT_DIR/Report.pdf"

log "output root  : $OUT_ROOT"
log "batch A      : $BATCH_A_ID"
log "batch B      : $BATCH_B_ID"
log "tri folder   : $TRI_DIR_NAME"
log "max N%       : $MAX_N_PCT"

# -------------------- helpers -----------------------------------------------
# Extract ref= token from the first FASTA header in a file.
extract_ref() {
  awk '
    /^>/ {
      for (i=1; i<=NF; i++) if ($i ~ /^ref=/) { sub(/^ref=/,"",$i); print $i; exit }
      exit
    }
  ' "$1"
}

# Compute N percentage of a FASTA file (ignoring headers).
n_pct() {
  "$PY" - "$1" <<'EOF'
import sys
fa = sys.argv[1]
seq = []
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            continue
        seq.append(line.strip())
s = "".join(seq).upper()
if not s:
    print("0.00")
else:
    n = sum(1 for c in s if c == "N")
    print(f"{100.0*n/len(s):.2f}")
EOF
}

# -------------------- index Tri_Consensus FASTAs ----------------------------
declare -A TRI_BY_REF=()
while IFS= read -r -d '' fa; do
  base="$(basename "$fa")"
  ref="${base%_tri_consensus.fa}"
  TRI_BY_REF["$ref"]="$fa"
done < <(find "$PIPELINE_OUTPUT/$TRI_DIR_NAME" -maxdepth 1 -type f -name '*_tri_consensus.fa' -print0)

log "tri_consensus found: ${#TRI_BY_REF[@]}"

# -------------------- index batch consensus FASTAs --------------------------
declare -A A_BY_REF=()
declare -A B_BY_REF=()

while IFS= read -r -d '' fa; do
  ref="$(extract_ref "$fa")"
  [[ -z "$ref" ]] && continue
  A_BY_REF["$ref"]="$fa"
done < <(find "$PIPELINE_OUTPUT/$BATCH_A_CONSENSUS" -maxdepth 1 -type f -name '*.fa' -print0)

while IFS= read -r -d '' fa; do
  ref="$(extract_ref "$fa")"
  [[ -z "$ref" ]] && continue
  B_BY_REF["$ref"]="$fa"
done < <(find "$PIPELINE_OUTPUT/$BATCH_B_CONSENSUS" -maxdepth 1 -type f -name '*.fa' -print0)

log "batch A consensus: ${#A_BY_REF[@]}"
log "batch B consensus: ${#B_BY_REF[@]}"

# -------------------- collate per policy ------------------------------------
: > "$COMBINED_FA"
printf 'ref_accession\tsource\tinput_file\tlength\tN_pct\tdecision\n' > "$MANIFEST"

# Collect the unique set of ref accessions across all three sources.
declare -A ALL_REFS=()
for r in "${!TRI_BY_REF[@]}"; do ALL_REFS["$r"]=1; done
for r in "${!A_BY_REF[@]}";   do ALL_REFS["$r"]=1; done
for r in "${!B_BY_REF[@]}";   do ALL_REFS["$r"]=1; done

kept=0
rejected=0

# Sort refs for deterministic output.
mapfile -t REF_LIST < <(printf '%s\n' "${!ALL_REFS[@]}" | sort)

declare -A EXCLUDED_REFS=()
for r in "${EXCLUDE_REFS[@]:-}"; do
  [[ -z "$r" ]] && continue
  EXCLUDED_REFS["$r"]=1
done

# Temp dir for per-ref compute results (parallel n_pct + seq_len jobs).
COLLATE_TMP="$OUT_ROOT/.collate_rows"
mkdir -p "$COLLATE_TMP"

# Pass 1: dispatch per-ref compute jobs in parallel.
# Each job writes a one-line result file: <chosen_fa> <source> <npct> <seq_len>
for ref in "${REF_LIST[@]}"; do
  chosen_fa=""
  source_label=""

  if [[ -n "${EXCLUDED_REFS[$ref]:-}" ]]; then
    log "  [skip] $ref  (in EXCLUDE_REFS)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ref" "excluded" "-" "-" "-" "excluded_by_config" \
      >> "$MANIFEST"
    continue
  fi

  if [[ -n "${TRI_BY_REF[$ref]:-}" ]]; then
    chosen_fa="${TRI_BY_REF[$ref]}"
    source_label="tri_consensus"
  elif [[ -n "${A_BY_REF[$ref]:-}" ]]; then
    chosen_fa="${A_BY_REF[$ref]}"
    source_label="batch_a_only"
  elif [[ -n "${B_BY_REF[$ref]:-}" ]]; then
    chosen_fa="${B_BY_REF[$ref]}"
    source_label="batch_b_only"
  fi

  [[ -z "$chosen_fa" ]] && continue

  _ref="$ref"; _cfa="$chosen_fa"; _src="$source_label"
  _row_file="$COLLATE_TMP/${ref}.dat"

  wait_for_slot "$MAX_PARALLEL"
  (
    npct=$(n_pct "$_cfa")
    seq_len=$("$PY" - "$_cfa" <<'EOF'
import sys
fa = sys.argv[1]
seq = []
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            continue
        seq.append(line.strip())
print(len("".join(seq)))
EOF
)
    printf '%s\t%s\t%s\t%s\t%s\n' "$_ref" "$_src" "$_cfa" "$npct" "$seq_len" > "$_row_file"
  ) &
done
wait  # collect all compute jobs

# Pass 2: sequential merge in sorted ref order (preserves COMBINED_FA order).
for ref in "${REF_LIST[@]}"; do
  row_file="$COLLATE_TMP/${ref}.dat"
  [[ -f "$row_file" ]] || continue
  IFS=$'\t' read -r _ref source_label chosen_fa npct seq_len < "$row_file"

  # Tri_consensus always accepted; single-batch must pass N% threshold.
  decision="kept"
  if [[ "$source_label" != "tri_consensus" ]]; then
    above=$(awk -v n="$npct" -v t="$MAX_N_PCT" 'BEGIN{print (n+0 > t+0)?"yes":"no"}')
    if [[ "$above" == "yes" ]]; then
      decision="rejected_high_N"
    fi
  fi

  if [[ "$decision" == "kept" ]]; then
    cp -f "$chosen_fa" "$OUT_PER_DIR/${ref}.fa"
    cat "$chosen_fa" >> "$COMBINED_FA"
    printf '\n' >> "$COMBINED_FA"
    kept=$((kept+1))
    log "  [keep] $ref  src=$source_label  N%=${npct}  len=${seq_len}"
  else
    rejected=$((rejected+1))
    log "  [drop] $ref  src=$source_label  N%=${npct}  len=${seq_len}  ($decision)"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ref" "$source_label" "$(basename "$chosen_fa")" "$seq_len" "$npct" "$decision" \
    >> "$MANIFEST"
done
rm -rf "$COLLATE_TMP"

log "kept=$kept  rejected=$rejected"

if (( kept == 0 )); then
  die "no good sequences collected; aborting"
fi

# -------------------- BLAST -------------------------------------------------
BLAST_DB_DIR="$OUT_BLAST_DIR/db"
mkdir -p "$BLAST_DB_DIR"
BLAST_DB="$BLAST_DB_DIR/ref"
if [[ ! -f "${BLAST_DB}.nhr" ]]; then
  log "[blast] makeblastdb -> $BLAST_DB"
  makeblastdb -in "$BLAST_REFERENCE" -dbtype nucl -out "$BLAST_DB" >/dev/null
fi

log "[blast] pairwise alignment -> $BLAST_PAIRWISE"
blastn \
  -query "$COMBINED_FA" -db "$BLAST_DB" \
  -evalue "$BLAST_EVALUE" -word_size "$BLAST_WORD_SIZE" \
  -num_threads "$BLAST_THREADS" \
  -max_target_seqs 3 \
  -outfmt 0 -out "$BLAST_PAIRWISE"

log "[blast] tabular -> $BLAST_TABULAR"
{
  printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\n'
  blastn \
    -query "$COMBINED_FA" -db "$BLAST_DB" \
    -evalue "$BLAST_EVALUE" -word_size "$BLAST_WORD_SIZE" \
    -num_threads "$BLAST_THREADS" \
    -max_target_seqs 1 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'
} > "$BLAST_TABULAR"

# -------------------- Markdown report ---------------------------------------
log "[report] writing $REPORT_MD"

# Build a per-amplicon BLAST verdict table from the tabular file.
"$PY" - "$MANIFEST" "$BLAST_TABULAR" "$REPORT_MD" \
  "$BLAST_MIN_PIDENT" "$BLAST_MIN_QCOV" \
  "$BATCH_A_ID" "$BATCH_B_ID" "$MAX_N_PCT" "$kept" "$rejected" \
  "$BLAST_PAIRWISE" <<'EOF'
import csv, sys, datetime, os

(manifest_p, blast_p, report_p, min_pid, min_qcov,
 batch_a, batch_b, max_n_pct, kept, rejected, pairwise_p) = sys.argv[1:]

min_pid = float(min_pid)
min_qcov = float(min_qcov)

# Load manifest.
rows = []
with open(manifest_p) as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    for r in rd:
        rows.append(r)

# Load best BLAST hit per query (top row per qseqid in the tabular file).
best = {}
with open(blast_p) as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    for r in rd:
        q = r['qseqid']
        if q not in best:
            best[q] = r

# Map manifest ref -> first query id whose FASTA we built. The combined fasta
# headers come straight from the source FASTAs, so qseqid is the first whitespace
# token of that header. For tri_consensus that token is "<ref>|HC..._tri_consensus";
# for single-batch it is "I_<batch>|<amplicon>". Match by substring of ref.
def find_q_for_ref(ref):
    # 1. Exact prefix match (tri_consensus header begins with the ref accession).
    for q in best:
        if q.startswith(ref + '|') or q == ref:
            return q
    # 2. Fallback: look up the per-amplicon FASTA we copied and grab its header.
    per_fa = os.path.join(
        os.path.dirname(manifest_p), 'per_amplicon', f'{ref}.fa')
    if os.path.exists(per_fa):
        with open(per_fa) as fh:
            for line in fh:
                if line.startswith('>'):
                    tok = line[1:].split()[0]
                    return tok if tok in best else None
    return None

now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

md = []
md.append("# Good Sequences Collation Report")
md.append("")
md.append(f"Generated: {now}")
md.append("")
md.append("## Purpose")
md.append("")
md.append(
    "Collates the best per-amplicon consensus across the three pipeline output "
    "folders (batch A, batch B, and the cross-batch Tri_Consensus) into a single "
    f"4th output folder, BLASTs every selected sequence against the multi-gene DMP "
    "reference panel, and reports per-amplicon verdicts."
)
md.append("")
md.append("## Selection policy")
md.append("")
md.append(
    f"1. **Tri_Consensus** entries (both batches agree, reference-confirmed) "
    f"are always retained.\n"
    f"2. For reference accessions not covered by a tri-consensus, the single-batch "
    f"consensus from `{batch_a}` or `{batch_b}` is retained only if its N "
    f"percentage is at or below **{max_n_pct}%**."
)
md.append("")
md.append("## Counts")
md.append("")
md.append(f"- Kept: **{kept}**")
md.append(f"- Rejected (high N%): **{rejected}**")
md.append(f"- Total amplicons considered: **{len(rows)}**")
md.append("")
md.append("## Per-amplicon manifest with BLAST verdict")
md.append("")
md.append("| Ref accession | Source | Input file | Length | N% | Decision | Best BLAST hit | pident | qcov | evalue | Verdict |")
md.append("|---|---|---|---:|---:|---|---|---:|---:|---|---|")

for r in rows:
    ref = r['ref_accession']
    src = r['source']
    fn = r['input_file']
    length = r['length']
    npct = r['N_pct']
    decision = r['decision']
    qid = find_q_for_ref(ref) if decision == 'kept' else None
    if qid and qid in best:
        h = best[qid]
        best_ref = h['sseqid']
        pid = float(h['pident'])
        qcov = float(h['qcovs'])
        ev = h['evalue']
        if pid >= min_pid and qcov >= min_qcov:
            verdict = 'match' if best_ref == ref else 'mismatch_ref'
        else:
            verdict = 'fail_threshold'
        md.append(
            f"| `{ref}` | {src} | {fn} | {length} | {npct} | {decision} "
            f"| `{best_ref}` | {pid:.2f} | {qcov:.0f} | {ev} | {verdict} |"
        )
    elif decision == 'kept':
        md.append(
            f"| `{ref}` | {src} | {fn} | {length} | {npct} | {decision} "
            f"| NO_HIT | - | - | - | no_hit |"
        )
    else:
        md.append(
            f"| `{ref}` | {src} | {fn} | {length} | {npct} | {decision} "
            f"| (not BLASTed) | - | - | - | - |"
        )

md.append("")
md.append("## BLAST configuration")
md.append("")
md.append(
    f"- Reference: `II_QualityClassified/PF05078_DMP_transcripts_cdhit.fa`\n"
    f"- Min %identity: {min_pid}\n"
    f"- Min query coverage: {min_qcov}\n"
)
md.append("")
md.append("## Outputs")
md.append("")
md.append("- `01_FASTA/all_good_sequences.fa` -- combined multi-FASTA of kept sequences")
md.append("- `01_FASTA/per_amplicon/<ref>.fa` -- one FASTA per kept sequence")
md.append("- `01_FASTA/good_sequences_manifest.tsv` -- provenance and quality metrics")
md.append("- `02_BLAST/blast_pairwise.aln` -- blastn pairwise alignment (`-outfmt 0`)")
md.append("- `02_BLAST/blast_tabular.tsv` -- blastn tabular (`-outfmt 6` + qcovs)")
md.append("- `03_Report/Report.md` and `Report.pdf` -- this report")
md.append("")
md.append("## BLAST pairwise alignment (`-outfmt 0`)")
md.append("")
md.append(
    "Verbatim contents of `02_BLAST/blast_pairwise.aln`. Each query is the "
    "kept consensus listed above; the database is the multi-gene DMP "
    "reference panel."
)
md.append("")
md.append("```text")
with open(pairwise_p) as fh:
    md.append(fh.read().rstrip())
md.append("```")
md.append("")

with open(report_p, 'w') as fh:
    fh.write('\n'.join(md))
print(f"wrote {report_p}")
EOF

# -------------------- Markdown -> PDF ---------------------------------------
if [[ "${BUILD_PDF,,}" != "true" ]]; then
  log "[pdf] BUILD_PDF=false, skipping PDF"
  log "DONE. report (md only): $REPORT_MD"
  exit 0
fi

if ! command -v pandoc >/dev/null 2>&1; then
  log "[pdf] pandoc not installed; skipping PDF (apt install pandoc)"
  exit 0
fi

ENGINE_NAME=""
if [[ "$PDF_ENGINE" == "weasyprint" ]]; then
  if "$PY" -c 'import weasyprint' >/dev/null 2>&1; then
    WRAPPER_DIR="$LOG_DIR/.pdf_engine_bin"
    mkdir -p "$WRAPPER_DIR"
    WRAPPER="$WRAPPER_DIR/weasyprint"
    cat > "$WRAPPER" <<EOF2
#!/usr/bin/env bash
exec "$PY" -m weasyprint "\$@"
EOF2
    chmod +x "$WRAPPER"
    export PATH="$WRAPPER_DIR:$PATH"
    ENGINE_NAME="weasyprint"
    log "[pdf] using weasyprint via: $PY -m weasyprint"
  elif command -v weasyprint >/dev/null 2>&1 && weasyprint --version >/dev/null 2>&1; then
    ENGINE_NAME="weasyprint"
  else
    log "[pdf] weasyprint not usable; install with: $PY -m pip install weasyprint"
    exit 1
  fi
elif [[ -n "$PDF_ENGINE" ]]; then
  ENGINE_NAME="$PDF_ENGINE"
fi

REPORT_CSS="$OUT_REPORT_DIR/.report_style.css"
cat > "$REPORT_CSS" <<EOF
@page { size: A4 landscape; margin: ${PDF_MARGIN}; }
html { font-size: ${PDF_FONT_SIZE}; }
body { font-family: "DejaVu Sans", "Helvetica", sans-serif; line-height: 1.3; }
h1 { font-size: 1.6em; }
h2 { font-size: 1.3em; }
h3 { font-size: 1.1em; }
pre, code { font-family: "DejaVu Sans Mono", monospace; font-size: 0.85em; }
pre { white-space: pre-wrap; overflow-wrap: break-word; }
table { font-size: 0.82em; border-collapse: collapse; margin: 0.5em 0; width: 100%; table-layout: auto; font-variant-numeric: tabular-nums; }
th, td { padding: 4px 6px; border: 1px solid #888; overflow-wrap: anywhere; white-space: normal; vertical-align: top; }
th { background: #d9d9d9; text-align: left; font-weight: bold; }
tr:nth-child(even) td { background: #f4f4f4; }
td code, th code { font-size: inherit; white-space: normal; overflow-wrap: anywhere; word-break: break-all; }
table th:nth-child(1), table td:nth-child(1),
table th:nth-child(7), table td:nth-child(7) { min-width: 1.55in; white-space: nowrap; }
table th:nth-child(1) code, table td:nth-child(1) code,
table th:nth-child(7) code, table td:nth-child(7) code { white-space: nowrap; overflow-wrap: normal; }
table th:nth-child(3), table td:nth-child(3) { min-width: 1.8in; }
table th:nth-child(4),  table td:nth-child(4),
table th:nth-child(5),  table td:nth-child(5),
table th:nth-child(8),  table td:nth-child(8),
table th:nth-child(9),  table td:nth-child(9),
table th:nth-child(10), table td:nth-child(10),
table th:nth-child(11), table td:nth-child(11) { width: 1%; white-space: nowrap; padding-left: 1px; padding-right: 1px; }
table th:nth-child(2),  table td:nth-child(2),
table th:nth-child(6),  table td:nth-child(6) { white-space: nowrap; }
img { max-width: 100%; }
EOF

REPORT_TPL="$OUT_REPORT_DIR/.report_template.html"
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
  "$(basename "$REPORT_MD")" -o "$(basename "$REPORT_PDF")"
  --resource-path="$OUT_REPORT_DIR"
  --metadata pagetitle="Good Sequences Collation Report"
  --template="$(basename "$REPORT_TPL")"
  --css="$(basename "$REPORT_CSS")"
)
[[ -n "$ENGINE_NAME" ]] && PANDOC_ARGS+=(--pdf-engine="$ENGINE_NAME")

log "[pdf] (cd $OUT_REPORT_DIR && pandoc ${PANDOC_ARGS[*]})"
( cd "$OUT_REPORT_DIR" && pandoc "${PANDOC_ARGS[@]}" )

log "DONE."
log "outputs under: $OUT_ROOT"
log "  combined FASTA: $COMBINED_FA"
log "  BLAST pairwise: $BLAST_PAIRWISE"
log "  BLAST tabular:  $BLAST_TABULAR"
log "  Report (md):    $REPORT_MD"
log "  Report (pdf):   $REPORT_PDF"
