#!/usr/bin/env bash
# =============================================================================
# b_low_quality_saver_pipeline.sh - main driver for the PCR_Sequence
# reference-guided basecalling pipeline (designed to recover usable consensus
# sequence from medium/low-quality Sanger reads).
#
# Reads its TOML config, then runs the steps enabled in [run] using the
# method(s) selected in [method].choice. Each step writes to its own
# subfolder under [project].output_dir.
# =============================================================================
set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Resource configuration: auto-detect cores; TOML [project].threads may
# override the per-tool thread count, but MAX_PARALLEL is always nproc-based.
# ---------------------------------------------------------------------------
CPU=$(nproc 2>/dev/null || echo 4)
MAX_PARALLEL="$CPU"

# wait_for_slot: block until fewer than $1 background jobs are running.
wait_for_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

CONFIG="${1:-b_low_quality_saver_pipeline_config.toml}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$ROOT/scripts"

# -------------------- pick a Python 3.11+ interpreter -----------------------
# The active conda env (e.g. pcr_sanger) may ship an old Python (3.5) that
# cannot parse f-strings or import tomllib. Find a usable interpreter.
_pick_python() {
  local cand
  for cand in "${PYTHON:-}" \
              "$HOME/miniconda3/bin/python3" "$HOME/anaconda3/bin/python3" \
              python3.13 python3.12 python3.11 python3 /usr/bin/python3; do
    [[ -z "$cand" ]] && continue
    command -v "$cand" >/dev/null 2>&1 || continue
    "$cand" -c 'import sys, tomllib, Bio' >/dev/null 2>&1 || continue
    "$cand" -c 'import sys; sys.exit(0 if sys.version_info>=(3,11) else 1)' 2>/dev/null \
      && { echo "$cand"; return 0; }
  done
  return 1
}
PYTHON="$(_pick_python)" || {
  echo "ERROR: need Python >= 3.11 with tomllib + Biopython. None found." >&2
  echo "  Set PYTHON=/path/to/python3.11+ (with biopython) and re-run." >&2
  exit 2
}
echo "[python] using $PYTHON ($("$PYTHON" --version 2>&1))"

# Resolve a relative config path against the script's directory.
[[ "$CONFIG" = /* ]] || CONFIG="$ROOT/$CONFIG"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config not found: $CONFIG" >&2
  exit 2
fi

# -------------------- logging -----------------------------------------------
# Mirror all stdout/stderr to a timestamped file under logs/.
LOG_DIR="${LOG_DIR:-$ROOT/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/low_quality_saver_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[log] writing to $LOG_FILE"

# -------------------- load config -> shell vars -----------------------------
# Produces variables like PROJECT_NAME, METHOD_CHOICE, RUN_MODE, BLAST_EVALUE...
eval "$("$PYTHON" "$SCRIPTS/read_toml.py" "$CONFIG" --export)"

# If TOML did not export PROJECT_THREADS (or it is zero), fall back to nproc.
# Per-tool thread counts passed to Python steps use PROJECT_THREADS; parallel
# batch-level dispatch uses MAX_PARALLEL (always = CPU / tool_threads, min 1).
if [[ -z "${PROJECT_THREADS:-}" || "${PROJECT_THREADS:-0}" -eq 0 ]] 2>/dev/null; then
  PROJECT_THREADS="$CPU"
fi
MAX_PARALLEL=$(( CPU / PROJECT_THREADS < 1 ? 1 : CPU / PROJECT_THREADS ))

# [trim].enabled = false -> force q_cutoff to 0 (which disables trimming in
# every consumer script). Default to disabled if the section is missing.
if [[ "${TRIM_ENABLED:-false}" != "true" ]]; then
  TRIM_Q_CUTOFF=0
fi
TRIM_Q_CUTOFF="${TRIM_Q_CUTOFF:-0}"
TRIM_MIN_LEN="${TRIM_MIN_LEN:-50}"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
# WSL DrvFs (/mnt/c) sometimes makes `mkdir -p` fail with "File exists" even
# when the dir is already present. Use this wrapper everywhere instead.
ensure_dir() { [[ -d "$1" ]] || mkdir -p "$1" 2>/dev/null || [[ -d "$1" ]] || die "cannot create directory: $1"; }

# -------------------- step selection ----------------------------------------
# If [run].operations is provided in the TOML, it is the authoritative list:
# each step runs iff its name appears in the array (and any commented-out
# entry in the TOML is silently dropped, since TOML never sees it). When the
# array is absent we fall back to the legacy [run].mode preset + per-step
# booleans for backward compatibility with older configs.
if [[ -n "${RUN_OPERATIONS:-}" ]]; then
  in_ops() { case " $RUN_OPERATIONS " in *" $1 "*) return 0;; *) return 1;; esac; }
  # `call` is the generic operation name (it dispatches to whichever methods
  # are in [method].choice). For self-documenting configs the method names
  # themselves -- trace_review, consensus_assembly, tracy_decompose -- may
  # appear in operations as aliases for `call` and have the same effect.
  _call_in_ops() {
    in_ops call && return 0
    local m
    for m in $METHOD_CHOICE; do in_ops "$m" && return 0; done
    return 1
  }
  RUN_CLASSIFY=$(       in_ops classify        && echo true || echo false)
  RUN_BLAST=$(          in_ops blast           && echo true || echo false)
  RUN_CALL=$(           _call_in_ops           && echo true || echo false)
  RUN_CONSENSUS=$(      in_ops consensus       && echo true || echo false)
  RUN_CONSENSUS_BLAST=$(in_ops consensus_blast && echo true || echo false)
  RUN_REPORT=$(         in_ops report          && echo true || echo false)
  RUN_BUILD_REPORT=$(   in_ops build_report    && echo true || echo false)
  log "step selection : from [run].operations -> $RUN_OPERATIONS"
else
  case "${RUN_MODE:-custom}" in
    qc)     RUN_CLASSIFY=true;  RUN_BLAST=false; RUN_CALL=false; RUN_CONSENSUS=false; RUN_CONSENSUS_BLAST=false; RUN_REPORT=true;  RUN_BUILD_REPORT=false ;;
    assign) RUN_CLASSIFY=true;  RUN_BLAST=true;  RUN_CALL=false; RUN_CONSENSUS=false; RUN_CONSENSUS_BLAST=false; RUN_REPORT=true;  RUN_BUILD_REPORT=false ;;
    call)   RUN_CLASSIFY=true;  RUN_BLAST=true;  RUN_CALL=true;  RUN_CONSENSUS=false; RUN_CONSENSUS_BLAST=false; RUN_REPORT=true;  RUN_BUILD_REPORT=false ;;
    full)   RUN_CLASSIFY=true;  RUN_BLAST=true;  RUN_CALL=true;  RUN_CONSENSUS=true;  RUN_CONSENSUS_BLAST=true;  RUN_REPORT=true;  RUN_BUILD_REPORT=true ;;
    custom) ;;  # honour the booleans loaded from TOML
    *) die "Unknown run.mode: $RUN_MODE (use qc|assign|call|full|custom)" ;;
  esac
  RUN_CONSENSUS_BLAST="${RUN_CONSENSUS_BLAST:-false}"
  RUN_BUILD_REPORT="${RUN_BUILD_REPORT:-false}"
fi

# -------------------- resolve common paths ---------------------------------
SAMPLE_ROOT="$(cd -- "$ROOT/$PROJECT_SAMPLE_ROOT" && pwd)"
OUTPUT_ROOT="$ROOT/$PROJECT_OUTPUT_ROOT"
PIPELINE_DIRNAME="$PROJECT_PIPELINE_DIRNAME"
REFERENCE="$ROOT/$PROJECT_REFERENCE"
[[ -f "$REFERENCE" ]] || die "Reference FASTA not found: $REFERENCE"

ensure_dir "$OUTPUT_ROOT"

# Read the [[batches]] table. Each row carries the per-batch overrides; defaults
# are pre-filled from [reads] by read_toml.py so every column is populated.
BATCHES_TSV="$("$PYTHON" "$SCRIPTS/read_toml.py" "$CONFIG" --batches)" \
  || die "Failed to parse [[batches]] from $CONFIG"

log "project        : $PROJECT_NAME"
log "methods        : $METHOD_CHOICE"
log "sample_root    : $SAMPLE_ROOT"
log "reference      : $REFERENCE"
log "output_root    : $OUTPUT_ROOT"
log "pipeline_dir   : $PIPELINE_DIRNAME"
log "tiers (filter) : ${PROJECT_TIERS:-<all>}"
log "steps          : classify=$RUN_CLASSIFY blast=$RUN_BLAST call=$RUN_CALL consensus=$RUN_CONSENSUS consensus_blast=$RUN_CONSENSUS_BLAST report=$RUN_REPORT build_report=$RUN_BUILD_REPORT"
log "batches        :"
printf '%s\n' "$BATCHES_TSV" | awk 'NR>1 {printf "                   %s  (%s)\n", $1, $2}'

# Track all merged-consensus output paths across batches (for any future
# cross-batch aggregation). Currently each batch self-contains its outputs.
ALL_CALL_DIRS=()
ALL_CONSENSUS_DIRS=()

# =============================================================================
# Per-batch loop. Each iteration sets BATCH_* env vars, computes batch-local
# OUT_ROOT under III_PipelineOutput/<batch_id>/<pipeline_dirname>/, then runs
# every enabled step inside this batch's scope. Outputs from different batches
# never share a folder.
# =============================================================================
# Read TSV header to map column names to indexes (defensive against schema drift).
IFS=$'\t' read -r -a BATCH_COLS <<<"$(printf '%s\n' "$BATCHES_TSV" | head -n1)"
col_idx() { local target="$1" i; for i in "${!BATCH_COLS[@]}"; do [[ "${BATCH_COLS[$i]}" == "$target" ]] && { echo "$i"; return; }; done; echo -1; }
IDX_ID="$(col_idx id)";              IDX_SDIR="$(col_idx sample_dir)"
IDX_FGLOB="$(col_idx forward_glob)"; IDX_RGLOB="$(col_idx reverse_glob)"
IDX_FSFX="$(col_idx forward_suffix)";IDX_RSFX="$(col_idx reverse_suffix)"

while IFS=$'\t' read -r -a BR; do
  BATCH_ID="${BR[$IDX_ID]}"
  BATCH_SDIR="${BR[$IDX_SDIR]}"
  BATCH_FGLOB="${BR[$IDX_FGLOB]:-*_forward.ab1}"
  BATCH_RGLOB="${BR[$IDX_RGLOB]:-*_reverse.ab1}"
  BATCH_FSUFFIX="${BR[$IDX_FSFX]:-_forward}"
  BATCH_RSUFFIX="${BR[$IDX_RSFX]:-_reverse}"

  BATCH_SAMPLE_PATH="$SAMPLE_ROOT/$BATCH_SDIR"
  [[ -d "$BATCH_SAMPLE_PATH" ]] || { log "  batch '$BATCH_ID' sample_dir not found: $BATCH_SAMPLE_PATH -- skipping"; continue; }

  OUT_ROOT="$OUTPUT_ROOT/$BATCH_ID/$PIPELINE_DIRNAME"
  ensure_dir "$OUT_ROOT"

  CLASSIFY_DIR="$OUT_ROOT/$OUTPUT_CLASSIFY"
  BLAST_DIR="$OUT_ROOT/$OUTPUT_BLAST"
  REPORT_DIR="$OUT_ROOT/$OUTPUT_REPORT"
  # Per-method merged-consensus dirs are resolved inside the method loop below.

  log ""
  log "============================================================"
  log "BATCH $BATCH_ID  (sample_dir=$BATCH_SDIR)"
  log "  out_root     : $OUT_ROOT"
  log "  forward_glob : $BATCH_FGLOB   forward_suffix=$BATCH_FSUFFIX"
  log "  reverse_glob : $BATCH_RGLOB   reverse_suffix=$BATCH_RSUFFIX"
  log "============================================================"

  # -------------------- helpers (batch-scoped) -----------------------------
  sample_dirs() {
    # For per-batch processing, only the batch's own sample folder is visible
    # to downstream steps; this keeps outputs scoped to one batch.
    if compgen -G "$BATCH_SAMPLE_PATH/*.ab1" > /dev/null; then
      echo "$BATCH_SAMPLE_PATH"
    fi
  }

  # =============================================================================
  # STEP 1: classify by Phred quality
  # =============================================================================
  if [[ "$RUN_CLASSIFY" == "true" ]]; then
    log "[1/6] classify -> $CLASSIFY_DIR"
    ensure_dir "$CLASSIFY_DIR"
    CLASSIFY_OUT_DIR="$CLASSIFY_DIR" \
    CLASSIFY_HIGH="$CLASSIFY_HIGH_THRESHOLD" \
    CLASSIFY_MED="$CLASSIFY_MEDIUM_THRESHOLD" \
      "$PYTHON" "$SCRIPTS/step_classify.py" \
        --sample-root    "$SAMPLE_ROOT" \
        --include        "$BATCH_SDIR" \
        --out-dir        "$CLASSIFY_DIR" \
        --high           "$CLASSIFY_HIGH_THRESHOLD" \
        --medium         "$CLASSIFY_MEDIUM_THRESHOLD" \
        --forward-suffix="$BATCH_FSUFFIX" \
        --reverse-suffix="$BATCH_RSUFFIX" \
        --read-glob="$BATCH_FGLOB" \
        --read-glob="$BATCH_RGLOB"
  fi

  # =============================================================================
  # STEP 2: BLAST assignment of each forward AND reverse read to its best DMP
  # transcript. blastn searches both strands by default, so reverse reads can
  # be blasted directly without reverse-complementing the FASTA.
  # =============================================================================
  if [[ "$RUN_BLAST" == "true" ]]; then
    log "[2/6] blast -> $BLAST_DIR"
    ensure_dir "$BLAST_DIR"
    command -v makeblastdb >/dev/null || die "makeblastdb not found (run setup.sh)"
    command -v blastn      >/dev/null || die "blastn not found (run setup.sh)"

    DB_DIR="$BLAST_DIR/db"
    ensure_dir "$DB_DIR"
    if [[ ! -f "$DB_DIR/ref.nhr" ]]; then
      makeblastdb -in "$REFERENCE" -dbtype nucl -out "$DB_DIR/ref" >/dev/null
    fi

    ASSIGN_TSV="$BLAST_DIR/assignments.tsv"
    printf 'sample\tread\tbest_ref\tpident\tqcov\tevalue\n' > "$ASSIGN_TSV"

    # Temp directory: each background job writes one row file to avoid
    # concurrent writes to the shared assignments TSV.
    BLAST_TMP="$BLAST_DIR/.blast_rows"
    ensure_dir "$BLAST_TMP"

    # _blast_one_read <sample> <stem> <ab1> <fa> <db_dir>
    # Runs ab1_to_fasta + blastn for a single read and writes one TSV row.
    _blast_one_read() {
      local sample="$1" stem="$2" ab1="$3" fa="$4" db_dir="$5"
      local row_file="$BLAST_TMP/${sample}_${stem}.row"
      "$PYTHON" "$SCRIPTS/ab1_to_fasta.py" --ab1 "$ab1" --out "$fa" --seq-id "$stem" \
        --trim-q "${TRIM_Q_CUTOFF:-0}" --trim-min-len "${TRIM_MIN_LEN:-50}"
      local hit
      hit=$(blastn -query "$fa" -db "$db_dir/ref" \
              -evalue "$BLAST_EVALUE" -word_size "$BLAST_WORD_SIZE" \
              -outfmt '6 sseqid pident qcovs evalue' \
              -max_target_seqs 1 -num_threads 1 | head -n1 || true)
      if [[ -n "$hit" ]]; then
        local sseq pid qcv ev pass
        sseq=$(awk '{print $1}' <<<"$hit")
        pid=$(awk  '{print $2}' <<<"$hit")
        qcv=$(awk  '{print $3}' <<<"$hit")
        ev=$(awk   '{print $4}' <<<"$hit")
        pass=$(awk -v p="$pid" -v c="$qcv" -v mp="$BLAST_MIN_PIDENT" -v mc="$BLAST_MIN_QCOV" \
               'BEGIN{print (p+0>=mp+0 && c+0>=mc+0)?"yes":"no"}')
        if [[ "$pass" == "yes" ]]; then
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$stem" "$sseq" "$pid" "$qcv" "$ev" > "$row_file"
        else
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$stem" "FAIL_THRESHOLD" "$pid" "$qcv" "$ev" > "$row_file"
        fi
      else
        printf '%s\t%s\tNO_HIT\t-\t-\t-\n' "$sample" "$stem" > "$row_file"
      fi
    }
    export -f _blast_one_read

    while read -r sdir; do
      sample="$(basename "$sdir")"
      while read -r ab1; do
        stem="$(basename "$ab1" .ab1)"
        fa="$BLAST_DIR/$sample/$stem.fa"
        ensure_dir "$(dirname "$fa")"
        # Dispatch each read in background; cap concurrent jobs to MAX_PARALLEL.
        # blastn uses -num_threads 1 per job; concurrency comes from MAX_PARALLEL.
        wait_for_slot "$MAX_PARALLEL"
        _blast_one_read "$sample" "$stem" "$ab1" "$fa" "$DB_DIR" &
      done < <(find "$sdir" -maxdepth 1 \( -name "$BATCH_FGLOB" -o -name "$BATCH_RGLOB" \) | sort)
    done < <(sample_dirs)
    wait  # collect all BLAST jobs before merging

    # Merge per-read row files into the shared TSV (deterministic order).
    while IFS= read -r row_file; do
      cat "$row_file" >> "$ASSIGN_TSV"
    done < <(find "$BLAST_TMP" -maxdepth 1 -name '*.row' | sort)
    rm -rf "$BLAST_TMP"

    # Pair-aware rescue: if one strand of a forward/reverse pair has a valid
    # ref assignment but its partner is NO_HIT / FAIL_THRESHOLD, propagate the
    # partner's assignment so the read is not silently dropped downstream.
    FWD_SUFFIX="$BATCH_FSUFFIX" REV_SUFFIX="$BATCH_RSUFFIX" \
    "$PYTHON" - "$ASSIGN_TSV" <<'PY'
import os, sys, re
from pathlib import Path

tsv = Path(sys.argv[1])
fwd_sfx = os.environ.get("FWD_SUFFIX", "_forward")
rev_sfx = os.environ.get("REV_SUFFIX", "_reverse")

lines = tsv.read_text().splitlines()
header = lines[0].split("\t")
rows = [l.split("\t") for l in lines[1:] if l.strip()]

BAD = {"NO_HIT", "FAIL_THRESHOLD", "-", ""}
strip_re = re.compile(rf"(?:{re.escape(fwd_sfx)}|{re.escape(rev_sfx)})$")

pairs = {}
for i, r in enumerate(rows):
    sample, read = r[0], r[1]
    if read.endswith(fwd_sfx):
        strand = "F"
    elif read.endswith(rev_sfx):
        strand = "R"
    else:
        continue
    base = strip_re.sub("", read)
    pairs.setdefault((sample, base), {})[strand] = i

rescued = 0
for (sample, base), strands in pairs.items():
    fi, ri = strands.get("F"), strands.get("R")
    if fi is None or ri is None:
        continue
    f, r = rows[fi], rows[ri]
    f_bad = f[2] in BAD
    r_bad = r[2] in BAD
    if f_bad and not r_bad:
        f[2] = r[2]  # propagate ref only; keep this strand's own pident/qcov/evalue marks
        rescued += 1
    elif r_bad and not f_bad:
        r[2] = f[2]
        rescued += 1

with tsv.open("w") as fh:
    fh.write("\t".join(header) + "\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

print(f"[blast] pair-aware rescue: propagated ref to {rescued} partner read(s)")
PY

    # Tier filter: drop assignments whose read tier (per quality_report.csv) is
    # not in PROJECT_TIERS. Only applied if both the report and tier list exist.
    QREPORT="$CLASSIFY_DIR/quality_report.csv"
    if [[ -n "${PROJECT_TIERS:-}" && -f "$QREPORT" ]]; then
      log "      tier filter: keeping reads in tiers [$PROJECT_TIERS]"
      PROJECT_TIERS="$PROJECT_TIERS" QREPORT="$QREPORT" \
        "$PYTHON" - "$ASSIGN_TSV" <<'PY'
import csv, os, sys
from pathlib import Path

tsv = Path(sys.argv[1])
qreport = Path(os.environ["QREPORT"])
allowed = set(os.environ["PROJECT_TIERS"].split())

tier_of: dict[tuple[str, str], str] = {}
with qreport.open(newline="") as fh:
    for row in csv.DictReader(fh):
        tier_of[(row["sample"], row["read"])] = row["tier"]

lines = tsv.read_text().splitlines()
header = lines[0]
kept, dropped = [], 0
for line in lines[1:]:
    if not line.strip():
        continue
    parts = line.split("\t")
    sample, read = parts[0], parts[1]
    tier = tier_of.get((sample, read))
    if tier is None or tier in allowed:
        kept.append(line)
    else:
        dropped += 1

with tsv.open("w") as fh:
    fh.write(header + "\n")
    for line in kept:
        fh.write(line + "\n")

print(f"[blast] tier filter: kept {len(kept)} reads, dropped {dropped}")
PY
    fi

    log "      -> $ASSIGN_TSV"
  fi

  # =============================================================================
  # STEP 2b: synthetic assignments (force_assign)
  # When BLAST is skipped (e.g., target paralog is already known) and the TOML
  # provides [blast].force_assign = "<ref_id>", synthesize an assignments.tsv
  # pinning every read in this batch to that reference. The downstream call /
  # consensus steps then run as if BLAST had perfectly classified everything.
  # The reference id must match a header in [project].reference verbatim
  # (sseqid = first whitespace-separated token of the FASTA header).
  # No-op when a real BLAST run wrote $ASSIGN_TSV moments earlier.
  # =============================================================================
  ASSIGN_TSV="$BLAST_DIR/assignments.tsv"
  if [[ "$RUN_BLAST" != "true" && -n "${BLAST_FORCE_ASSIGN:-}" && ! -f "$ASSIGN_TSV" ]]; then
    log "[2b/6] force_assign -> $BLAST_FORCE_ASSIGN (synthesizing $ASSIGN_TSV)"
    ensure_dir "$BLAST_DIR"
    printf 'sample\tread\tbest_ref\tpident\tqcov\tevalue\n' > "$ASSIGN_TSV"
    sample="$(basename "$BATCH_SAMPLE_PATH")"
    n_rows=0
    while IFS= read -r ab1; do
      stem="$(basename "$ab1" .ab1)"
      printf '%s\t%s\t%s\t100.000\t100\t0.0\n' "$sample" "$stem" "$BLAST_FORCE_ASSIGN" >> "$ASSIGN_TSV"
      n_rows=$((n_rows + 1))
    done < <(find "$BATCH_SAMPLE_PATH" -maxdepth 1 \( -name "$BATCH_FGLOB" -o -name "$BATCH_RGLOB" \) | sort)
    log "      -> wrote $n_rows synthetic row(s)"
  fi

# =============================================================================
# STEPS 3-4: loop over each method in METHOD_CHOICE (per-batch)
# =============================================================================
CALL_DIRS=()
CONSENSUS_DIRS=()
CONSENSUS_BLAST_DIRS=()
for method in $METHOD_CHOICE; do
  # Resolve per-method output folders (call dir + merged-consensus dir +
  # consensus-blast dir). tracy_decompose has no consensus / consensus-blast
  # folders because it emits a per-read consensus directly.
  CONSENSUS_DIR=""
  CONSENSUS_BLAST_DIR=""
  case "$method" in
    trace_review)
      CALL_DIR="$OUT_ROOT/$OUTPUT_CALL_DIRS_TRACE_REVIEW"
      CONSENSUS_DIR="$OUT_ROOT/$OUTPUT_CONSENSUS_DIRS_TRACE_REVIEW"
      CONSENSUS_BLAST_DIR="$OUT_ROOT/${OUTPUT_CONSENSUS_BLAST_DIRS_TRACE_REVIEW:-04a_BlastConsensus_TraceReview}"
      ;;
    consensus_assembly)
      CALL_DIR="$OUT_ROOT/$OUTPUT_CALL_DIRS_CONSENSUS_ASSEMBLY"
      CONSENSUS_DIR="$OUT_ROOT/$OUTPUT_CONSENSUS_DIRS_CONSENSUS_ASSEMBLY"
      CONSENSUS_BLAST_DIR="$OUT_ROOT/${OUTPUT_CONSENSUS_BLAST_DIRS_CONSENSUS_ASSEMBLY:-04b_BlastConsensus_ConsensusAssembly}"
      ;;
    tracy_decompose)
      CALL_DIR="$OUT_ROOT/$OUTPUT_CALL_DIRS_TRACY_DECOMPOSE"
      ;;
    *) die "Unknown method.choice: $method" ;;
  esac
  CALL_DIRS+=("$CALL_DIR")
  [[ -n "$CONSENSUS_DIR" ]] && CONSENSUS_DIRS+=("$CONSENSUS_DIR")
  [[ -n "$CONSENSUS_BLAST_DIR" ]] && CONSENSUS_BLAST_DIRS+=("$CONSENSUS_BLAST_DIR")

  # ---------------------------------------------------------------------------
  # STEP 3: call - run the chosen method
  # ---------------------------------------------------------------------------
  if [[ "$RUN_CALL" == "true" ]]; then
    ASSIGN_TSV="$BLAST_DIR/assignments.tsv"
    if [[ ! -f "$ASSIGN_TSV" ]]; then
      log "[3/6] call ($method) -> SKIPPED (no $ASSIGN_TSV; add 'blast' to [run].operations, or run blast once before)"
      continue
    fi
    log "[3/6] call ($method) -> $CALL_DIR"
    ensure_dir "$CALL_DIR"

    case "$method" in
      trace_review)
        "$PYTHON" "$SCRIPTS/step_trace_review.py" \
          --sample-root  "$SAMPLE_ROOT" \
          --include      "$BATCH_SDIR" \
          --reference    "$REFERENCE" \
          --assignments  "$ASSIGN_TSV" \
          --out-dir      "$CALL_DIR" \
          --min-quality          "$TRACE_REVIEW_MIN_QUALITY" \
          --secondary-peak-ratio "$TRACE_REVIEW_SECONDARY_PEAK_RATIO" \
          --allow-ref-overwrite  "$TRACE_REVIEW_ALLOW_REF_OVERWRITE" \
          --require-trace-support "$TRACE_REVIEW_REQUIRE_TRACE_SUPPORT" \
          --mark-uncalled        "$TRACE_REVIEW_MARK_UNCALLED" \
          --tracy-for-indels     "$TRACE_REVIEW_TRACY_FOR_INDELS" \
          --het-indel-trigger    "$TRACE_REVIEW_HET_INDEL_TRIGGER" \
          --max-low-q-frac       "${TRACE_REVIEW_MAX_LOW_Q_FRAC:-0.85}" \
          --threads              "${TRACE_REVIEW_THREADS:-$PROJECT_THREADS}" \
          --trim-q               "${TRIM_Q_CUTOFF:-0}" \
          --trim-min-len         "${TRIM_MIN_LEN:-50}" \
          --forward-suffix="$BATCH_FSUFFIX" \
          --reverse-suffix="$BATCH_RSUFFIX"
        ;;
      consensus_assembly)
        "$PYTHON" "$SCRIPTS/step_consensus_assembly.py" \
          --sample-root  "$SAMPLE_ROOT" \
          --include      "$BATCH_SDIR" \
          --reference    "$REFERENCE" \
          --assignments  "$ASSIGN_TSV" \
          --out-dir      "$CALL_DIR" \
          --aligner      "$CONSENSUS_ASSEMBLY_ALIGNER" \
          --min-quality  "$CONSENSUS_ASSEMBLY_MIN_BASE_QUALITY" \
          --disagreement "$CONSENSUS_ASSEMBLY_DISAGREEMENT" \
          --min-coverage "$CONSENSUS_ASSEMBLY_MIN_COVERAGE" \
          --gap-open     "$CONSENSUS_ASSEMBLY_GAP_OPEN" \
          --gap-extend   "$CONSENSUS_ASSEMBLY_GAP_EXTEND" \
          --threads      "${CONSENSUS_ASSEMBLY_THREADS:-$PROJECT_THREADS}" \
          --trim-q       "${TRIM_Q_CUTOFF:-0}" \
          --trim-min-len "${TRIM_MIN_LEN:-50}" \
          --forward-suffix="$BATCH_FSUFFIX" \
          --reverse-suffix="$BATCH_RSUFFIX"
        ;;
      tracy_decompose)
        command -v tracy >/dev/null || die "tracy not found (run setup.sh)"
        "$PYTHON" "$SCRIPTS/step_tracy_decompose.py" \
          --sample-root  "$SAMPLE_ROOT" \
          --include      "$BATCH_SDIR" \
          --reference    "$REFERENCE" \
          --assignments  "$ASSIGN_TSV" \
          --out-dir      "$CALL_DIR" \
          --mode         "$TRACY_DECOMPOSE_MODE" \
          --pratio       "$TRACY_DECOMPOSE_PRATIO" \
          --trim-left    "$TRACY_DECOMPOSE_TRIM_LEFT" \
          --trim-right   "$TRACY_DECOMPOSE_TRIM_RIGHT" \
          --trim-quality "$TRACY_DECOMPOSE_TRIM_QUALITY" \
          --call-variants "$TRACY_DECOMPOSE_CALL_VARIANTS" \
          --kmer         "$TRACY_DECOMPOSE_KMER" \
          --auto-decompose-q "$TRACY_DECOMPOSE_AUTO_DECOMPOSE_Q" \
          --max-low-q-frac   "${TRACY_DECOMPOSE_MAX_LOW_Q_FRAC:-0.85}" \
          --threads          "${TRACY_DECOMPOSE_THREADS:-$PROJECT_THREADS}" \
          --forward-suffix="$BATCH_FSUFFIX" \
          --reverse-suffix="$BATCH_RSUFFIX"
        ;;
    esac
  fi

  # ---------------------------------------------------------------------------
  # STEP 4: consensus (merge forward + reverse strand calls)
  # ---------------------------------------------------------------------------
  if [[ "$RUN_CONSENSUS" == "true" ]]; then
    if [[ "$method" == "tracy_decompose" ]]; then
      log "[4/6] consensus ($method) -> SKIPPED (tracy_decompose already produces a per-read consensus)"
    elif [[ ! -f "$BLAST_DIR/assignments.tsv" ]]; then
      log "[4/6] consensus ($method) -> SKIPPED (no $BLAST_DIR/assignments.tsv; add 'blast' to [run].operations, or run blast once before)"
    elif [[ ! -d "$CALL_DIR" ]]; then
      log "[4/6] consensus ($method) -> SKIPPED (no call output at $CALL_DIR; add 'call' to [run].operations, or run it once before)"
    else
      log "[4/6] consensus ($method) -> $CONSENSUS_DIR"
      ensure_dir "$CONSENSUS_DIR"
      "$PYTHON" "$SCRIPTS/step_merge_consensus.py" \
        --reference    "$REFERENCE" \
        --call-dir     "$CALL_DIR" \
        --assignments  "$BLAST_DIR/assignments.tsv" \
        --out-dir      "$CONSENSUS_DIR" \
        --min-q        "$CONSENSUS_MIN_STRAND_Q" \
        --policy       "$CONSENSUS_DISAGREEMENT" \
        --min-cov      "$CONSENSUS_MIN_COVERAGE" \
        --forward-suffix="$BATCH_FSUFFIX" \
        --reverse-suffix="$BATCH_RSUFFIX"
    fi
  fi

  # ---------------------------------------------------------------------------
  # STEP 5: consensus_blast (blast each merged-consensus FASTA back against
  # the DMP reference to confirm paralog assignment and emit pairwise alignments)
  # ---------------------------------------------------------------------------
  if [[ "$RUN_CONSENSUS_BLAST" == "true" ]]; then
    if [[ "$method" == "tracy_decompose" ]]; then
      log "[5/6] consensus_blast ($method) -> SKIPPED (no merged-consensus step for tracy_decompose)"
    elif [[ -z "$CONSENSUS_BLAST_DIR" ]]; then
      log "[5/6] consensus_blast ($method) -> SKIPPED (no consensus_blast dir configured)"
    elif [[ ! -d "$CONSENSUS_DIR" ]]; then
      log "[5/6] consensus_blast ($method) -> SKIPPED (consensus dir missing: $CONSENSUS_DIR)"
    else
      log "[5/6] consensus_blast ($method) -> $CONSENSUS_BLAST_DIR"
      ensure_dir "$CONSENSUS_BLAST_DIR"
      command -v makeblastdb >/dev/null || die "makeblastdb not found (run setup.sh)"
      command -v blastn      >/dev/null || die "blastn not found (run setup.sh)"
      THREADS="$PROJECT_THREADS" \
      EVALUE="$BLAST_EVALUE" \
      WORD_SIZE="$BLAST_WORD_SIZE" \
        bash "$SCRIPTS/blast_consensus.sh" "$CONSENSUS_DIR" "$CONSENSUS_BLAST_DIR" "$REFERENCE"
    fi
  fi


done

  # =============================================================================
  # STEP 6: report (per-batch aggregate covering all methods in this batch)
  # =============================================================================
  if [[ "$RUN_REPORT" == "true" ]]; then
    log "[6/6] report -> $REPORT_DIR"
    ensure_dir "$REPORT_DIR"
    "$PYTHON" "$SCRIPTS/step_report.py" \
      --output-root  "$OUT_ROOT" \
      --classify-dir "$CLASSIFY_DIR" \
      --blast-dir    "$BLAST_DIR" \
      --call-dir     "${CALL_DIRS[*]}" \
      --consensus-dir "${CONSENSUS_DIRS[*]:-}" \
      --report-dir   "$REPORT_DIR" \
      --method       "$METHOD_CHOICE" \
      --config       "$CONFIG"
  fi

  # Accumulate batch outputs for any future cross-batch summarization.
  ALL_CALL_DIRS+=("${CALL_DIRS[@]}")
  ALL_CONSENSUS_DIRS+=("${CONSENSUS_DIRS[@]}")

done < <(printf '%s\n' "$BATCHES_TSV" | tail -n +2)

# =============================================================================
# STEP 7: build_report (cross-batch Markdown + PDF report)
# Runs once after the per-batch loop; c_build_report.sh iterates [[batches]] itself.
# =============================================================================
if [[ "$RUN_BUILD_REPORT" == "true" ]]; then
  BUILD_REPORT_SH="$ROOT/c_build_report.sh"
  if [[ ! -x "$BUILD_REPORT_SH" ]]; then
    log "build_report -> SKIPPED (not executable or missing: $BUILD_REPORT_SH)"
  else
    log ""
    log "============================================================"
    log "build_report -> $BUILD_REPORT_SH"
    log "============================================================"
    PYTHON="$PYTHON" LQS_CONFIG="$CONFIG" bash "$BUILD_REPORT_SH"
  fi
fi

log "DONE."
