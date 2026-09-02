#!/usr/bin/env bash
# =============================================================================
# a_merge_consensus_pipeline.sh
#
# Standalone F+R consensus-merging pipeline.
#
# Walks <io.input_root>/<tier>/<sample>/, pairs each *_forward.ab1 with its
# matching *_reverse.ab1, and writes a single merged consensus FASTA per
# amplicon to <io.output_root>/<tier>/<sample>/, preserving the by-quality
# organization of the input tree.
#
# Usage:
#   ./a_merge_consensus_pipeline.sh [a_merge_consensus_config.toml]
# =============================================================================
set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Resource configuration: auto-detect cores; pull override from TOML after
# config is loaded (PROJECT_THREADS). Cap parallel jobs so total CPU load
# stays bounded: one Python worker per core is fine for I/O-bound consensus.
# ---------------------------------------------------------------------------
CPU=$(nproc 2>/dev/null || echo 4)
MAX_PARALLEL="$CPU"

CONFIG="${1:-a_merge_consensus_config.toml}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
COMMON_SCRIPTS="$ROOT/scripts/common"
MERGE_SCRIPTS="$ROOT/scripts/merge_consensus"
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
PYTHON="$(_pick_python)" || { echo "ERROR: need Python >= 3.11." >&2; exit 2; }

# Resolve a relative config path against the script's directory.
[[ "$CONFIG" = /* ]] || CONFIG="$ROOT/$CONFIG"

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }
[[ -f "$COMMON_SCRIPTS/01_read_toml.py" ]] || { echo "ERROR: missing scripts/common/01_read_toml.py" >&2; exit 2; }
[[ -f "$MERGE_SCRIPTS/01_pairwise_consensus.py" ]] || { echo "ERROR: missing scripts/merge_consensus/01_pairwise_consensus.py" >&2; exit 2; }

# wait_for_slot: block until fewer than $1 background jobs are running.
wait_for_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

# ------------------ logging ------------------------------------------------
# Mirror all stdout/stderr to a timestamped file under logs/.
LOG_DIR="${LOG_DIR:-$ROOT/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/merge_consensus_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[log] writing to $LOG_FILE"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ------------------ load config -> shell vars -------------------------------
# 01_read_toml.py --export emits KEY=value lines (uppercased, dots -> underscores)
eval "$("$PYTHON" "$COMMON_SCRIPTS/01_read_toml.py" "$CONFIG" --export)"

# Override MAX_PARALLEL with TOML value if set (PROJECT_THREADS key from [project]).
if [[ -n "${PROJECT_THREADS:-}" && "$PROJECT_THREADS" -gt 0 ]] 2>/dev/null; then
  MAX_PARALLEL="$PROJECT_THREADS"
fi

# [trim].enabled = false -> force q_cutoff to 0 (disables trimming in the
# 01_pairwise_consensus.py consumer). Default to disabled if section is missing.
if [[ "${TRIM_ENABLED:-false}" != "true" ]]; then
  TRIM_Q_CUTOFF=0
fi
TRIM_Q_CUTOFF="${TRIM_Q_CUTOFF:-0}"
TRIM_MIN_LEN="${TRIM_MIN_LEN:-50}"

INPUT_ROOT="$ROOT/$IO_INPUT_ROOT"
OUTPUT_ROOT="$ROOT/$IO_OUTPUT_ROOT"
PIPELINE_DIRNAME="$IO_PIPELINE_DIRNAME"
[[ -d "$INPUT_ROOT" ]] || die "input_root not found: $INPUT_ROOT"
mkdir -p "$OUTPUT_ROOT"

# Read the [[batches]] table. Each row has per-batch overrides; defaults come
# from [reads] (filled in by 01_read_toml.py --batches).
BATCHES_TSV="$("$PYTHON" "$COMMON_SCRIPTS/01_read_toml.py" "$CONFIG" --batches)" \
  || die "Failed to parse [[batches]] from $CONFIG"

log "project       : $PROJECT_NAME"
log "input_root    : $INPUT_ROOT"
log "output_root   : $OUTPUT_ROOT"
log "pipeline_dir  : $PIPELINE_DIRNAME"
log "tiers         : $IO_TIERS"
log "min_quality   : $CONSENSUS_MIN_QUALITY"
log "disagreement  : $CONSENSUS_DISAGREEMENT"
log "min_overlap   : $CONSENSUS_MIN_OVERLAP"
if [[ "$TRIM_Q_CUTOFF" -gt 0 ]]; then
  log "end-trim      : Mott Q=$TRIM_Q_CUTOFF  min_len=$TRIM_MIN_LEN"
else
  log "end-trim      : disabled"
fi
log "batches       :"
printf '%s\n' "$BATCHES_TSV" | awk 'NR>1 {printf "                 %s  (%s)\n", $1, $2}'

# Map TSV column names to indexes (defensive against schema drift).
IFS=$'\t' read -r -a BATCH_COLS <<<"$(printf '%s\n' "$BATCHES_TSV" | head -n1)"
col_idx() { local target="$1" i; for i in "${!BATCH_COLS[@]}"; do [[ "${BATCH_COLS[$i]}" == "$target" ]] && { echo "$i"; return; }; done; echo -1; }
IDX_ID="$(col_idx id)";             IDX_SDIR="$(col_idx sample_dir)"
IDX_FGLOB="$(col_idx forward_glob)"
IDX_FTOK="$(col_idx forward_token)"; IDX_RTOK="$(col_idx reverse_token)"

# ------------------ process each batch -------------------------------------
grand_pairs=0
grand_unpaired=0

while IFS=$'\t' read -r -a BR; do
  BATCH_ID="${BR[$IDX_ID]}"
  BATCH_SDIR="${BR[$IDX_SDIR]}"
  BATCH_FGLOB="${BR[$IDX_FGLOB]:-*_forward.ab1}"
  BATCH_FTOK="${BR[$IDX_FTOK]:-_forward.ab1}"
  BATCH_RTOK="${BR[$IDX_RTOK]:-_reverse.ab1}"

  batch_out_root="$OUTPUT_ROOT/$BATCH_ID/$PIPELINE_DIRNAME"
  mkdir -p "$batch_out_root"

  log ""
  log "============================================================"
  log "BATCH $BATCH_ID  (sample_dir=$BATCH_SDIR)"
  log "  out_root      : $batch_out_root"
  log "  forward_glob  : $BATCH_FGLOB"
  log "  forward_token : $BATCH_FTOK   reverse_token: $BATCH_RTOK"
  log "============================================================"

  batch_pairs=0
  batch_unpaired=0

  for tier in $IO_TIERS; do
    tier_in="$INPUT_ROOT/$tier"
    if [[ ! -d "$tier_in" ]]; then
      log "tier '$tier' not present under input_root: skipping"
      continue
    fi
    # Restrict to this batch's sample subdir only; cross-batch sample mixing
    # is what the [[batches]] split is designed to prevent.
    sdir="$tier_in/$BATCH_SDIR"
    if [[ ! -d "$sdir" ]]; then
      log "  tier '$tier' has no sample '$BATCH_SDIR': skipping"
      continue
    fi

    tier_out="$batch_out_root/$tier"
    mkdir -p "$tier_out"
    summary="$tier_out/$OUTPUT_SUMMARY_FILENAME"
    printf 'sample\tamplicon\tforward\treverse\tconsensus_fa\tstatus\n' > "$summary"

    log "----- tier: $tier  sample: $BATCH_SDIR -----"

    sample="$BATCH_SDIR"
    sample_out="$tier_out/$sample"
    mkdir -p "$sample_out"

    # find forward reads in this sample
    while IFS= read -r -d '' fwd; do
      fwd_name="$(basename "$fwd")"
      # derive matching reverse filename
      rev_name="${fwd_name/${BATCH_FTOK}/${BATCH_RTOK}}"
      rev="$sdir/$rev_name"
      # amplicon stem = forward filename without the forward token
      stem="${fwd_name%${BATCH_FTOK}}"

      out_fa="$sample_out/${stem}${OUTPUT_CONSENSUS_SUFFIX}"

      if [[ ! -f "$rev" ]]; then
        log "  [skip] $sample/$stem : no matching reverse ($rev_name)"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$sample" "$stem" "$fwd_name" "MISSING" "-" "no_reverse_mate" >> "$summary"
        batch_unpaired=$((batch_unpaired + 1))
        continue
      fi

      report_arg=()
      if [[ "$OUTPUT_WRITE_PER_PAIR_REPORT" == "true" ]]; then
        report_arg=(--report "$sample_out/${stem}${OUTPUT_REPORT_SUFFIX}")
      fi

      # Capture TSV row before backgrounding so summary is written from the
      # current (foreground) shell, avoiding concurrent write races.
      tsv_row="$(printf '%s\t%s\t%s\t%s\t%s\t%s' \
        "$sample" "$stem" "$fwd_name" "$rev_name" \
        "$tier/$sample/$(basename "$out_fa")" "merged")"

      wait_for_slot "$MAX_PARALLEL"
      {
        "$PYTHON" "$MERGE_SCRIPTS/01_pairwise_consensus.py" \
          --forward-ab1   "$fwd" \
          --reverse-ab1   "$rev" \
          --out-fa        "$out_fa" \
          --seq-id        "${sample}__${stem}" \
          --min-q         "$CONSENSUS_MIN_QUALITY" \
          --disagreement  "$CONSENSUS_DISAGREEMENT" \
          --min-overlap   "$CONSENSUS_MIN_OVERLAP" \
          --match-score    "$ALIGN_MATCH_SCORE" \
          --mismatch-score "$ALIGN_MISMATCH_SCORE" \
          --gap-open       "$ALIGN_GAP_OPEN" \
          --gap-extend     "$ALIGN_GAP_EXTEND" \
          --trim-q         "${TRIM_Q_CUTOFF:-0}" \
          --trim-min-len   "${TRIM_MIN_LEN:-50}" \
          "${report_arg[@]}"
      } &

      printf '%s\n' "$tsv_row" >> "$summary"
      batch_pairs=$((batch_pairs + 1))

    done < <(find "$sdir" -maxdepth 1 -type f -name "$BATCH_FGLOB" -print0)

    # Wait for all background consensus jobs for this tier/sample to finish
    # before inspecting the output directory or writing the tier summary.
    wait

    # remove empty sample_out (no forward reads in this sample)
    if [[ -z "$(ls -A "$sample_out" 2>/dev/null)" ]]; then
      rmdir "$sample_out" 2>/dev/null || true
    fi

    log "tier '$tier' summary -> $summary"
  done

  log "batch '$BATCH_ID' totals: merged_pairs=$batch_pairs  unpaired=$batch_unpaired"
  grand_pairs=$((grand_pairs + batch_pairs))
  grand_unpaired=$((grand_unpaired + batch_unpaired))

done < <(printf '%s\n' "$BATCHES_TSV" | tail -n +2)

log ""
log "DONE. merged_pairs=$grand_pairs  unpaired_forwards=$grand_unpaired"
log "outputs under: $OUTPUT_ROOT"
