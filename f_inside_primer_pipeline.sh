#!/usr/bin/env bash
# =============================================================================
# f_inside_primer_pipeline.sh
#
# Dedicated wrapper for I_HC01109784, whose traces use mixed outside/inside
# primer names:
#   PH_2_DMP8-OF + PH_2_DMP8-IR
#   PH_3_DMP8-IF + PH_3_DMP8-OR
#
# The underlying low-quality saver pipeline expects one forward suffix and one
# reverse suffix per batch, so the config exposes these as two virtual batches.
# =============================================================================
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$SCRIPT_DIR/f_inside_primer_pipeline_config.toml}"

# -------------------- logging -----------------------------------------------
# Mirror this wrapper's full stdout/stderr to a timestamped file under logs/.
# b_low_quality_saver_pipeline.sh still writes its own log file (it sets up
# its own tee on exec), so per-batch output ends up in both files; the
# wrapper's log additionally captures overlap_consensus and
# overlap_consensus_blast output that runs after b_*.sh returns.
LOG_DIR="${LOG_DIR:-$SCRIPT_DIR/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/inside_primer_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
printf '[%s] [f_inside_primer_pipeline.sh] writing to %s\n' "$(date +%H:%M:%S)" "$LOG_FILE"
printf '[%s] [f_inside_primer_pipeline.sh] config: %s\n'    "$(date +%H:%M:%S)" "$CONFIG"

if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastn >/dev/null 2>&1; then
  for tool_dir in \
    "${CONDA_PREFIX:-}/bin" \
    "$HOME/miniconda3/envs/pcr_sanger/bin" \
    "$HOME/miniconda3/bin" \
    "$HOME/anaconda3/envs/pcr_sanger/bin" \
    "$HOME/anaconda3/bin"; do
    if [[ -x "$tool_dir/makeblastdb" && -x "$tool_dir/blastn" ]]; then
      export PATH="$tool_dir:$PATH"
      break
    fi
  done
fi

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
PYTHON_BIN="$(_pick_python)" || {
  echo "ERROR: need Python >= 3.11 with tomllib + Biopython." >&2
  exit 2
}

bash "$SCRIPT_DIR/b_low_quality_saver_pipeline.sh" "$CONFIG"

# Read the operations whitelist (if any) so we can gate cross-batch steps the
# same way b_low_quality_saver_pipeline.sh gates its per-batch steps. An empty
# RUN_OPERATIONS means the config did not provide [run].operations — fall back
# to legacy behaviour (always run overlap_consensus when this wrapper is used).
RUN_OPERATIONS="$("$PYTHON_BIN" -c "
import sys, tomllib
with open(sys.argv[1], 'rb') as fh:
    cfg = tomllib.load(fh)
print(' '.join((cfg.get('run') or {}).get('operations') or []))
" "$CONFIG")"

if [[ -z "$RUN_OPERATIONS" ]] || [[ " $RUN_OPERATIONS " == *" overlap_consensus "* ]]; then
  "$PYTHON_BIN" "$SCRIPT_DIR/scripts/overlap_consensus.py" \
    --workspace "$SCRIPT_DIR" \
    --config "$CONFIG"
else
  echo "[overlap_consensus] SKIPPED (not in [run].operations)"
fi

# -----------------------------------------------------------------------------
# overlap_consensus_blast: BLAST the final stitched full-gene FASTA against
# [project].reference to confirm paralog identity end-to-end. Runs only when
# 'overlap_consensus_blast' is in [run].operations. Reads its parameters from
# the [overlap_consensus_blast] section (with sensible fall-backs to [blast]).
# -----------------------------------------------------------------------------
if [[ " $RUN_OPERATIONS " == *" overlap_consensus_blast "* ]]; then
  eval "$("$PYTHON_BIN" -c "
import sys, tomllib
with open(sys.argv[1], 'rb') as fh:
    cfg = tomllib.load(fh)
oc  = cfg.get('overlap_consensus') or {}
ocb = cfg.get('overlap_consensus_blast') or {}
bl  = cfg.get('blast') or {}
proj = cfg.get('project') or {}
out_dir = oc.get('output_dir', '')
seq_id  = oc.get('seq_id', 'overlap_consensus')
ref     = ocb.get('reference', proj.get('reference', ''))
evalue  = ocb.get('evalue', bl.get('evalue', '1e-5'))
ws      = ocb.get('word_size', bl.get('word_size', 9))
which   = ocb.get('which', 'reference_length')  # reference_length | covered_span | both
print(f\"OC_OUT_DIR='{out_dir}'\")
print(f\"OC_SEQ_ID='{seq_id}'\")
print(f\"OCB_REF='{ref}'\")
print(f\"OCB_EVALUE='{evalue}'\")
print(f\"OCB_WORD_SIZE='{ws}'\")
print(f\"OCB_WHICH='{which}'\")
" "$CONFIG")"

  if [[ -z "${OC_OUT_DIR:-}" || -z "${OCB_REF:-}" ]]; then
    echo "[overlap_consensus_blast] SKIPPED (missing overlap_consensus.output_dir or reference)" >&2
  else
    OC_OUT_DIR_ABS="$OC_OUT_DIR"
    [[ "$OC_OUT_DIR_ABS" = /* ]] || OC_OUT_DIR_ABS="$SCRIPT_DIR/$OC_OUT_DIR_ABS"
    OCB_REF_ABS="$OCB_REF"
    [[ "$OCB_REF_ABS" = /* ]] || OCB_REF_ABS="$SCRIPT_DIR/$OCB_REF_ABS"
    [[ -f "$OCB_REF_ABS" ]] || { echo "[overlap_consensus_blast] ERROR: reference not found: $OCB_REF_ABS" >&2; exit 1; }

    OCB_OUT_DIR="$OC_OUT_DIR_ABS/blast"
    mkdir -p "$OCB_OUT_DIR/db"
    if [[ ! -f "$OCB_OUT_DIR/db/ref.nhr" ]]; then
      makeblastdb -in "$OCB_REF_ABS" -dbtype nucl -out "$OCB_OUT_DIR/db/ref" >/dev/null
    fi

    _blast_one_fa() {
      local fa="$1" tag="$2"
      [[ -f "$fa" ]] || { echo "[overlap_consensus_blast] SKIPPED $tag (no FASTA: $fa)"; return 0; }
      local outfmt7="$OCB_OUT_DIR/${OC_SEQ_ID}_${tag}.blast.tsv"
      local outpair="$OCB_OUT_DIR/${OC_SEQ_ID}_${tag}.pairwise.txt"
      blastn -query "$fa" -db "$OCB_OUT_DIR/db/ref" \
        -evalue "$OCB_EVALUE" -word_size "$OCB_WORD_SIZE" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
        -max_target_seqs 5 -num_threads 1 > "$outfmt7"
      blastn -query "$fa" -db "$OCB_OUT_DIR/db/ref" \
        -evalue "$OCB_EVALUE" -word_size "$OCB_WORD_SIZE" \
        -max_target_seqs 1 -num_threads 1 > "$outpair"
      local top
      top="$(head -n1 "$outfmt7" || true)"
      if [[ -n "$top" ]]; then
        echo "[overlap_consensus_blast] $tag top hit: $top"
      else
        echo "[overlap_consensus_blast] $tag NO_HIT against $OCB_REF_ABS"
      fi
    }

    case "$OCB_WHICH" in
      reference_length) _blast_one_fa "$OC_OUT_DIR_ABS/${OC_SEQ_ID}_reference_length.fa" "reference_length" ;;
      covered_span)     _blast_one_fa "$OC_OUT_DIR_ABS/${OC_SEQ_ID}_covered_span.fa"     "covered_span" ;;
      both)             _blast_one_fa "$OC_OUT_DIR_ABS/${OC_SEQ_ID}_reference_length.fa" "reference_length"
                        _blast_one_fa "$OC_OUT_DIR_ABS/${OC_SEQ_ID}_covered_span.fa"     "covered_span" ;;
      *) echo "[overlap_consensus_blast] ERROR: unknown overlap_consensus_blast.which='$OCB_WHICH' (use reference_length|covered_span|both)" >&2; exit 1 ;;
    esac
    echo "[overlap_consensus_blast] -> $OCB_OUT_DIR"
  fi
fi

# -----------------------------------------------------------------------------
# overlap_consensus_report: build a Markdown summary of the stitching results,
# coverage, and BLAST paralog confirmation. Runs whenever overlap_consensus ran
# (output_dir is set), regardless of whether overlap_consensus_blast also ran.
# -----------------------------------------------------------------------------
if [[ -n "${OC_OUT_DIR:-}" ]]; then
  OC_OUT_DIR_ABS="${OC_OUT_DIR_ABS:-$OC_OUT_DIR}"
  [[ "$OC_OUT_DIR_ABS" = /* ]] || OC_OUT_DIR_ABS="$SCRIPT_DIR/$OC_OUT_DIR_ABS"
  "$PYTHON_BIN" "$SCRIPT_DIR/scripts/overlap_consensus_report.py" \
    --output-dir "$OC_OUT_DIR_ABS" \
    --workspace  "$SCRIPT_DIR"
else
  # OC_OUT_DIR not set means overlap_consensus was skipped; re-read from config.
  _OC_DIR="$("$PYTHON_BIN" -c "
import sys, tomllib
with open(sys.argv[1], 'rb') as fh:
    cfg = tomllib.load(fh)
print((cfg.get('overlap_consensus') or {}).get('output_dir', ''))
" "$CONFIG")"
  if [[ -n "$_OC_DIR" ]]; then
    [[ "$_OC_DIR" = /* ]] || _OC_DIR="$SCRIPT_DIR/$_OC_DIR"
    [[ -d "$_OC_DIR" ]] && "$PYTHON_BIN" "$SCRIPT_DIR/scripts/overlap_consensus_report.py" \
      --output-dir "$_OC_DIR" \
      --workspace  "$SCRIPT_DIR"
  fi
fi

printf '[%s] [f_inside_primer_pipeline.sh] DONE. log: %s\n' "$(date +%H:%M:%S)" "$LOG_FILE"
