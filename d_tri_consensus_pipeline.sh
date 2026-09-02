#!/usr/bin/env bash
# =============================================================================
# d_tri_consensus_pipeline.sh
#
# Builds a majority-rule tri-consensus for every amplicon sequenced in BOTH
# Sanger batches. For each shared amplicon (matched by the `ref=` accession in
# each consensus FASTA header), the orchestrator calls scripts/tri_consensus.py
# with the two batch consensuses and the shared reference, and writes the
# result to a single output folder under III_PipelineOutput/.
#
# Inputs (per batch, configured in d_tri_consensus_pipeline_config.toml):
#   {output_root}/{batch}/low_quality_saver_pipeline/03a_TraceReview/<sample>/<glob>
#     -- default glob *_forward.reviewed.fa restricts input to the forward
#        trace-reviewed read per amplicon (reverse direction excluded).
#   {output_root}/{batch}/low_quality_saver_pipeline/03a_TraceReview/_ref/<ref>.fa
#
# Output:
#   {output_root}/{io.output_dirname}/<ref>_tri_consensus.fa
#   {output_root}/{io.output_dirname}/<ref>_tri_consensus.report.tsv
#   {output_root}/{io.output_dirname}/tri_consensus_summary.tsv
#
# Usage:
#   ./d_tri_consensus_pipeline.sh [d_tri_consensus_pipeline_config.toml]
# =============================================================================
set -Eeuo pipefail

# ---------------------------------------------------------------------------
# Resource configuration: auto-detect cores; override from TOML after load.
# ---------------------------------------------------------------------------
CPU=$(nproc 2>/dev/null || echo 4)
MAX_PARALLEL="$CPU"

# wait_for_slot: block until fewer than $1 background jobs are running.
wait_for_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

CONFIG="${1:-d_tri_consensus_pipeline_config.toml}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$ROOT/scripts"

_pick_python() {
  local cand
  for cand in "${PYTHON:-}" \
              "$HOME/miniconda3/bin/python3" "$HOME/anaconda3/bin/python3" \
              python3.13 python3.12 python3.11 python3 /usr/bin/python3; do
    [[ -z "$cand" ]] && continue
    command -v "$cand" >/dev/null 2>&1 || continue
    "$cand" -c 'import sys, tomllib' >/dev/null 2>&1 || continue
    "$cand" -c 'import sys; sys.exit(0 if sys.version_info>=(3,11) else 1)' 2>/dev/null \
      && { echo "$cand"; return 0; }
  done
  return 1
}
PYTHON="$(_pick_python)" || { echo "ERROR: need Python >= 3.11." >&2; exit 2; }

[[ "$CONFIG" = /* ]] || CONFIG="$ROOT/$CONFIG"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }
[[ -f "$SCRIPTS/read_toml.py" ]]    || { echo "ERROR: missing scripts/read_toml.py" >&2; exit 2; }
[[ -f "$SCRIPTS/tri_consensus.py" ]] || { echo "ERROR: missing scripts/tri_consensus.py" >&2; exit 2; }

LOG_DIR="${LOG_DIR:-$ROOT/logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/tri_consensus_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[log] writing to $LOG_FILE"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# Load TOML into shell vars: KEY=value (uppercased, dots -> _).
eval "$("$PYTHON" "$SCRIPTS/read_toml.py" "$CONFIG" --export)"

# Use TOML blast.threads if set, otherwise fall back to nproc.
BLAST_THREADS="${BLAST_THREADS:-$CPU}"

# Resolve roots and per-batch directories.
OUTPUT_ROOT="$ROOT/$IO_OUTPUT_ROOT"
OUTPUT_DIR="$OUTPUT_ROOT/$IO_OUTPUT_DIRNAME"
INPUTS_DIR="$OUTPUT_DIR/${INPUTS_STASH_DIRNAME:-_inputs}"
[[ -d "$OUTPUT_ROOT" ]] || die "output_root not found: $OUTPUT_ROOT"
mkdir -p "$OUTPUT_DIR" "$INPUTS_DIR"

A_ID="$BATCH_A_ID"
B_ID="$BATCH_B_ID"
A_CONS_DIR="$OUTPUT_ROOT/$BATCH_A_CONSENSUS_DIR"
B_CONS_DIR="$OUTPUT_ROOT/$BATCH_B_CONSENSUS_DIR"
A_REF_DIR="$OUTPUT_ROOT/$BATCH_A_REFERENCE_DIR"
B_REF_DIR="$OUTPUT_ROOT/$BATCH_B_REFERENCE_DIR"
[[ -d "$A_CONS_DIR" ]] || die "batch_a consensus_dir not found: $A_CONS_DIR"
[[ -d "$B_CONS_DIR" ]] || die "batch_b consensus_dir not found: $B_CONS_DIR"
[[ -d "$A_REF_DIR"  ]] || die "batch_a reference_dir not found: $A_REF_DIR"

A_FASTA_GLOB="${BATCH_A_FASTA_GLOB:-${CONSENSUS_FASTA_GLOB:-*.fa}}"
B_FASTA_GLOB="${BATCH_B_FASTA_GLOB:-${CONSENSUS_FASTA_GLOB:-*.fa}}"
FASTA_SUFFIX="${CONSENSUS_FASTA_SUFFIX:-_tri_consensus.fa}"
REPORT_SUFFIX="${CONSENSUS_REPORT_SUFFIX:-_tri_consensus.report.tsv}"
ALIGN_SUFFIX="${CONSENSUS_ALIGN_SUFFIX:-_tri_consensus.ref_align.txt}"
SUMMARY="$OUTPUT_DIR/${CONSENSUS_SUMMARY_FILENAME:-tri_consensus_summary.tsv}"

log "project      : $PROJECT_NAME"
log "output_dir   : $OUTPUT_DIR"
log "batch A      : $A_ID  ($A_CONS_DIR)"
log "batch B      : $B_ID  ($B_CONS_DIR)"
log "reference    : $A_REF_DIR"

# ---------------------------------------------------------------------------
# Index each batch by reference accession parsed from the consensus header.
# A consensus FASTA header looks like:
#   >I_HC01463503|<amplicon> ref=<ref_accession> len=... ...
# We extract <ref_accession> via the read_toml.py-companion helper below.
# ---------------------------------------------------------------------------
extract_ref() {
  # Print the value of the `ref=` token in the FIRST header line, or "" if absent.
  awk '
    /^>/ {
      for (i=1; i<=NF; i++) if ($i ~ /^ref=/) { sub(/^ref=/,"",$i); print $i; exit }
      exit
    }
  ' "$1"
}

declare -A A_BY_REF=()
declare -A B_BY_REF=()

while IFS= read -r -d '' fa; do
  ref="$(extract_ref "$fa")"
  [[ -z "$ref" ]] && { log "  [skip A] no ref= in $(basename "$fa")"; continue; }
  A_BY_REF["$ref"]="$fa"
done < <(find "$A_CONS_DIR" -maxdepth 1 -type f -name "$A_FASTA_GLOB" -print0)

while IFS= read -r -d '' fa; do
  ref="$(extract_ref "$fa")"
  [[ -z "$ref" ]] && { log "  [skip B] no ref= in $(basename "$fa")"; continue; }
  B_BY_REF["$ref"]="$fa"
done < <(find "$B_CONS_DIR" -maxdepth 1 -type f -name "$B_FASTA_GLOB" -print0)

log "amplicons    : A=${#A_BY_REF[@]}  B=${#B_BY_REF[@]}"

printf 'ref_accession\tsanger_a\tsanger_b\treference\tconsensus_fa\tagree_all\tsanger_variant\tsanger_assist\tref_fallback\tref_tiebreak\tstatus\n' > "$SUMMARY"

# Manifest of copied input FASTAs (records the original filepath of every
# sanger_a / sanger_b / reference consumed by tri_consensus.py).
INPUTS_MANIFEST="$INPUTS_DIR/${INPUTS_MANIFEST_FILENAME:-inputs_manifest.tsv}"
printf 'ref_accession\trole\tcopied_to\toriginal_path\n' > "$INPUTS_MANIFEST"

# stash_input <ref> <role> <src_fa>
#   Copies src_fa into _inputs/ as <ref>.<role>.<basename> and appends a row to
#   the manifest. Role is one of: sanger_a, sanger_b, reference.
stash_input() {
  local ref="$1" role="$2" src="$3"
  [[ -z "$src" || ! -f "$src" ]] && return 0
  local dst="$INPUTS_DIR/${ref}.${role}.$(basename "$src")"
  cp -f "$src" "$dst"
  printf '%s\t%s\t%s\t%s\n' "$ref" "$role" "${dst#$ROOT/}" "${src#$ROOT/}" >> "$INPUTS_MANIFEST"
}

shared=0
skipped=0

# Temp directory for per-amplicon summary rows (parallel writes).
TRI_TMP="$OUTPUT_DIR/.tri_rows"
mkdir -p "$TRI_TMP"

for ref in "${!A_BY_REF[@]}"; do
  a_fa="${A_BY_REF[$ref]}"
  b_fa="${B_BY_REF[$ref]:-}"
  if [[ -z "$b_fa" ]]; then
    log "  [skip] $ref present in A only"
    printf '%s\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t%s\n' \
      "$ref" "$(basename "$a_fa")" "a_only" >> "$SUMMARY"
    stash_input "$ref" "sanger_a" "$a_fa"
    skipped=$((skipped+1))
    continue
  fi
  ref_fa="$A_REF_DIR/${ref}.fa"
  if [[ ! -f "$ref_fa" && -d "$B_REF_DIR" ]]; then
    ref_fa="$B_REF_DIR/${ref}.fa"
  fi
  if [[ ! -f "$ref_fa" ]]; then
    log "  [skip] $ref reference FASTA not found in $A_REF_DIR or $B_REF_DIR"
    printf '%s\t%s\t%s\t-\t-\t-\t-\t-\t-\t-\t%s\n' \
      "$ref" "$(basename "$a_fa")" "$(basename "$b_fa")" "no_reference" >> "$SUMMARY"
    stash_input "$ref" "sanger_a" "$a_fa"
    stash_input "$ref" "sanger_b" "$b_fa"
    skipped=$((skipped+1))
    continue
  fi

  out_fa="$OUTPUT_DIR/${ref}${FASTA_SUFFIX}"
  out_rep="$OUTPUT_DIR/${ref}${REPORT_SUFFIX}"
  out_aln="$OUTPUT_DIR/${ref}${ALIGN_SUFFIX}"
  row_file="$TRI_TMP/${ref}.row"

  log "  [build] $ref"
  stash_input "$ref" "sanger_a"  "$a_fa"
  stash_input "$ref" "sanger_b"  "$b_fa"
  stash_input "$ref" "reference" "$ref_fa"

  # Capture variables used inside the subshell before dispatching.
  _ref="$ref"; _a_fa="$a_fa"; _b_fa="$b_fa"; _ref_fa="$ref_fa"
  _out_fa="$out_fa"; _out_rep="$out_rep"; _out_aln="$out_aln"
  _A_ID="$A_ID"; _B_ID="$B_ID"; _row_file="$row_file"

  wait_for_slot "$MAX_PARALLEL"
  (
    stats="$("$PYTHON" "$SCRIPTS/tri_consensus.py" \
      --sanger-a   "$_a_fa" \
      --sanger-b   "$_b_fa" \
      --reference  "$_ref_fa" \
      --out-fasta  "$_out_fa" \
      --out-report "$_out_rep" \
      --out-align  "$_out_aln" \
      --name       "${_ref}|${_A_ID}+${_B_ID}_tri_consensus" | tail -n +4)"

    # tri_consensus.py prints "  CLASS: N" lines after the two "wrote" lines.
    agree_all=$(printf '%s\n' "$stats" | awk '$1=="AGREE_ALL:" {print $2}')
    sanger_var=$(printf '%s\n' "$stats" | awk '$1=="SANGER_VARIANT:" {print $2}')
    sanger_ast=$(printf '%s\n' "$stats" | awk '$1=="SANGER_ASSIST_REF:" {print $2}')
    ref_fb=$(printf '%s\n' "$stats" | awk '$1=="REF_FALLBACK:" {print $2}')
    ref_tb=$(printf '%s\n' "$stats" | awk '$1=="REF_TIEBREAK:" {print $2}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$_ref" "$(basename "$_a_fa")" "$(basename "$_b_fa")" "$(basename "$_ref_fa")" \
      "$(basename "$_out_fa")" \
      "${agree_all:-0}" "${sanger_var:-0}" "${sanger_ast:-0}" \
      "${ref_fb:-0}" "${ref_tb:-0}" "merged" > "$_row_file"
  ) &
  shared=$((shared+1))
done

# Wait for all background tri_consensus jobs before writing the summary.
wait

# Merge per-amplicon rows into SUMMARY (sorted for determinism).
while IFS= read -r row_file; do
  cat "$row_file" >> "$SUMMARY"
done < <(find "$TRI_TMP" -maxdepth 1 -name '*.row' | sort)
rm -rf "$TRI_TMP"

# Report references that exist only in batch B (A-only handled above).
for ref in "${!B_BY_REF[@]}"; do
  if [[ -z "${A_BY_REF[$ref]:-}" ]]; then
    log "  [skip] $ref present in B only"
    printf '%s\t-\t%s\t-\t-\t-\t-\t-\t-\t-\t%s\n' \
      "$ref" "$(basename "${B_BY_REF[$ref]}")" "b_only" >> "$SUMMARY"
    stash_input "$ref" "sanger_b" "${B_BY_REF[$ref]}"
    skipped=$((skipped+1))
  fi
done

# ---------------------------------------------------------------------------
# BLAST step -- re-confirm each tri-consensus best-hits its expected ref
# against the multi-gene DMP reference, mirroring the blastn step in
# b_low_quality_saver_pipeline.sh.
# ---------------------------------------------------------------------------
if [[ "${BLAST_ENABLED:-false}" == "true" && "$shared" -gt 0 ]]; then
  BLAST_REF="$ROOT/$BLAST_REFERENCE"
  [[ -f "$BLAST_REF" ]] || die "blast reference not found: $BLAST_REF"
  command -v makeblastdb >/dev/null || die "makeblastdb not found (run setup.sh)"
  command -v blastn      >/dev/null || die "blastn not found (run setup.sh)"

  BLAST_DIR="$OUTPUT_DIR/blast"
  DB_DIR="$BLAST_DIR/db"
  mkdir -p "$DB_DIR"
  if [[ ! -f "$DB_DIR/ref.nhr" ]]; then
    log "[blast] makeblastdb -> $DB_DIR/ref"
    makeblastdb -in "$BLAST_REF" -dbtype nucl -out "$DB_DIR/ref" >/dev/null
  fi

  ASSIGN_TSV="$OUTPUT_DIR/${BLAST_ASSIGNMENTS_FILENAME:-tri_consensus_blast.tsv}"
  printf 'consensus\texpected_ref\tbest_ref\tpident\tqcov\tevalue\tstatus\n' > "$ASSIGN_TSV"

  while IFS= read -r -d '' qfa; do
    base="$(basename "$qfa")"
    expected_ref="${base%${FASTA_SUFFIX}}"
    hit=$(blastn -query "$qfa" -db "$DB_DIR/ref" \
            -evalue "$BLAST_EVALUE" -word_size "$BLAST_WORD_SIZE" \
            -outfmt '6 sseqid pident qcovs evalue' \
            -max_target_seqs 1 -num_threads "${BLAST_THREADS:-$CPU}" | head -n1 || true)
    if [[ -z "$hit" ]]; then
      printf '%s\t%s\tNO_HIT\t-\t-\t-\tno_hit\n' "$base" "$expected_ref" >> "$ASSIGN_TSV"
      log "  [blast] $base : NO_HIT"
      continue
    fi
    sseq=$(awk '{print $1}' <<<"$hit")
    pid=$(awk  '{print $2}' <<<"$hit")
    qcv=$(awk  '{print $3}' <<<"$hit")
    ev=$(awk   '{print $4}' <<<"$hit")
    pass=$(awk -v p="$pid" -v c="$qcv" -v mp="$BLAST_MIN_PIDENT" -v mc="$BLAST_MIN_QCOV" \
           'BEGIN{print (p+0>=mp+0 && c+0>=mc+0)?"yes":"no"}')
    if [[ "$pass" != "yes" ]]; then
      status="fail_threshold"
    elif [[ "$sseq" == "$expected_ref" ]]; then
      status="match"
    else
      status="mismatch_ref"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$base" "$expected_ref" "$sseq" "$pid" "$qcv" "$ev" "$status" >> "$ASSIGN_TSV"
    log "  [blast] $base : $sseq pident=$pid qcov=$qcv ($status)"
  done < <(find "$OUTPUT_DIR" -maxdepth 1 -type f -name "*${FASTA_SUFFIX}" -print0)

  log "[blast] assignments -> $ASSIGN_TSV"
fi

log ""
log "DONE. shared_amplicons=$shared  skipped=$skipped"
log "outputs under: $OUTPUT_DIR"
log "summary      : $SUMMARY"
