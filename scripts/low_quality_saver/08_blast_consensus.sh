#!/usr/bin/env bash
# Runs blastn on every *_consensus.fa under <in_dir> against the multi-gene
# DMP reference and emits per-query pairwise alignment + tabular reports,
# plus a single aggregated assignments TSV.
#
# Usage:
#   scripts/low_quality_saver/08_blast_consensus.sh <in_dir> <out_dir> [reference.fa]
#
# Designed to be invoked from inside the pcr_sanger conda env (which provides
# makeblastdb / blastn). Safe to re-run; reuses any DB already built under
# <out_dir>/db.
set -Eeuo pipefail

IN_DIR="${1:?usage: 08_blast_consensus.sh <in_dir> <out_dir> [reference.fa]}"
OUT_DIR="${2:?usage: 08_blast_consensus.sh <in_dir> <out_dir> [reference.fa]}"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
REF="${3:-$ROOT/II_QualityClassified/PF05078_DMP_transcripts_cdhit.fa}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
EVALUE="${EVALUE:-1e-5}"
WORD_SIZE="${WORD_SIZE:-9}"

[[ -d "$IN_DIR" ]] || { echo "ERROR: in_dir not found: $IN_DIR" >&2; exit 1; }
[[ -f "$REF"    ]] || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
command -v blastn      >/dev/null || { echo "ERROR: blastn not on PATH" >&2; exit 2; }
command -v makeblastdb >/dev/null || { echo "ERROR: makeblastdb not on PATH" >&2; exit 2; }

DB_DIR="$OUT_DIR/db"
mkdir -p "$DB_DIR"
if [[ ! -f "$DB_DIR/ref.nhr" ]]; then
  echo "[blast] makeblastdb -> $DB_DIR/ref"
  makeblastdb -in "$REF" -dbtype nucl -out "$DB_DIR/ref" >/dev/null
fi

ASSIGN="$OUT_DIR/assignments.tsv"
printf 'sample\tconsensus\tbest_ref\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcov\n' > "$ASSIGN"

# Discover consensus FASTAs at two layouts:
#   1) flat:           <IN_DIR>/*_consensus.fa
#   2) per-sample:     <IN_DIR>/<sample>/*_consensus.fa   (LQS pipeline layout)
shopt -s nullglob globstar
fastas=("$IN_DIR"/**/*_consensus.fa)
shopt -u nullglob globstar
[[ ${#fastas[@]} -gt 0 ]] || { echo "WARN: no *_consensus.fa under $IN_DIR" >&2; exit 0; }

for fa in "${fastas[@]}"; do
  stem="$(basename "$fa" .fa)"
  # Mirror sample subdirectory into OUT_DIR. If consensus is directly under
  # IN_DIR (flat layout), sample is "" and outputs land directly in OUT_DIR.
  rel_dir="$(dirname "${fa#$IN_DIR/}")"
  if [[ "$rel_dir" == "." ]]; then
    sample=""
    out_sub="$OUT_DIR"
  else
    sample="$rel_dir"
    out_sub="$OUT_DIR/$rel_dir"
    mkdir -p "$out_sub"
  fi
  echo "[blast] ${sample:+$sample/}$stem"
  tab="$out_sub/${stem}.blast.tsv"
  aln="$out_sub/${stem}.blast_align.txt"
  blastn -query "$fa" -db "$DB_DIR/ref" \
    -evalue "$EVALUE" -word_size "$WORD_SIZE" \
    -outfmt '6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
    -max_target_seqs 5 -num_threads "$THREADS" > "$tab"
  blastn -query "$fa" -db "$DB_DIR/ref" \
    -evalue "$EVALUE" -word_size "$WORD_SIZE" \
    -outfmt 0 -max_target_seqs 3 -num_threads "$THREADS" > "$aln"
  top="$(head -n1 "$tab" || true)"
  if [[ -n "$top" ]]; then
    printf '%s\t%s\t%s\n' "$sample" "$stem" "$top" >> "$ASSIGN"
  else
    printf '%s\t%s\tNO_HIT\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n' "$sample" "$stem" >> "$ASSIGN"
  fi
done

echo "[blast] done. outputs in $OUT_DIR"
echo "[blast] assignments: $ASSIGN"
