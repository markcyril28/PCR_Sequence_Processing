#!/usr/bin/env bash
# =============================================================================
# setup.sh: install the PCR_Sequence pipeline dependencies.
#
# Creates (or reuses) a conda env named `pcr_sanger` with:
#   - python>=3.11, biopython, tomli (fallback for py<3.11)
#   - blast (provides makeblastdb/blastn)
#   - tracy (for Method 3 and the optional het-indel rescue in Method 1)
#
# Usage:
#   bash setup.sh                  # create env + check tools
#   bash setup.sh --env myenv      # custom env name
#   bash setup.sh --no-conda       # pip-only install into the active interpreter
# =============================================================================
set -Eeuo pipefail

ENV_NAME="pcr_sanger"
USE_CONDA=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env)      ENV_NAME="$2"; shift 2 ;;
    --no-conda) USE_CONDA=0; shift ;;
    -h|--help)  sed -n '2,15p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

have() { command -v "$1" >/dev/null 2>&1; }

if [[ "$USE_CONDA" -eq 1 ]] && (have mamba || have conda); then
  PM="$(have mamba && echo mamba || echo conda)"
  echo "[setup] Using $PM env: $ENV_NAME"

  # Locate the env prefix (if it exists) so we can probe for required tools
  # without invoking `mamba install` (which prints noisy 'unknown command
  # "inspect"' errors on some mamba versions even though the transaction
  # actually succeeds).
  ENV_PREFIX="$(conda env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n {print $NF}' | head -n1)"

  env_has_all_tools() {
    [[ -n "$ENV_PREFIX" && -d "$ENV_PREFIX/bin" ]] || return 1
    for t in python blastn makeblastdb tracy; do
      [[ -x "$ENV_PREFIX/bin/$t" ]] || return 1
    done
    "$ENV_PREFIX/bin/python" -c "import Bio" >/dev/null 2>&1 || return 1
    return 0
  }

  if [[ -z "$ENV_PREFIX" ]]; then
    echo "[setup] Creating env $ENV_NAME"
    "$PM" create -y -n "$ENV_NAME" -c bioconda -c conda-forge \
        python=3.11 biopython blast tracy
  elif env_has_all_tools; then
    echo "[setup] Env $ENV_NAME already has all required tools; skipping install"
  else
    echo "[setup] Env $ENV_NAME exists but is missing tools; installing"
    "$PM" install -y -n "$ENV_NAME" -c bioconda -c conda-forge \
        biopython blast tracy
  fi
  echo
  echo "[setup] Activate with:  conda activate $ENV_NAME"
else
  echo "[setup] No conda detected (or --no-conda): falling back to pip"
  python -m pip install --upgrade pip
  python -m pip install biopython tomli
  echo
  echo "[setup] WARNING: BLAST+ and Tracy are NOT pip-installable."
  echo "        Install them manually:"
  echo "          - BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
  echo "          - Tracy : https://github.com/gear-genomics/tracy/releases"
fi

echo
echo "[setup] Tool check:"
for t in python blastn makeblastdb tracy; do
  if have "$t"; then
    printf "  OK   %-12s -> %s\n" "$t" "$(command -v "$t")"
  else
    printf "  MISS %-12s (needed for some pipeline steps)\n" "$t"
  fi
done

# Make every script executable (no-op on Windows but harmless).
chmod +x "$(dirname "$0")"/*.sh 2>/dev/null || true
find "$(dirname "$0")/scripts" -type f \( -name '*.py' -o -name '*.sh' \) \
  -exec chmod +x {} + 2>/dev/null || true

echo
echo "[setup] Done. Next:"
echo "  1. Edit *_config.toml to choose options"
echo "  2. Run a pipeline, e.g.:"
echo "       bash b_low_quality_saver_pipeline.sh b_low_quality_saver_pipeline_config.toml"
echo "       bash a_merge_consensus_pipeline.sh   a_merge_consensus_config.toml"
