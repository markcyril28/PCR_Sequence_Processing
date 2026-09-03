#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/a_merge_consensus_pipeline.sh"
bash "$SCRIPT_DIR/b_low_quality_saver_pipeline.sh"
bash "$SCRIPT_DIR/c_build_report.sh"
bash "$SCRIPT_DIR/d_tri_consensus_pipeline.sh"
bash "$SCRIPT_DIR/e_collate_good_sequences.sh"
