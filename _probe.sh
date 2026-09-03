#!/usr/bin/env bash
echo "=== candidate pythons (need >=3.11 + tomllib + Bio) ==="
for p in \
  "$HOME/miniconda3/bin/python3" \
  "$HOME/miniconda3/envs/pcr_sanger/bin/python3" \
  "$HOME/miniconda3/envs/egg/bin/python3" \
  "$HOME/miniconda3/envs/nat/bin/python3" \
  "$HOME/miniconda3/envs/gea/bin/python3" \
  /usr/bin/python3; do
  [ -x "$p" ] || continue
  v=$("$p" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null)
  bio=$("$p" -c 'import tomllib,Bio;print("OK")' 2>/dev/null)
  echo "$p -> py=${v:-?} bio=${bio:-NO}"
done
echo "=== pcr_sanger tools ==="
for t in "$HOME/miniconda3/envs/pcr_sanger/bin"/{blastn,makeblastdb,tracy}; do
  [ -x "$t" ] && echo "OK $t" || echo "MISSING $t"
done
