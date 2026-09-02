#!/usr/bin/env python3
"""Tiny TOML reader for bash pipelines.

Usage:
    python read_toml.py <config.toml> <dotted.key>
    python read_toml.py <config.toml> --export        # KEY=value lines (uppercased, dots->_)
    python read_toml.py <config.toml> --batches       # TSV of [[batches]] table with defaults

The flat exporter skips any array-of-tables sections (e.g. `[[batches]]`) so
they do not produce malformed KEY='...' lines. Use --batches to consume them.

Requires Python 3.11+ (uses stdlib `tomllib`). Falls back to `tomli` on older Pythons.
"""
from __future__ import annotations

import sys
from pathlib import Path

try:
    import tomllib  # py311+
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


def _flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten(v, key))
        else:
            # Skip array-of-tables (list of dicts); they are consumed by --batches
            # or other table-aware modes, not by the flat exporter.
            if isinstance(v, list) and v and all(isinstance(item, dict) for item in v):
                continue
            out[key] = v
    return out


# Columns emitted by --batches, in order. Anything not provided in a batch
# entry falls back to the corresponding default from the listed sections.
_BATCH_COLUMNS = [
    ("id",             []),                                   # required, no default
    ("sample_dir",     []),                                   # required, no default
    ("forward_glob",   ["reads.forward_glob"]),
    ("reverse_glob",   ["reads.reverse_glob"]),
    ("forward_suffix", ["reads.forward_suffix"]),
    ("reverse_suffix", ["reads.reverse_suffix"]),
    ("forward_token",  ["reads.forward_token"]),
    ("reverse_token",  ["reads.reverse_token"]),
]


def _lookup(flat: dict, dotted: str) -> str:
    val = flat.get(dotted)
    if val is None:
        return ""
    if isinstance(val, bool):
        return "true" if val else "false"
    return str(val)


def _emit_batches(data: dict, flat: dict) -> int:
    batches = data.get("batches")
    if not batches:
        print("ERROR: no [[batches]] table in config", file=sys.stderr)
        return 1
    if not isinstance(batches, list) or not all(isinstance(b, dict) for b in batches):
        print("ERROR: [[batches]] must be an array of tables", file=sys.stderr)
        return 1

    # active_batches may live under [project] (LQS pipeline) or [io] (merge
    # pipeline); fall back to a top-level key for ad-hoc configs. An empty
    # list (or missing key) means "process every defined batch".
    active = (
        (data.get("project") or {}).get("active_batches")
        or (data.get("io") or {}).get("active_batches")
        or data.get("active_batches")
        or []
    )
    if active and not isinstance(active, list):
        print("ERROR: active_batches must be a TOML array", file=sys.stderr)
        return 1
    active_set = {str(x) for x in active} if active else None

    header = [c for c, _ in _BATCH_COLUMNS]
    print("\t".join(header))
    for b in batches:
        if active_set is not None and str(b.get("id", "")) not in active_set:
            continue
        row: list[str] = []
        for col, defaults in _BATCH_COLUMNS:
            if col in b and b[col] is not None:
                val = b[col]
                val = "true" if val is True else "false" if val is False else str(val)
            else:
                val = ""
                for dotted in defaults:
                    val = _lookup(flat, dotted)
                    if val:
                        break
            row.append(val)
        # Required fields must be non-empty.
        if not row[0] or not row[1]:
            print(f"ERROR: batch missing required 'id' or 'sample_dir': {b}", file=sys.stderr)
            return 1
        print("\t".join(row))
    return 0


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2

    cfg_path = Path(sys.argv[1])
    arg = sys.argv[2]

    with cfg_path.open("rb") as fh:
        data = tomllib.load(fh)

    flat = _flatten(data)

    if arg == "--export":
        for k, v in flat.items():
            shell_key = k.upper().replace(".", "_")
            if isinstance(v, bool):
                shell_val = "true" if v else "false"
            elif isinstance(v, list):
                shell_val = " ".join(str(item) for item in v)
            else:
                shell_val = str(v)
            # Quote with single quotes; escape any single quotes in the value.
            shell_val = shell_val.replace("'", "'\\''")
            print(f"{shell_key}='{shell_val}'")
        return 0

    if arg == "--batches":
        return _emit_batches(data, flat)

    if arg not in flat:
        print(f"Key not found: {arg}", file=sys.stderr)
        return 1

    val = flat[arg]
    if isinstance(val, bool):
        print("true" if val else "false")
    else:
        print(val)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
