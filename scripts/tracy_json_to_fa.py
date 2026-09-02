#!/usr/bin/env python3
"""Extract per-base FASTA + quality scores from a Tracy JSON output.

`tracy basecall -f json` (or the JSON emitted alongside `tracy align`) contains
the trimmed primary basecall and per-base Phred quality. This helper writes:
    <stem>.fa   single-record FASTA of the basecalled sequence
    <stem>.qual whitespace-separated Phred integers (one per base)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, type=Path)
    ap.add_argument("--out-fa", required=True, type=Path)
    ap.add_argument("--out-qual", required=True, type=Path)
    ap.add_argument("--seq-id", required=True)
    args = ap.parse_args()

    data = json.loads(args.json.read_text())

    # Tracy JSON layout has evolved; try the common keys.
    seq = (
        data.get("primarySeq")
        or data.get("basecalls")
        or data.get("sequence")
    )
    quals = (
        data.get("qual")
        or data.get("quality")
        or data.get("phred")
    )

    if not seq or not quals:
        # Fall back to nested location used by recent Tracy releases.
        bc = data.get("basecall", {}) if isinstance(data.get("basecall"), dict) else {}
        seq = seq or bc.get("primarySeq") or bc.get("seq")
        quals = quals or bc.get("qual") or bc.get("quality")

    if not seq or not quals:
        print(f"ERROR: could not find primary sequence/quality in {args.json}",
              file=sys.stderr)
        return 1

    if isinstance(quals, str):
        quals = quals.split()
    quals = [int(q) for q in quals]

    args.out_fa.parent.mkdir(parents=True, exist_ok=True)
    with args.out_fa.open("w") as fh:
        fh.write(f">{args.seq_id}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    args.out_qual.write_text(" ".join(str(q) for q in quals) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
