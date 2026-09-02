#!/usr/bin/env python3
"""Convert an ABI .ab1 trace to a single-record FASTA (primary basecall).

Optionally end-trims the read with Mott's algorithm before writing the FASTA
(--trim-q > 0). Trimming improves BLAST query quality by removing low-Q
flanks that the basecaller did not mark as N.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trim_reads import slice_trim


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ab1",    required=True, type=Path)
    ap.add_argument("--out",    required=True, type=Path)
    ap.add_argument("--seq-id", required=True)
    ap.add_argument("--trim-q",       type=int, default=0,
                    help="Mott trim quality cutoff (0 = disable trimming).")
    ap.add_argument("--trim-min-len", type=int, default=50,
                    help="Skip trimming if it would leave fewer than this many bases.")
    args = ap.parse_args()

    rec = next(SeqIO.parse(str(args.ab1), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in ABI file: {args.ab1}")
    s = str(rec.seq)
    if args.trim_q > 0:
        q = list(rec.letter_annotations.get("phred_quality", [0] * len(s)))
        s, _q, lo, hi = slice_trim(s, q, q_cutoff=args.trim_q,
                                   min_len=args.trim_min_len)
        trim_tag = f" trim=[{lo},{hi})/{len(rec.seq)}"
    else:
        trim_tag = ""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write(f">{args.seq_id}{trim_tag}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
