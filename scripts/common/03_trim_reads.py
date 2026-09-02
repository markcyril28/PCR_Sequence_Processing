#!/usr/bin/env python3
"""End-trim Sanger reads using Mott's modified algorithm.

Sanger basecallers do not always emit `N` for low-quality 5' / 3' tails; they
also call low-confidence A/C/G/T bases that propagate into consensus calls,
BLAST hits, and alignment frames. Internal Ns are harmless (alignment treats
them as wildcards), but noisy END regions can:
    * shift global-alignment frames around junk tails,
    * cause merged consensus to pick a low-Q miscall over the high-Q strand,
    * inflate length / coverage stats with unusable sequence.

This module exposes a single function `mott_trim()` that returns the
[start, end) interval of a read to KEEP. Consumers should then either:
    (a) slice the seq + qual to that interval, OR
    (b) MASK the bases outside the interval (set base to N, Q to 0) so that
        peak-location indices and per-channel trace arrays remain aligned to
        the original basecall positions.

The trim is conservative: if the read is short or the trimmed interval would
fall below `min_len`, the FULL read is returned unchanged.

CLI (diagnostic):
    python scripts/common/03_trim_reads.py read.ab1 --q 20
"""
from __future__ import annotations

# =============================================================================
# Defaults: kept at the top so they are easy to find/edit.
# =============================================================================
DEFAULT_Q_CUTOFF = 20      # Mott's quality cutoff (Phred)
DEFAULT_MIN_LEN  = 50      # never trim below this many bases (read kept as-is)
DEFAULT_MASK_BASE = "N"    # symbol used when masking trimmed flanks

from pathlib import Path
from typing import Sequence


def mott_trim(qualities: Sequence[int],
              q_cutoff: int = DEFAULT_Q_CUTOFF,
              min_len: int = DEFAULT_MIN_LEN) -> tuple[int, int]:
    """Return (start, end): the half-open interval to KEEP.

    Implements Mott's modified trimming algorithm: assign each base a score
    of `q_cutoff - q[i]`, then find the contiguous subsequence whose score
    SUM is most negative (i.e. highest summed quality above the cutoff).
    Equivalently: maximum-sum subarray of `q[i] - q_cutoff`.

    If no positive-score subarray exists or the resulting window is shorter
    than `min_len`, the full read is returned unchanged.
    """
    n = len(qualities)
    if n == 0:
        return (0, 0)
    if n <= min_len:
        return (0, n)

    best_sum = 0.0
    cur_sum = 0.0
    best_lo = best_hi = 0
    cur_lo = 0
    for i, q in enumerate(qualities):
        cur_sum += (q - q_cutoff)
        if cur_sum < 0:
            cur_sum = 0
            cur_lo = i + 1
        elif cur_sum > best_sum:
            best_sum = cur_sum
            best_lo = cur_lo
            best_hi = i + 1

    if best_hi - best_lo < min_len:
        # Trimmed window collapsed; keep the full read so the downstream
        # caller can still try (it may still produce a consensus or be
        # filtered by a different rule).
        return (0, n)
    return (best_lo, best_hi)


def slice_trim(seq: str, quals: Sequence[int],
               q_cutoff: int = DEFAULT_Q_CUTOFF,
               min_len: int = DEFAULT_MIN_LEN
               ) -> tuple[str, list[int], int, int]:
    """Return (trimmed_seq, trimmed_quals, start, end).

    Use for consumers that DO NOT need to preserve original base indices
    (e.g. writing a FASTA for BLAST). For consumers that index into trace
    arrays / peak locations, prefer `mask_trim()` instead.
    """
    lo, hi = mott_trim(quals, q_cutoff=q_cutoff, min_len=min_len)
    return seq[lo:hi], list(quals[lo:hi]), lo, hi


def mask_trim(seq: str, quals: Sequence[int],
              q_cutoff: int = DEFAULT_Q_CUTOFF,
              min_len: int = DEFAULT_MIN_LEN,
              mask_base: str = DEFAULT_MASK_BASE
              ) -> tuple[str, list[int], int, int]:
    """Return (masked_seq, masked_quals, start, end).

    Bases OUTSIDE [start, end) are replaced with `mask_base` and their
    quality is set to 0, so any downstream Q-floor (`min_q`) check will
    treat them as N without disturbing peak-location / trace indices.
    """
    lo, hi = mott_trim(quals, q_cutoff=q_cutoff, min_len=min_len)
    n = len(seq)
    if (lo, hi) == (0, n):
        return seq, list(quals), lo, hi
    out_seq = mask_base * lo + seq[lo:hi] + mask_base * (n - hi)
    out_q   = [0] * lo + list(quals[lo:hi]) + [0] * (n - hi)
    return out_seq, out_q, lo, hi


# ---------------------------------------------------------------------------
# CLI: quick inspection of how a single .ab1 would be trimmed.
# ---------------------------------------------------------------------------
def _cli() -> int:
    import argparse
    from Bio import SeqIO

    p = argparse.ArgumentParser(description="Inspect Mott trim of an .ab1.")
    p.add_argument("ab1", type=Path)
    p.add_argument("--q",       type=int, default=DEFAULT_Q_CUTOFF)
    p.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN)
    args = p.parse_args()

    rec = next(SeqIO.parse(str(args.ab1), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in {args.ab1}")
    seq = str(rec.seq)
    q   = list(rec.letter_annotations.get("phred_quality", [0] * len(seq)))
    lo, hi = mott_trim(q, q_cutoff=args.q, min_len=args.min_len)
    kept = hi - lo
    mean_full = sum(q) / max(1, len(q))
    mean_kept = sum(q[lo:hi]) / max(1, kept)
    print(f"file       : {args.ab1.name}")
    print(f"length     : {len(seq)}")
    print(f"trim window: [{lo}, {hi})  -> kept={kept}  trimmed_5'={lo}  trimmed_3'={len(seq) - hi}")
    print(f"mean Q full / kept : {mean_full:.1f}  /  {mean_kept:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
