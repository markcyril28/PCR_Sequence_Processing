#!/usr/bin/env python3
"""Shared utility: Tracy-style visual pairwise alignment formatter.

Usage (from another script in the same directory):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from alignment_viz import format_alignment
"""
from __future__ import annotations

# ── configurable defaults ────────────────────────────────────────────────────
BLOCK_WIDTH = 50   # sequence characters per alignment block (matches Tracy)
# ────────────────────────────────────────────────────────────────────────────


def format_alignment(
    name1: str,
    seq1: str,
    name2: str,
    seq2: str,
    score: float | None = None,
    width: int = BLOCK_WIDTH,
) -> str:
    """Return a Tracy-style visual alignment string.

    *seq1* and *seq2* must be the already-gapped (aligned) sequences of equal
    length.  Position counters are 1-based counts over non-gap characters.

    Output format (one block shown):

        PCR_Seq           1 AAATAAGAATCTAGAATGATACT...
                                                      ||||||
        Reference         1 ----------------------ATGGAC...
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            f"aligned sequences must have the same length "
            f"({len(seq1)} vs {len(seq2)})"
        )

    lw = max(len(name1), len(name2), 3)   # label column width
    pw = 10                                # position column width
    prefix = " " * (lw + pw + 1)          # leading blanks for the match line

    lines: list[str] = []
    if score is not None:
        lines.append(f"Alignment score: {score:.0f}")
    lines.append("#" + "-" * 63)
    lines.append("")

    pos1 = pos2 = 1
    for i in range(0, len(seq1), width):
        chunk1 = seq1[i : i + width]
        chunk2 = seq2[i : i + width]

        match = "".join(
            "|" if c1 != "-" and c2 != "-" and c1.upper() == c2.upper() else " "
            for c1, c2 in zip(chunk1, chunk2)
        )

        lines.append(f"{name1:<{lw}}{pos1:>{pw}} {chunk1}")
        lines.append(f"{prefix}{match}")
        lines.append(f"{name2:<{lw}}{pos2:>{pw}} {chunk2}")
        lines.append("")

        pos1 += sum(1 for c in chunk1 if c != "-")
        pos2 += sum(1 for c in chunk2 if c != "-")

    lines.append("#" + "-" * 63)
    return "\n".join(lines) + "\n"
