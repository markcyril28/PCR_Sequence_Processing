#!/usr/bin/env python3
"""Build a tri-consensus from two Sanger F+R consensuses and the reference.

For each position, take a majority vote across:
  - Sanger batch A (HC01463503)  -- unambiguous A/C/G/T contributes 1 vote
  - Sanger batch B (HC01468262)  -- unambiguous A/C/G/T contributes 1 vote
  - Reference (always contributes 1 vote)

Rules
-----
1. N in a Sanger track casts no vote.
2. IUPAC ambiguity (W/S/M/K/R/Y/B/D/H/V) in a Sanger track casts no vote
   (we keep only confidently elucidated Sanger calls).
3. If both Sanger tracks contribute a vote and they agree with each other,
   the result is that base regardless of the reference (potential variant).
4. Otherwise the reference base wins ties.
5. If both Sanger tracks are N/ambig, the reference base is used.

The output FASTA carries the same length as the inputs (no realignment).
A per-position TSV records the three inputs, vote tallies, the chosen base,
and a class label: AGREE_ALL, SANGER_VARIANT, SANGER_ASSIST_REF, REF_FALLBACK.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alignment_viz import format_alignment  # noqa: E402

UNAMBIG = set("ACGT")


def read_single_fasta(path: Path) -> tuple[str, str]:
    header = ""
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                break
            header = line[1:].strip()
            continue
        parts.append(line.strip())
    return header, "".join(parts).upper()


def classify(a: str, b: str, r: str, chosen: str) -> str:
    a_vote = a in UNAMBIG
    b_vote = b in UNAMBIG
    if a_vote and b_vote and a == b == r:
        return "AGREE_ALL"
    if a_vote and b_vote and a == b and a != r:
        return "SANGER_VARIANT"
    if (a_vote ^ b_vote) and chosen == (a if a_vote else b) == r:
        return "SANGER_ASSIST_REF"
    if not a_vote and not b_vote:
        return "REF_FALLBACK"
    if chosen == r:
        return "REF_TIEBREAK"
    return "OTHER"


def tri_consensus(a: str, b: str, r: str) -> tuple[str, list[tuple[int, str, str, str, str, str]]]:
    assert len(a) == len(b) == len(r), "all three sequences must share length"
    out: list[str] = []
    rows: list[tuple[int, str, str, str, str, str]] = []
    for i, (ai, bi, ri) in enumerate(zip(a, b, r), start=1):
        votes: dict[str, int] = {}
        if ai in UNAMBIG:
            votes[ai] = votes.get(ai, 0) + 1
        if bi in UNAMBIG:
            votes[bi] = votes.get(bi, 0) + 1
        votes[ri] = votes.get(ri, 0) + 1
        # Sanger-supported variant overrides ref
        if ai in UNAMBIG and ai == bi and ai != ri:
            chosen = ai
        else:
            top = max(votes.values())
            top_bases = [base for base, n in votes.items() if n == top]
            chosen = ri if ri in top_bases else sorted(top_bases)[0]
        out.append(chosen)
        rows.append((i, ai, bi, ri, chosen, classify(ai, bi, ri, chosen)))
    return "".join(out), rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sanger-a", required=True, type=Path)
    ap.add_argument("--sanger-b", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-report", required=True, type=Path)
    ap.add_argument(
        "--out-align",
        type=Path,
        help="Optional path for a Tracy-style consensus-vs-reference block alignment "
             "(matches the '*.ref_align.txt' files emitted by 03a_TraceReview).",
    )
    ap.add_argument("--name", default="tri_consensus")
    args = ap.parse_args()

    ha, sa = read_single_fasta(args.sanger_a)
    hb, sb = read_single_fasta(args.sanger_b)
    hr, sr = read_single_fasta(args.reference)
    if not (len(sa) == len(sb) == len(sr)):
        raise SystemExit(
            f"length mismatch: sanger_a={len(sa)} sanger_b={len(sb)} ref={len(sr)}"
        )

    seq, rows = tri_consensus(sa, sb, sr)
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    summary: dict[str, int] = {}
    for _, _, _, _, _, cls in rows:
        summary[cls] = summary.get(cls, 0) + 1

    header = (
        f">{args.name} ref={hr} sanger_a={ha} sanger_b={hb} "
        f"len={len(seq)} agree_all={summary.get('AGREE_ALL', 0)} "
        f"sanger_variant={summary.get('SANGER_VARIANT', 0)} "
        f"sanger_assist={summary.get('SANGER_ASSIST_REF', 0)} "
        f"ref_fallback={summary.get('REF_FALLBACK', 0)} "
        f"ref_tiebreak={summary.get('REF_TIEBREAK', 0)}"
    )
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    args.out_fasta.write_text(f"{header}\n{body}\n")

    with args.out_report.open("w") as fh:
        fh.write("pos\tsanger_a\tsanger_b\treference\tconsensus\tclass\n")
        for pos, ai, bi, ri, ci, cls in rows:
            fh.write(f"{pos}\t{ai}\t{bi}\t{ri}\t{ci}\t{cls}\n")

    if args.out_align is not None:
        args.out_align.parent.mkdir(parents=True, exist_ok=True)
        agree = summary.get("AGREE_ALL", 0) + summary.get("SANGER_ASSIST_REF", 0)
        header = (
            f"# tri_consensus vs reference  ({args.name})\n"
            f"# sanger_a={ha}\n"
            f"# sanger_b={hb}\n"
            f"# reference={hr}\n"
            f"# length={len(seq)}  agree_with_ref={agree}  "
            f"sanger_variant={summary.get('SANGER_VARIANT', 0)}  "
            f"ref_fallback={summary.get('REF_FALLBACK', 0)}  "
            f"ref_tiebreak={summary.get('REF_TIEBREAK', 0)}\n"
        )
        block = format_alignment("PCR_Seq", seq, "Reference", sr)
        args.out_align.write_text(header + block)

    print(f"wrote {args.out_fasta}")
    print(f"wrote {args.out_report}")
    if args.out_align is not None:
        print(f"wrote {args.out_align}")
    for k, v in sorted(summary.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
