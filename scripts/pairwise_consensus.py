#!/usr/bin/env python3
"""Reference-free F+R Sanger consensus merger.

Reads a forward .ab1 and a reverse .ab1, reverse-complements the reverse read
(and reverses its quality vector), pairwise-aligns the two strands, and emits
a single consensus FASTA. Used by merge_pipeline.sh.

Consensus rules (per aligned column):
    * Both gaps                    -> skip
    * One gap, other base above min_q -> take the base
    * One gap, other base below min_q -> take the base, lowercase flag (kept as N)
    * Both bases agree             -> take that base
    * Both bases disagree          -> resolve by `--disagreement`:
          higher_q -> base from the strand with higher Phred Q
          iupac    -> IUPAC ambiguity code
          n        -> 'N'
    * If both strands' Q < min_q   -> 'N'
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alignment_viz import format_alignment
from trim_reads import mask_trim

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}


def _read_ab1(path: Path) -> tuple[str, list[int]]:
    rec = next(SeqIO.parse(str(path), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in ABI file: {path}")
    seq = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality", [0] * len(seq)))
    if len(quals) != len(seq):
        quals = (quals + [0] * len(seq))[: len(seq)]
    return seq, quals


def _rev_comp(seq: str, quals: list[int]) -> tuple[str, list[int]]:
    return str(Seq(seq).reverse_complement()), list(reversed(quals))


def _align(f_seq: str, r_seq: str,
           match: int, mismatch: int, gap_open: int, gap_extend: int):
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = match
    a.mismatch_score = mismatch
    a.open_gap_score = gap_open
    a.extend_gap_score = gap_extend
    return a.align(f_seq, r_seq)[0]


def _merge(aln,
           f_seq: str, f_q: list[int],
           r_seq: str, r_q: list[int],
           min_q: int, policy: str) -> tuple[str, dict]:
    a_f, a_r = str(aln[0]), str(aln[1])
    fi = ri = 0
    cons: list[str] = []
    stats = {"len": 0, "agree": 0, "disagree": 0, "f_only": 0, "r_only": 0,
             "n_calls": 0, "overlap": 0}

    for cf, cr in zip(a_f, a_r):
        fb = cf if cf != "-" else None
        rb = cr if cr != "-" else None
        fq = f_q[fi] if fb is not None and fi < len(f_q) else 0
        rq = r_q[ri] if rb is not None and ri < len(r_q) else 0
        if fb is not None:
            fi += 1
        if rb is not None:
            ri += 1

        if fb is None and rb is None:
            continue

        if fb is not None and rb is not None:
            stats["overlap"] += 1
            f_ok = fq >= min_q
            r_ok = rq >= min_q
            if not f_ok and not r_ok:
                cons.append("N")
                stats["n_calls"] += 1
            elif fb == rb:
                cons.append(fb)
                stats["agree"] += 1
            else:
                stats["disagree"] += 1
                if policy == "higher_q":
                    cons.append(fb if fq >= rq else rb)
                elif policy == "n":
                    cons.append("N")
                    stats["n_calls"] += 1
                else:  # iupac
                    cons.append(IUPAC.get(frozenset(fb + rb), "N"))
        elif fb is not None:
            stats["f_only"] += 1
            cons.append(fb if fq >= min_q else "N")
            if fq < min_q:
                stats["n_calls"] += 1
        else:
            stats["r_only"] += 1
            cons.append(rb if rq >= min_q else "N")
            if rq < min_q:
                stats["n_calls"] += 1

    stats["len"] = len(cons)
    return "".join(cons), stats


def _write_fasta(path: Path, seq_id: str, header_extra: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f">{seq_id} {header_extra}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--forward-ab1", required=True, type=Path)
    p.add_argument("--reverse-ab1", required=True, type=Path)
    p.add_argument("--out-fa",       required=True, type=Path)
    p.add_argument("--seq-id",       required=True)
    p.add_argument("--min-q",        type=int, default=20)
    p.add_argument("--disagreement", default="higher_q",
                   choices=["higher_q", "iupac", "n"])
    p.add_argument("--min-overlap",  type=int, default=20)
    p.add_argument("--match-score",    type=int, default=2)
    p.add_argument("--mismatch-score", type=int, default=-1)
    p.add_argument("--gap-open",       type=int, default=-4)
    p.add_argument("--gap-extend",     type=int, default=-1)
    p.add_argument("--report",         type=Path, default=None)
    p.add_argument("--trim-q",         type=int, default=0,
                   help="Mott trim quality cutoff applied to BOTH reads "
                        "before alignment (0 = disable). Trimmed flanks are "
                        "masked to N with Q=0.")
    p.add_argument("--trim-min-len",   type=int, default=50,
                   help="Skip trimming a read if it would leave fewer bases.")
    args = p.parse_args()

    f_seq, f_q = _read_ab1(args.forward_ab1)
    r_seq_raw, r_q_raw = _read_ab1(args.reverse_ab1)

    if args.trim_q > 0:
        f_seq, f_q, f_lo, f_hi = mask_trim(f_seq, f_q,
                                           q_cutoff=args.trim_q,
                                           min_len=args.trim_min_len)
        r_seq_raw, r_q_raw, r_lo, r_hi = mask_trim(r_seq_raw, r_q_raw,
                                                   q_cutoff=args.trim_q,
                                                   min_len=args.trim_min_len)
        print(f"[merge] {args.seq_id} trim_q={args.trim_q}  "
              f"fwd kept=[{f_lo},{f_hi})/{len(f_seq)}  "
              f"rev kept=[{r_lo},{r_hi})/{len(r_seq_raw)}", flush=True)

    # Try the two biologically meaningful orientations of the reverse read:
    #   1. revcomp  -> standard case: reverse trace reads the opposite strand,
    #                  so reverse-complementing puts both reads on the + strand.
    #   2. as-is    -> safety check in case the file labeled "reverse" was
    #                  actually sequenced from the same strand as the forward
    #                  (mis-labeled primer / wrong file naming).
    # Pick whichever orientation gives the higher global alignment score.
    r_seq_rc, r_q_rc = _rev_comp(r_seq_raw, r_q_raw)
    r_seq_fw, r_q_fw = r_seq_raw, r_q_raw

    aln_rc = _align(f_seq, r_seq_rc,
                    args.match_score, args.mismatch_score,
                    args.gap_open, args.gap_extend)
    aln_fw = _align(f_seq, r_seq_fw,
                    args.match_score, args.mismatch_score,
                    args.gap_open, args.gap_extend)

    if aln_rc.score >= aln_fw.score:
        aln, r_seq, r_q, orientation = aln_rc, r_seq_rc, r_q_rc, "revcomp"
    else:
        aln, r_seq, r_q, orientation = aln_fw, r_seq_fw, r_q_fw, "as_is"

    print(f"[merge] {args.seq_id} orientation={orientation} "
          f"score_revcomp={aln_rc.score:.1f} score_as_is={aln_fw.score:.1f}",
          flush=True)

    cons, stats = _merge(aln, f_seq, f_q, r_seq, r_q,
                         args.min_q, args.disagreement)
    stats["orientation"] = orientation
    stats["score_revcomp"] = float(aln_rc.score)
    stats["score_as_is"] = float(aln_fw.score)

    align_txt = args.out_fa.with_suffix(".fwd_rev_align.txt")
    align_txt.parent.mkdir(parents=True, exist_ok=True)
    rev_label = "Reverse_RC" if orientation == "revcomp" else "Reverse_asis"
    align_txt.write_text(
        f"# orientation={orientation}  "
        f"score_revcomp={aln_rc.score:.1f}  score_as_is={aln_fw.score:.1f}\n"
        + format_alignment("Forward", str(aln[0]), rev_label, str(aln[1]),
                           score=aln.score)
    )

    if stats["overlap"] < args.min_overlap:
        print(f"[merge] WARN {args.seq_id}: overlap={stats['overlap']} "
              f"< min_overlap={args.min_overlap}", flush=True)

    header = (f"len={stats['len']} overlap={stats['overlap']} "
              f"agree={stats['agree']} disagree={stats['disagree']} "
              f"f_only={stats['f_only']} r_only={stats['r_only']} "
              f"N={stats['n_calls']} min_q={args.min_q} policy={args.disagreement} "
              f"orientation={orientation}")
    _write_fasta(args.out_fa, args.seq_id, header, cons)

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            f"seq_id\t{args.seq_id}\n"
            f"forward\t{args.forward_ab1.name}\n"
            f"reverse\t{args.reverse_ab1.name}\n"
            f"length\t{stats['len']}\n"
            f"overlap\t{stats['overlap']}\n"
            f"agree\t{stats['agree']}\n"
            f"disagree\t{stats['disagree']}\n"
            f"forward_only\t{stats['f_only']}\n"
            f"reverse_only\t{stats['r_only']}\n"
            f"N_calls\t{stats['n_calls']}\n"
            f"min_q\t{args.min_q}\n"
            f"policy\t{args.disagreement}\n"
            f"orientation\t{orientation}\n"
            f"score_revcomp\t{aln_rc.score:.1f}\n"
            f"score_as_is\t{aln_fw.score:.1f}\n"
        )

    print(f"[merge] {args.seq_id} -> {args.out_fa.name}  "
          f"len={stats['len']} overlap={stats['overlap']} "
          f"disagree={stats['disagree']} N={stats['n_calls']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
