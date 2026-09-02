#!/usr/bin/env python3
"""Build F+R consensus from trace-reviewed FASTA/QUAL pairs.

For each sample directory under <in-dir>, pair every
    <stem>_forward.reviewed.fa  +  <stem>_forward.reviewed.qual
with
    <stem>_reverse.reviewed.fa  +  <stem>_reverse.reviewed.qual
reverse-complement the reverse read, pairwise-align the two strands, and
write a consensus FASTA next to the inputs:
    <stem>.consensus.fa
    <stem>.consensus.report.tsv

Consensus rules per aligned column (mirrors
``scripts/merge_consensus/01_pairwise_consensus.py``):
    * both gaps                       -> skip
    * one gap, base passes min-Q      -> take the base
    * one gap, base below min-Q       -> N
    * agreement (incl. matching IUPAC)-> take that base
    * disagreement                    -> resolve by --disagreement
        higher_q -> base from the strand with higher Phred Q
        iupac    -> IUPAC ambiguity code spanning both calls
        n        -> 'N'
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
format_alignment = importlib.import_module(
    "scripts.common.02_alignment_viz"
).format_alignment

IUPAC_PAIR = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}
IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _iupac_merge(b1: str, b2: str) -> str:
    bases = set(IUPAC_EXPAND.get(b1, "N")) | set(IUPAC_EXPAND.get(b2, "N"))
    if not bases or bases == {"N"}:
        return "N"
    return IUPAC_PAIR.get(frozenset(bases), "N") if len(bases) > 1 else next(iter(bases))


def _read_fa(path: Path) -> tuple[str, str]:
    rec = next(SeqIO.parse(str(path), "fasta"), None)
    if rec is None:
        raise SystemExit(f"ERROR: empty FASTA: {path}")
    return rec.id, str(rec.seq).upper()


def _read_qual(path: Path, expected_len: int) -> list[int]:
    quals = [int(x) for x in path.read_text().split()]
    if len(quals) < expected_len:
        quals = quals + [0] * (expected_len - len(quals))
    elif len(quals) > expected_len:
        quals = quals[:expected_len]
    return quals


def _rev_comp(seq: str, quals: list[int]) -> tuple[str, list[int]]:
    return str(Seq(seq).reverse_complement()), list(reversed(quals))


def _orient_pair(f_seq: str, r_seq_raw: str, r_q_raw: list[int],
                 match: int, mismatch: int, gap_open: int, gap_extend: int
                 ) -> tuple[str, list[int], str, float, float]:
    """Pick the orientation of the reverse read (as-is vs reverse-complement)
    that aligns better to the forward read. Returns (r_seq, r_q, orientation,
    score_as_is, score_revcomp). Robust to both raw-strand inputs and inputs
    that have already been oriented upstream (e.g. trace_review output).
    """
    a = PairwiseAligner()
    a.mode = "global"; a.match_score = match; a.mismatch_score = mismatch
    a.open_gap_score = gap_open; a.extend_gap_score = gap_extend
    score_asis = float(a.align(f_seq, r_seq_raw)[0].score)
    r_seq_rc, r_q_rc = _rev_comp(r_seq_raw, r_q_raw)
    score_rc = float(a.align(f_seq, r_seq_rc)[0].score)
    if score_rc >= score_asis:
        return r_seq_rc, r_q_rc, "revcomp", score_asis, score_rc
    return r_seq_raw, list(r_q_raw), "as_is", score_asis, score_rc


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
                # Treat IUPAC overlap (e.g. R vs A) as agreement on the
                # unambiguous base when one side is unambiguous.
                fset = set(IUPAC_EXPAND.get(fb, "N"))
                rset = set(IUPAC_EXPAND.get(rb, "N"))
                inter = fset & rset
                if inter and (len(fset) == 1 or len(rset) == 1):
                    cons.append(fb if len(fset) == 1 else rb)
                    stats["agree"] += 1
                    continue
                stats["disagree"] += 1
                if policy == "higher_q":
                    cons.append(fb if fq >= rq else rb)
                elif policy == "n":
                    cons.append("N")
                    stats["n_calls"] += 1
                else:  # iupac
                    cons.append(_iupac_merge(fb, rb))
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


def _pairs(sample_dir: Path,
           forward_suffix: str = "_forward",
           reverse_suffix: str = "_reverse"
           ) -> list[tuple[str, Path, Path, Path, Path]]:
    """Return list of (base_stem, fwd_fa, fwd_qual, rev_fa, rev_qual).

    Strand pairing strips ``forward_suffix`` / ``reverse_suffix`` from the
    stem (e.g. ``_forward`` / ``_reverse`` for the Macrogen-style default, or
    ``F`` / ``R`` for batches that name files ``2_2F.reviewed.fa`` and
    ``2_2R.reviewed.fa``). ``base_stem`` is the per-pair name with the strand
    suffix removed, used for the output ``<base>.consensus.fa`` next to the
    inputs.
    """
    fwd_pattern = f"*{forward_suffix}.reviewed.fa"
    fwd_tail = f"{forward_suffix}.reviewed.fa"
    out = []
    for fwd_fa in sorted(sample_dir.glob(fwd_pattern)):
        if not fwd_fa.name.endswith(fwd_tail):
            continue
        base = fwd_fa.name[: -len(fwd_tail)]
        fwd_qual = sample_dir / f"{base}{forward_suffix}.reviewed.qual"
        rev_fa   = sample_dir / f"{base}{reverse_suffix}.reviewed.fa"
        rev_qual = sample_dir / f"{base}{reverse_suffix}.reviewed.qual"
        if not (fwd_qual.exists() and rev_fa.exists() and rev_qual.exists()):
            print(f"  SKIP {base}: missing mate or .qual file")
            continue
        out.append((base, fwd_fa, fwd_qual, rev_fa, rev_qual))
    return out


def _process_pair(base: str, fwd_fa: Path, fwd_qual: Path,
                  rev_fa: Path, rev_qual: Path,
                  args) -> dict:
    f_id, f_seq = _read_fa(fwd_fa)
    r_id, r_seq_raw = _read_fa(rev_fa)
    f_q = _read_qual(fwd_qual, len(f_seq))
    r_q_raw = _read_qual(rev_qual, len(r_seq_raw))
    r_seq, r_q, orientation, score_asis, score_rc = _orient_pair(
        f_seq, r_seq_raw, r_q_raw,
        args.match_score, args.mismatch_score,
        args.gap_open, args.gap_extend)

    aln = _align(f_seq, r_seq,
                 args.match_score, args.mismatch_score,
                 args.gap_open, args.gap_extend)
    cons, stats = _merge(aln, f_seq, f_q, r_seq, r_q,
                         args.min_q, args.disagreement)

    if stats["overlap"] < args.min_overlap:
        print(f"  WARN {base}: overlap={stats['overlap']} < min_overlap={args.min_overlap}")

    out_fa = fwd_fa.parent / f"{base}.consensus.fa"
    header = (f"len={stats['len']} overlap={stats['overlap']} "
              f"agree={stats['agree']} disagree={stats['disagree']} "
              f"f_only={stats['f_only']} r_only={stats['r_only']} "
              f"N={stats['n_calls']} min_q={args.min_q} policy={args.disagreement} "
              f"orientation={orientation} source=trace_reviewed")
    _write_fasta(out_fa, base, header, cons)

    rev_label = "Reverse_RC" if orientation == "revcomp" else "Reverse_asis"
    (fwd_fa.parent / f"{base}.fwd_rev_align.txt").write_text(
        f"# orientation={orientation}  "
        f"score_as_is={score_asis:.1f}  score_revcomp={score_rc:.1f}\n"
        + format_alignment("Forward", str(aln[0]), rev_label, str(aln[1]),
                           score=aln.score)
    )

    report = fwd_fa.parent / f"{base}.consensus.report.tsv"
    report.write_text(
        f"seq_id\t{base}\n"
        f"forward\t{fwd_fa.name}\n"
        f"reverse\t{rev_fa.name}\n"
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
        f"score_as_is\t{score_asis:.1f}\n"
        f"score_revcomp\t{score_rc:.1f}\n"
    )

    print(f"  OK   {base}  len={stats['len']} overlap={stats['overlap']} "
          f"disagree={stats['disagree']} N={stats['n_calls']} "
          f"orient={orientation} -> {out_fa.name}")
    return {"base": base, **stats}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", required=True, type=Path,
                   help="03a_TraceReview directory (containing per-sample subdirs) "
                        "OR a single sample directory.")
    p.add_argument("--min-q",        type=int, default=20)
    p.add_argument("--disagreement", default="higher_q",
                   choices=["higher_q", "iupac", "n"])
    p.add_argument("--min-overlap",  type=int, default=20)
    p.add_argument("--match-score",    type=int, default=2)
    p.add_argument("--mismatch-score", type=int, default=-1)
    p.add_argument("--gap-open",       type=int, default=-4)
    p.add_argument("--gap-extend",     type=int, default=-1)
    p.add_argument("--forward-suffix", default="_forward",
                   help="Stem suffix marking forward-strand reviewed FASTAs "
                        "(default: _forward; use 'F' for batches that emit "
                        "files like '2_2F.reviewed.fa').")
    p.add_argument("--reverse-suffix", default="_reverse",
                   help="Stem suffix marking reverse-strand reviewed FASTAs "
                        "(default: _reverse; use 'R' for batches that emit "
                        "files like '2_2R.reviewed.fa').")
    args = p.parse_args()

    in_dir: Path = args.in_dir
    if not in_dir.is_dir():
        raise SystemExit(f"ERROR: not a directory: {in_dir}")

    # Decide whether in_dir is itself a sample dir or a parent of sample dirs.
    fwd_glob = f"*{args.forward_suffix}.reviewed.fa"
    has_pairs_here = any(in_dir.glob(fwd_glob))
    sample_dirs = [in_dir] if has_pairs_here else sorted(
        d for d in in_dir.iterdir()
        if d.is_dir() and any(d.glob(fwd_glob))
    )

    if not sample_dirs:
        print(f"No reviewed F/R pairs found under {in_dir}")
        return 0

    summary: list[dict] = []
    for sd in sample_dirs:
        print(f"[sample] {sd.name}")
        for pair in _pairs(sd, args.forward_suffix, args.reverse_suffix):
            summary.append({"sample": sd.name, **_process_pair(*pair, args)})

    if summary:
        out_summary = in_dir / "consensus_summary.tsv"
        cols = ["sample", "base", "len", "overlap", "agree", "disagree",
                "f_only", "r_only", "n_calls"]
        with out_summary.open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for row in summary:
                fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
        print(f"-> {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
