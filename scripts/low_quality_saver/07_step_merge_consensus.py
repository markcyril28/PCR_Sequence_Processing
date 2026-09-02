#!/usr/bin/env python3
"""Step 4: merge per-strand calls into one consensus FASTA per amplicon.

Discovers forward + reverse pairs in <call-dir>/<sample>/. Pair files are
expected to be named <stem>.<ext> where ext is either `.reviewed.fa` (Method 1)
or `.aligned.fa` (Method 2), both accompanied by `<stem>.<ext base>.qual`.

Pairs are matched by stripping `_forward` / `_reverse` from the stem.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}

EXT_PATTERNS = [(".reviewed.fa", ".reviewed.qual"),
                (".aligned.fa",  ".aligned.qual")]


def _read_assignments(path: Path) -> dict[tuple[str, str], str]:
    out = {}
    with path.open() as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2] not in ("NO_HIT", "FAIL_THRESHOLD", "-", ""):
                out[(p[0], p[1])] = p[2]
    return out


def _read_seq(p: Path) -> str:
    rec = next(SeqIO.parse(str(p), "fasta"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in FASTA file: {p}")
    return str(rec.seq).upper()


def _read_qual(p: Path, n: int) -> list[int]:
    if not p.exists():
        return [0] * n
    qs = [int(x) for x in p.read_text().split()]
    if len(qs) > n:
        print(f"  WARN: quality file longer than sequence ({len(qs)} > {n}): {p}")
    if len(qs) < n:
        qs += [0] * (n - len(qs))
    return qs[:n]


def _merge(ref: str, fb: list[str], fq: list[int], rb: list[str], rq: list[int],
           min_q: int, policy: str, min_cov: int) -> str:
    out = []
    L = len(ref)
    for i in range(L):
        a = fb[i] if i < len(fb) else "-"
        b = rb[i] if i < len(rb) else "-"
        qa = fq[i] if i < len(fq) else 0
        qb = rq[i] if i < len(rq) else 0
        a_ok = a in "ACGTRYSWKMBDHV" and qa >= min_q
        b_ok = b in "ACGTRYSWKMBDHV" and qb >= min_q
        if not a_ok and not b_ok:
            out.append("N"); continue
        if a_ok and not b_ok:
            out.append(a if min_cov == 1 else "N"); continue
        if b_ok and not a_ok:
            out.append(b if min_cov == 1 else "N"); continue
        if a == b:
            out.append(a)
        elif policy == "higher_q":
            out.append(a if qa >= qb else b)
        elif policy == "n":
            out.append("N")
        else:
            out.append(IUPAC.get(frozenset(a + b), "N"))
    return "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference",   required=True, type=Path)
    ap.add_argument("--call-dir",    required=True, type=Path)
    ap.add_argument("--assignments", required=True, type=Path)
    ap.add_argument("--out-dir",     required=True, type=Path)
    ap.add_argument("--min-q",   type=int, default=20)
    ap.add_argument("--policy",  default="iupac")
    ap.add_argument("--min-cov", type=int, default=1)
    ap.add_argument("--forward-suffix", default="_forward",
                    help="Stem suffix marking forward-strand reads (default: _forward).")
    ap.add_argument("--reverse-suffix", default="_reverse",
                    help="Stem suffix marking reverse-strand reads (default: _reverse).")
    args = ap.parse_args()
    _strand_re = re.compile(rf"(?:{re.escape(args.forward_suffix)}|{re.escape(args.reverse_suffix)})$")

    refs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(args.reference), "fasta")}
    assignments = _read_assignments(args.assignments)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Group by (sample, base_stem) where base_stem strips _forward/_reverse.
    pairs: dict[tuple[str, str], dict[str, Path]] = {}
    for sample_dir in args.call_dir.iterdir():
        if not sample_dir.is_dir() or sample_dir.name.startswith("_"):
            continue
        for fa in sample_dir.iterdir():
            for fa_ext, q_ext in EXT_PATTERNS:
                if fa.name.endswith(fa_ext):
                    stem = fa.name[: -len(fa_ext)]
                    base = _strand_re.sub("", stem)
                    strand = "F" if stem.endswith(args.forward_suffix) else "R"
                    key = (sample_dir.name, base)
                    pairs.setdefault(key, {})[strand] = fa
                    pairs[key][f"{strand}_q"] = sample_dir / f"{stem}{q_ext}"
                    pairs[key]["sample"] = sample_dir.name
                    pairs[key]["full_stem_" + strand] = stem
                    break

    for (sample, base), info in sorted(pairs.items()):
        have_f = "F" in info
        have_r = "R" in info
        if not (have_f or have_r):
            continue
        if not (have_f and have_r) and args.min_cov >= 2:
            print(f"  SKIP {sample}/{base} (need both forward and reverse for min_cov>=2)")
            continue
        # Resolve which reference to use (prefer the forward read's assignment).
        ref_id = None
        if have_f:
            ref_id = assignments.get((sample, info["full_stem_F"]))
        if not ref_id and have_r:
            ref_id = assignments.get((sample, info["full_stem_R"]))
        if not ref_id or ref_id not in refs:
            print(f"  SKIP {sample}/{base} (no reference assignment)")
            continue

        ref = refs[ref_id]
        ext_is_reviewed = ((info["F"].name.endswith(".reviewed.fa") if have_f
                            else info["R"].name.endswith(".reviewed.fa")))

        def _load_strand(fa_key: str, q_key: str, is_reverse: bool):
            if fa_key not in info:
                return ["-"] * len(ref), [0] * len(ref)
            seq = _read_seq(info[fa_key])
            qs = _read_qual(info[q_key], len(seq))
            # Method 1's reviewed seqs are read-coordinates and may be on
            # either strand depending on what trace_review picked. Score-test
            # both orientations against the reference and use the better one.
            # This is robust to upstream orientation policy and to mislabeled
            # _forward/_reverse files.
            if ext_is_reviewed:
                from Bio.Align import PairwiseAligner
                a = PairwiseAligner(); a.mode = "global"
                a.match_score = 2; a.mismatch_score = -1
                a.open_gap_score = -4; a.extend_gap_score = -1
                score_asis = float(a.align(ref, seq)[0].score)
                rc_seq = str(Seq(seq).reverse_complement())
                score_rc = float(a.align(ref, rc_seq)[0].score)
                if score_rc > score_asis:
                    seq = rc_seq
                    qs = list(reversed(qs))
                aln = a.align(ref, seq)[0]
                ar, aq = str(aln[0]), str(aln[1])
                bs, qq, ri = [], [], 0
                for r, c in zip(ar, aq):
                    if r == "-":
                        if c != "-":
                            ri += 1
                        continue
                    if c == "-":
                        bs.append("-"); qq.append(0)
                    else:
                        bs.append(c); qq.append(qs[ri] if ri < len(qs) else 0)
                        ri += 1
                return bs, qq
            # consensus_assembly outputs are already on reference coordinates.
            return list(seq), qs

        fb, fqp = _load_strand("F", "F_q", is_reverse=False)
        rb, rqp = _load_strand("R", "R_q", is_reverse=True)

        cons = _merge(ref, fb, fqp, rb, rqp,
                      args.min_q, args.policy, args.min_cov)
        cov_note = ("both" if (have_f and have_r)
                    else ("forward_only" if have_f else "reverse_only"))
        out_fa = args.out_dir / sample / f"{base}_consensus.fa"
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        with out_fa.open("w") as fh:
            fh.write(f">{sample}|{base} ref={ref_id} len={len(cons)} "
                     f"min_q={args.min_q} policy={args.policy} "
                     f"min_cov={args.min_cov} strands={cov_note}\n")
            for i in range(0, len(cons), 60):
                fh.write(cons[i:i+60] + "\n")
        n = cons.count("N")
        print(f"  {sample}/{base}  ref={ref_id}  len={len(cons)}  "
              f"N={n} ({100*n/max(len(cons),1):.1f}%)  strands={cov_note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
