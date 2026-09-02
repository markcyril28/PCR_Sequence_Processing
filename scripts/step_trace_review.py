#!/usr/bin/env python3
"""Method 1: Reference-guided trace review / editing.

For each read assigned to a reference transcript:
  1. Read primary basecall + per-base Phred + per-base trace peaks from the .ab1
  2. Pairwise-align the basecall to the reference
  3. At every mismatch OR low-Q position (Q < min_quality), inspect the trace:
        - if the secondary peak / primary peak ratio at that scan >=
          secondary_peak_ratio AND one of the two peaks matches the reference
          allele, set the consensus to either:
              * the reference allele  (if allow_ref_overwrite=true and the ref
                allele is supported by the trace), OR
              * the IUPAC ambiguity code for the two competing peaks
        - otherwise leave the original call (or mark_uncalled if Q==0)
  4. Optionally trigger Tracy `decompose` for reads whose fraction of low-Q
     bases exceeds het_indel_trigger (likely heterozygous indel).

Outputs (per read) under <out-dir>/<sample>/:
    <stem>.reviewed.fa     edited per-strand sequence
    <stem>.reviewed.qual   per-base Phred (post-edit; edited bases get q=0 sentinel preserved)
    <stem>.edits.tsv       log of every base that was changed
    <stem>.ref_align.txt   Tracy-style visual alignment of read vs reference
    <stem>.tracy.vcf       (optional) variants from Tracy decompose
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alignment_viz import format_alignment
from trim_reads import mask_trim
from consensus_reviewed import _pairs as _consensus_pairs
from consensus_reviewed import _process_pair as _consensus_process_pair

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
    frozenset("ACG"): "V", frozenset("ACGT"): "N",
}


def _bool(v: str) -> bool:
    return str(v).lower() in {"1", "true", "yes", "y"}


def _read_assignments(path: Path) -> dict[tuple[str, str], str]:
    out = {}
    with path.open() as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sample, read, ref = parts[0], parts[1], parts[2]
            if ref in ("NO_HIT", "FAIL_THRESHOLD", "-", ""):
                continue
            out[(sample, read)] = ref
    return out


def _load_reference(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(path), "fasta")}


def _read_ab1(ab1: Path) -> tuple[str, list[int], list[int], dict[str, list[int]]]:
    """Return (seq, quals, peak_locations, traces) from an ABI file."""
    rec = next(SeqIO.parse(str(ab1), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in ABI file: {ab1}")
    ann = rec.annotations.get("abif_raw", {})

    seq = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality", []))
    if not quals:
        quals = [0] * len(seq)

    peaks_raw = ann.get("PLOC1") or ann.get("PLOC2") or b""
    if isinstance(peaks_raw, (bytes, bytearray)):
        # SHORT array, big-endian, 2 bytes per peak.
        import struct
        n = len(peaks_raw) // 2
        peaks = list(struct.unpack(f">{n}H", peaks_raw))
    else:
        peaks = list(peaks_raw)

    # Per-channel processed traces. ABIF stores them as DATA9..DATA12 in the
    # order of FWO_1 (the channel order), e.g. b"GATC".
    order = ann.get("FWO_1", b"GATC")
    if isinstance(order, bytes):
        order = order.decode("ascii", "ignore")
    chans = {}
    for ch, key in zip(order, ("DATA9", "DATA10", "DATA11", "DATA12")):
        data = ann.get(key, b"")
        if isinstance(data, (bytes, bytearray)):
            import struct
            n = len(data) // 2
            chans[ch] = list(struct.unpack(f">{n}h", data))
        else:
            chans[ch] = list(data)
    return seq, quals, peaks, chans


def _trace_peak_supports(chans: dict[str, list[int]], peaks: list[int],
                         pos: int, ref_base: str, ratio_thresh: float
                         ) -> tuple[bool, str | None, str | None]:
    """At basecall index `pos`, does the chromatogram support `ref_base`?
    Returns (supported, primary_base, secondary_base)."""
    if pos >= len(peaks) or ref_base not in "ACGT":
        return (False, None, None)
    scan = peaks[pos]
    intensities = {b: (chans[b][scan] if b in chans and scan < len(chans[b]) else 0)
                   for b in "ACGT"}
    ranked = sorted(intensities.items(), key=lambda kv: kv[1], reverse=True)
    primary, p_int = ranked[0]
    secondary, s_int = ranked[1]
    if p_int <= 0:
        return (False, None, None)
    ratio = s_int / p_int
    # Reference allele is supported if it's the primary OR a credible secondary.
    if primary == ref_base:
        return (True, primary, secondary if ratio >= ratio_thresh else None)
    if secondary == ref_base and ratio >= ratio_thresh:
        return (True, primary, secondary)
    return (False, primary, secondary if ratio >= ratio_thresh else None)


def _align(ref: str, read: str) -> tuple[str, str]:
    a = PairwiseAligner()
    a.mode = "global"; a.match_score = 2; a.mismatch_score = -1
    a.open_gap_score = -4; a.extend_gap_score = -1
    aln = a.align(ref, read)[0]
    return str(aln[0]), str(aln[1])


def _align_score(ref: str, read: str) -> float:
    a = PairwiseAligner()
    a.mode = "global"; a.match_score = 2; a.mismatch_score = -1
    a.open_gap_score = -4; a.extend_gap_score = -1
    return float(a.align(ref, read)[0].score)


_COMPLEMENT_CHAN = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _orient_to_reference(seq: str, quals: list[int], peaks: list[int],
                         chans: dict[str, list[int]], ref: str
                         ) -> tuple[str, list[int], list[int],
                                    dict[str, list[int]], str, float, float]:
    """Pick the orientation (as-is vs reverse-complement) that aligns better
    to ``ref``. Returns (seq, quals, peaks, chans, orientation, score_asis,
    score_rc). When orientation is ``revcomp``, the chromatogram channel dict
    is complement-remapped (A<->T, C<->G) so that downstream peak lookups read
    the correct strand at each scan, and the per-base ``peaks`` list is
    reversed so peaks[i] corresponds to the i-th base of the RC'd sequence.
    """
    score_asis = _align_score(ref, seq)
    rc_seq = str(Seq(seq).reverse_complement())
    score_rc = _align_score(ref, rc_seq)
    if score_rc > score_asis:
        rc_quals = list(reversed(quals))
        rc_peaks = list(reversed(peaks))
        rc_chans = {
            base: chans.get(_COMPLEMENT_CHAN[base], []) for base in "ACGT"
        }
        return rc_seq, rc_quals, rc_peaks, rc_chans, "revcomp", score_asis, score_rc
    return seq, quals, peaks, chans, "as_is", score_asis, score_rc


def _review_one(seq: str, quals: list[int], peaks: list[int],
                chans: dict[str, list[int]], ref: str,
                min_q: int, ratio_thresh: float,
                allow_overwrite: bool, require_support: bool,
                mark_uncalled: str
                ) -> tuple[str, list[int], list[dict], str, str]:
    aligned_ref, aligned_read = _align(ref, seq)
    new_bases: list[str] = []
    new_quals: list[int] = []
    edits: list[dict] = []
    read_idx = 0

    for ref_base, read_base in zip(aligned_ref, aligned_read):
        if ref_base == "-":
            # Insertion in read; keep as-is, no edit logged.
            if read_base != "-":
                new_bases.append(read_base)
                new_quals.append(quals[read_idx] if read_idx < len(quals) else 0)
                read_idx += 1
            continue
        if read_base == "-":
            # Deletion in read relative to ref. Look at the surrounding trace?
            # Conservative: leave a gap (will be merged later as N or ref).
            new_bases.append(mark_uncalled)
            new_quals.append(0)
            continue

        q = quals[read_idx] if read_idx < len(quals) else 0
        is_low_q  = q < min_q
        is_mismatch = (read_base != ref_base)

        if is_low_q or is_mismatch:
            supported, primary, secondary = _trace_peak_supports(
                chans, peaks, read_idx, ref_base, ratio_thresh)
            if supported and allow_overwrite:
                # Use IUPAC if both primary and ref differ but both peaks credible.
                if secondary and primary != ref_base:
                    code = IUPAC.get(frozenset(primary + ref_base), "N")
                    new_bases.append(code)
                else:
                    new_bases.append(ref_base)
                new_quals.append(max(q, min_q))  # bump quality post-edit
                edits.append({"pos": read_idx, "old": read_base, "new": new_bases[-1],
                              "ref": ref_base, "q": q, "primary": primary or "",
                              "secondary": secondary or "",
                              "reason": "low_q" if is_low_q else "mismatch"})
            elif require_support and is_low_q and read_base in "ACGT":
                # No trace support for ref AND no support for current call -> N
                if primary != read_base:
                    new_bases.append(mark_uncalled); new_quals.append(0)
                    edits.append({"pos": read_idx, "old": read_base, "new": mark_uncalled,
                                  "ref": ref_base, "q": q, "primary": primary or "",
                                  "secondary": secondary or "",
                                  "reason": "low_q_no_support"})
                else:
                    new_bases.append(read_base); new_quals.append(q)
            else:
                new_bases.append(read_base); new_quals.append(q)
        else:
            new_bases.append(read_base); new_quals.append(q)
        read_idx += 1

    return "".join(new_bases), new_quals, edits, aligned_ref, aligned_read


def _maybe_run_tracy(ab1: Path, ref_record_fa: Path, out_prefix: Path) -> Path | None:
    """Try Tracy decompose for het-indel reads. Returns the path to any
    successfully produced output file (json/bcf/decomp/...) or None.

    Tracy frequently fails to align very noisy traces with the default trim
    settings; on failure we retry once with no trimming before giving up
    quietly (no error message), since we already have the reviewed FASTA.
    """
    if shutil.which("tracy") is None:
        return None
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    attempts = [
        ["tracy", "decompose", "-r", str(ref_record_fa),
         "-o", str(out_prefix), str(ab1)],
        # Retry: disable quality/edge trimming (-q 0 -l 0 -t 0) so a usable
        # alignment region can still be found in low-Q reads.
        ["tracy", "decompose", "-r", str(ref_record_fa),
         "-o", str(out_prefix), "-q", "0", "-l", "0", "-t", "0", str(ab1)],
    ]
    last_err = ""
    for cmd in attempts:
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            last_err = ""
            break
        except subprocess.CalledProcessError as e:
            last_err = e.stderr.decode(errors="ignore")
    if last_err:
        # Soft-fail: trace_review already produced the reviewed FASTA, so a
        # tracy alignment failure is non-fatal: note it and move on.
        print(f"  note: tracy decompose unavailable for {ab1.name} "
              f"({last_err.strip().splitlines()[-1][:160] if last_err.strip() else 'failed'})")
        return None
    # Tracy writes several files prefixed with `out_prefix`; presence of any of
    # them (json/bcf/decomp/align*/abif) indicates success. Prefer json/bcf.
    candidates = sorted(out_prefix.parent.glob(out_prefix.name + "*"))
    for ext in (".json", ".bcf", ".vcf", ".decomp"):
        for c in candidates:
            if c.suffix == ext:
                return c
    return candidates[0] if candidates else None


def _process_one(sample: str, ab1_str: str, ref_id: str, ref_seq: str,
                 ref_fa_str: str, out_dir_str: str,
                 cfg: dict) -> tuple[str, dict]:
    ab1 = Path(ab1_str)
    stem = ab1.stem
    out_dir = Path(out_dir_str)
    out_sample = out_dir / sample
    out_sample.mkdir(parents=True, exist_ok=True)

    seq, quals, peaks, chans = _read_ab1(ab1)
    trim_tag = ""
    if cfg.get("trim_q", 0) > 0:
        seq, quals, t_lo, t_hi = mask_trim(
            seq, quals,
            q_cutoff=cfg["trim_q"],
            min_len=cfg.get("trim_min_len", 50),
        )
        trim_tag = f" trim=[{t_lo},{t_hi})"
    seq, quals, peaks, chans, orientation, score_asis, score_rc = \
        _orient_to_reference(seq, quals, peaks, chans, ref_seq)
    low_q_frac = sum(1 for q in quals if q < cfg["min_q"]) / max(len(quals), 1)

    if cfg["max_low_q_frac"] < 1.0 and low_q_frac > cfg["max_low_q_frac"]:
        row = {"sample": sample, "read": stem, "ref": ref_id, "edits": 0,
               "lowQ_frac": f"{low_q_frac:.3f}", "tracy_decompose": False,
               "status": "REJECTED_LOW_Q"}
        line = (f"  {sample}/{stem}  ref={ref_id}  REJECTED "
                f"(lowQ_frac={low_q_frac:.2f} > {cfg['max_low_q_frac']})")
        return line, row

    new_seq, new_quals, edits, aref_str, aread_str = _review_one(
        seq, quals, peaks, chans, ref_seq,
        cfg["min_q"], cfg["ratio"], cfg["allow_overwrite"],
        cfg["require_support"], cfg["mark_uncalled"])

    (out_sample / f"{stem}.reviewed.fa").write_text(
        f">{stem} ref={ref_id} reviewed=1 orientation={orientation}\n" +
        "\n".join(new_seq[i:i+60] for i in range(0, len(new_seq), 60)) + "\n")
    (out_sample / f"{stem}.reviewed.qual").write_text(
        " ".join(str(q) for q in new_quals) + "\n")
    with (out_sample / f"{stem}.edits.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pos","old","new","ref","q","primary","secondary","reason"], delimiter="\t")
        w.writeheader(); w.writerows(edits)
    (out_sample / f"{stem}.ref_align.txt").write_text(
        f"# orientation={orientation}  "
        f"score_as_is={score_asis:.1f}  score_revcomp={score_rc:.1f}\n"
        + format_alignment("PCR_Seq", aread_str, "Reference", aref_str))

    triggered_tracy = False
    if cfg["use_tracy"] and low_q_frac >= cfg["het_indel_trigger"]:
        vcf = _maybe_run_tracy(ab1, Path(ref_fa_str),
                               out_sample / f"{stem}.tracy")
        triggered_tracy = vcf is not None

    line = (f"  {sample}/{stem}  ref={ref_id}  edits={len(edits)} "
            f"lowQ_frac={low_q_frac:.2f} tracy={triggered_tracy} "
            f"orient={orientation}{trim_tag}")
    row = {"sample": sample, "read": stem, "ref": ref_id,
           "edits": len(edits), "lowQ_frac": f"{low_q_frac:.3f}",
           "tracy_decompose": triggered_tracy, "status": "OK"}
    return line, row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-root", required=True, type=Path)
    ap.add_argument("--reference",   required=True, type=Path)
    ap.add_argument("--assignments", required=True, type=Path)
    ap.add_argument("--out-dir",     required=True, type=Path)
    ap.add_argument("--min-quality",          type=int,   default=20)
    ap.add_argument("--secondary-peak-ratio", type=float, default=0.33)
    ap.add_argument("--allow-ref-overwrite",  default="true")
    ap.add_argument("--require-trace-support", default="true")
    ap.add_argument("--mark-uncalled",        default="N")
    ap.add_argument("--tracy-for-indels",     default="true")
    ap.add_argument("--het-indel-trigger",    type=float, default=0.30)
    ap.add_argument("--max-low-q-frac",       type=float, default=0.85,
                    help="Reject (do not emit reviewed FASTA) reads whose "
                         "low-Q fraction exceeds this threshold; 1.0 disables.")
    ap.add_argument("--threads",              type=int,   default=1,
                    help="Parallel worker processes (default: 1).")
    ap.add_argument("--trim-q",       type=int, default=0,
                    help="Mott trim quality cutoff applied before trace review "
                         "(0 = disable). Trimmed flanks are masked to N with "
                         "Q=0 so peak/trace indices remain valid.")
    ap.add_argument("--trim-min-len", type=int, default=50,
                    help="Skip trimming a read if it would leave fewer bases.")
    ap.add_argument("--include", action="append", default=None,
                    help="Restrict to these sample subdir names (repeatable).")
    ap.add_argument("--forward-suffix", default="_forward",
                    help="Stem suffix marking forward-strand reads (default: _forward).")
    ap.add_argument("--reverse-suffix", default="_reverse",
                    help="Stem suffix marking reverse-strand reads (default: _reverse).")
    ap.add_argument("--emit-consensus", default="true",
                    help="After review, write <base>.consensus.fa next to each "
                         "F/R reviewed FASTA pair (default: true).")
    ap.add_argument("--consensus-min-q",          type=int, default=5,
                    help="Per-strand Phred floor at F+R merge (default 5; "
                         "matches trace_review.min_quality).")
    ap.add_argument("--consensus-disagreement",   default="iupac",
                    choices=["higher_q", "iupac", "n"],
                    help="F/R conflict policy when both strands disagree "
                         "(default: iupac, matches pipeline STEP 4).")
    ap.add_argument("--consensus-min-overlap",    type=int, default=20)
    ap.add_argument("--consensus-match-score",    type=int, default=2)
    ap.add_argument("--consensus-mismatch-score", type=int, default=-1)
    ap.add_argument("--consensus-gap-open",       type=int, default=-4)
    ap.add_argument("--consensus-gap-extend",     type=int, default=-1)
    args = ap.parse_args()
    _include = set(args.include) if args.include else None

    refs = _load_reference(args.reference)
    assignments = _read_assignments(args.assignments)
    cfg = {
        "min_q": args.min_quality,
        "ratio": args.secondary_peak_ratio,
        "allow_overwrite": _bool(args.allow_ref_overwrite),
        "require_support": _bool(args.require_trace_support),
        "mark_uncalled": args.mark_uncalled,
        "use_tracy": _bool(args.tracy_for_indels),
        "het_indel_trigger": args.het_indel_trigger,
        "max_low_q_frac": args.max_low_q_frac,
        "trim_q": args.trim_q,
        "trim_min_len": args.trim_min_len,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Pre-materialize per-reference single-record FASTAs (Tracy needs them).
    ref_fa_cache: dict[str, Path] = {}
    def _ref_fa(rid: str) -> Path:
        if rid in ref_fa_cache:
            return ref_fa_cache[rid]
        p = args.out_dir / "_ref" / f"{rid}.fa"
        p.parent.mkdir(parents=True, exist_ok=True)
        seq = refs[rid]
        p.write_text(f">{rid}\n" + "\n".join(seq[i:i+60] for i in range(0, len(seq), 60)) + "\n")
        ref_fa_cache[rid] = p
        return p

    # Build the task list, resolving ref assignments + mate fallback up front.
    tasks: list[tuple[str, str, str, str, str]] = []
    for sample_dir in sorted(p for p in args.sample_root.iterdir() if p.is_dir()):
        if _include is not None and sample_dir.name not in _include:
            continue
        for ab1 in sorted(sample_dir.glob("*.ab1")):
            stem = ab1.stem
            key = (sample_dir.name, stem)
            ref_id = assignments.get(key)
            if ref_id is None:
                if args.forward_suffix in stem:
                    mate_stem = stem.replace(args.forward_suffix, args.reverse_suffix)
                elif args.reverse_suffix in stem:
                    mate_stem = stem.replace(args.reverse_suffix, args.forward_suffix)
                else:
                    continue
                ref_id = assignments.get((sample_dir.name, mate_stem))
                if ref_id is None:
                    continue
            if ref_id not in refs:
                print(f"  WARN: ref {ref_id} missing from FASTA, skipping {stem}", flush=True)
                continue
            ref_fa_path = _ref_fa(ref_id)
            tasks.append((sample_dir.name, str(ab1), ref_id,
                          refs[ref_id], str(ref_fa_path)))

    threads = max(1, int(args.threads))
    out_dir_str = str(args.out_dir)
    summary_rows: list[dict] = []

    if threads <= 1 or len(tasks) <= 1:
        print(f"[trace_review] dispatching {len(tasks)} reads sequentially", flush=True)
        for sample, ab1_str, ref_id, ref_seq, ref_fa_str in tasks:
            line, row = _process_one(sample, ab1_str, ref_id, ref_seq,
                                     ref_fa_str, out_dir_str, cfg)
            print(line, flush=True)
            summary_rows.append(row)
    else:
        print(f"[trace_review] dispatching {len(tasks)} reads across {threads} workers", flush=True)
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futs = [pool.submit(_process_one, sample, ab1_str, ref_id, ref_seq,
                                ref_fa_str, out_dir_str, cfg)
                    for sample, ab1_str, ref_id, ref_seq, ref_fa_str in tasks]
            for fut in as_completed(futs):
                line, row = fut.result()
                print(line, flush=True)
                summary_rows.append(row)

    with (args.out_dir / "trace_review_summary.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample","read","ref","edits","lowQ_frac","tracy_decompose","status"], delimiter="\t")
        w.writeheader(); w.writerows(summary_rows)
    print(f"  -> {args.out_dir / 'trace_review_summary.tsv'}", flush=True)

    if _bool(args.emit_consensus):
        _emit_pair_consensus(args)
    return 0


def _emit_pair_consensus(args) -> None:
    """For each per-sample subdir under out-dir, pair every reviewed F/R FASTA
    and write <base>.consensus.fa alongside the inputs.

    Reuses ``consensus_reviewed._pairs`` and ``_process_pair`` so the merge
    logic stays in one place and matches the standalone tool's behavior.
    """
    cons_args = argparse.Namespace(
        min_q=args.consensus_min_q,
        disagreement=args.consensus_disagreement,
        min_overlap=args.consensus_min_overlap,
        match_score=args.consensus_match_score,
        mismatch_score=args.consensus_mismatch_score,
        gap_open=args.consensus_gap_open,
        gap_extend=args.consensus_gap_extend,
    )

    summary: list[dict] = []
    sample_dirs = sorted(d for d in args.out_dir.iterdir()
                         if d.is_dir() and not d.name.startswith("_"))
    for sd in sample_dirs:
        pairs = _consensus_pairs(sd, args.forward_suffix, args.reverse_suffix)
        if not pairs:
            continue
        print(f"[consensus] {sd.name}", flush=True)
        for pair in pairs:
            summary.append({"sample": sd.name, **_consensus_process_pair(*pair, cons_args)})

    if summary:
        out_summary = args.out_dir / "consensus_summary.tsv"
        cols = ["sample", "base", "len", "overlap", "agree", "disagree",
                "f_only", "r_only", "n_calls"]
        with out_summary.open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for row in summary:
                fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
        print(f"  -> {out_summary}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
