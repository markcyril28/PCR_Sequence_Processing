#!/usr/bin/env python3
"""Method 3: Tracy reference-guided basecalling / decomposition.

For each assigned read, runs `tracy align` (single trace; clean reads) or
`tracy decompose` (heterozygous indel deconvolution). Writes Tracy outputs
under <out-dir>/<sample>/<stem>.* . When mode="auto", `decompose` is used if
the .ab1's mean Phred quality is below 25, otherwise `align`.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO


def _bool(v: str) -> bool: return str(v).lower() in {"1", "true", "yes", "y"}


def _read_assignments(path: Path) -> dict[tuple[str, str], str]:
    out = {}
    with path.open() as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2] not in ("NO_HIT", "FAIL_THRESHOLD", "-", ""):
                out[(p[0], p[1])] = p[2]
    return out


def _mean_q(ab1: Path) -> float:
    rec = next(SeqIO.parse(str(ab1), "abi"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no records in ABI file: {ab1}")
    qs = list(rec.letter_annotations.get("phred_quality", []))
    return (sum(qs) / len(qs)) if qs else 0.0


def _ab1_quals(ab1: Path) -> list[int]:
    rec = next(SeqIO.parse(str(ab1), "abi"), None)
    if rec is None:
        return []
    return list(rec.letter_annotations.get("phred_quality", []))


def _low_q_frac(ab1: Path, min_q: int) -> float:
    qs = _ab1_quals(ab1)
    return (sum(1 for q in qs if q < min_q) / len(qs)) if qs else 1.0


def _try_tracy(cmd_base: list[str], ab1: Path, prefix: Path,
               trim_left: int, trim_right: int, trim_quality: int,
               kmer: int | None) -> tuple[bool, str]:
    """Run tracy with the given trim/kmer params. Returns (ok, stderr)."""
    cmd = list(cmd_base) + [
        "-o", str(prefix),
        "-q", str(trim_quality),
        "-l", str(trim_left),
        "-t", str(trim_right),
    ]
    if kmer is not None:
        cmd += ["-k", str(kmer)]
    cmd.append(str(ab1))
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        return False, e.stderr.decode(errors="ignore")


def _process_one(sample: str, ab1: Path, ref_id: str, cfg: dict) -> str:
    """Worker: run tracy on a single ab1, with retries/fallbacks. Returns log line."""
    out_dir: Path = cfg["out_dir"]
    ref_fa = cfg["ref_dir"] / f"{ref_id}.fa"

    choose = cfg["mode"]
    if choose == "auto":
        choose = "decompose" if _mean_q(ab1) < cfg["auto_decompose_q"] else "align"

    if cfg["max_low_q_frac"] < 1.0:
        lqf = _low_q_frac(ab1, cfg["trim_quality"])
        if lqf > cfg["max_low_q_frac"]:
            return (f"  SKIP {sample}/{ab1.stem}  ref={ref_id}  "
                    f"lowQ_frac={lqf:.2f} > {cfg['max_low_q_frac']}")

    out_sample = out_dir / sample
    out_sample.mkdir(parents=True, exist_ok=True)
    stem = ab1.stem
    prefix = out_sample / f"{stem}.tracy_{choose}"

    def _base(action: str) -> list[str]:
        b = ["tracy", action, "-r", str(ref_fa), "-p", str(cfg["pratio"])]
        if not cfg["call_variants"]:
            b += ["-v", "0"]
        return b

    kmer = cfg["kmer"] if choose == "decompose" else None
    attempts: list[tuple[str, int, int, int | None, Path]] = [
        (choose, cfg["trim_left"], cfg["trim_right"], kmer, prefix),
        (choose, max(cfg["trim_left"] // 2, 0), max(cfg["trim_right"] // 2, 0),
         kmer, prefix.with_name(prefix.name + ".retry")),
        (choose, 0, 0, kmer, prefix.with_name(prefix.name + ".notrim")),
    ]
    if choose == "align":
        attempts.append(("decompose", 0, 0, cfg["kmer"],
                         out_sample / f"{stem}.tracy_decompose.fallback"))
    else:
        attempts.append(("align", 0, 0, None,
                         out_sample / f"{stem}.tracy_align.fallback"))

    ok = False
    last_err = ""
    used_action = choose
    for action, tl, tr, km, pfx in attempts:
        ok, last_err = _try_tracy(_base(action), ab1, pfx, tl, tr,
                                  cfg["trim_quality"], km)
        if ok:
            used_action = action
            break

    if ok:
        tag = used_action if used_action == choose else f"{used_action} (fallback)"
        # Replace Tracy's generic >Alt / Alt labels with >PCR_Seq so output
        # files use project-specific naming rather than Tracy's internal term.
        new_hdr = f">PCR_Seq sample={sample} read={stem} ref={ref_id}"
        for txt_file in pfx.parent.glob(pfx.name + "*.txt"):
            content = txt_file.read_text()
            if ">Alt\n" in content:
                txt_file.write_text(content.replace(">Alt\n", new_hdr + "\n", 1))
        for align_file in pfx.parent.glob(pfx.name + "*.align*"):
            content = align_file.read_text()
            if "Alt" in content:
                # FASTA headers: >Alt1 / >Alt2 -> >PCR_Seq1 / >PCR_Seq2
                new_content = re.sub(r'>Alt(\d+)', r'>PCR_Seq\1', content)
                # Alignment block labels at line start: "Alt " -> "PCR_Seq "
                new_content = re.sub(r'^Alt\b', 'PCR_Seq', new_content, flags=re.MULTILINE)
                if new_content != content:
                    align_file.write_text(new_content)
        return f"  {sample}/{stem}  -> {tag}"
    return f"  ERR {sample}/{stem}: {last_err[:300]}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-root", required=True, type=Path)
    ap.add_argument("--reference",   required=True, type=Path)
    ap.add_argument("--assignments", required=True, type=Path)
    ap.add_argument("--out-dir",     required=True, type=Path)
    ap.add_argument("--mode",         default="auto")  # align | decompose | auto
    ap.add_argument("--pratio",       type=float, default=0.33)
    ap.add_argument("--trim-left",    type=int,   default=50)
    ap.add_argument("--trim-right",   type=int,   default=50)
    ap.add_argument("--trim-quality", type=int,   default=20)
    ap.add_argument("--call-variants", default="true")
    ap.add_argument("--kmer",         type=int,   default=15)
    ap.add_argument("--auto-decompose-q", type=float, default=25.0,
                    help="In mode=auto, use 'decompose' if mean Phred Q < this; else 'align'.")
    ap.add_argument("--max-low-q-frac", type=float, default=0.85,
                    help="Skip reads whose fraction of bases below trim-quality "
                         "exceeds this threshold (1.0 disables the filter).")
    ap.add_argument("--threads", type=int, default=1,
                    help="Number of parallel tracy invocations. tracy itself is "
                         "single-threaded per trace; reads are processed in parallel.")
    ap.add_argument("--include", action="append", default=None,
                    help="Restrict to these sample subdir names (repeatable).")
    ap.add_argument("--forward-suffix", default="_forward",
                    help="Stem suffix marking forward-strand reads (default: _forward).")
    ap.add_argument("--reverse-suffix", default="_reverse",
                    help="Stem suffix marking reverse-strand reads (default: _reverse).")
    args = ap.parse_args()
    _include = set(args.include) if args.include else None

    if shutil.which("tracy") is None:
        raise SystemExit("tracy not found on PATH (run setup.sh)")

    refs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(args.reference), "fasta")}
    assignments = _read_assignments(args.assignments)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ref_cache: dict[str, Path] = {}

    def _ref_fa(rid: str) -> Path:
        if rid in ref_cache: return ref_cache[rid]
        p = args.out_dir / "_ref" / f"{rid}.fa"
        p.parent.mkdir(parents=True, exist_ok=True)
        seq = refs[rid]
        p.write_text(f">{rid}\n" + "\n".join(seq[i:i+60] for i in range(0, len(seq), 60)) + "\n")
        ref_cache[rid] = p
        return p

    tasks: list[tuple[str, Path, str]] = []
    for sample_dir in sorted(p for p in args.sample_root.iterdir() if p.is_dir()):
        if _include is not None and sample_dir.name not in _include:
            continue
        for ab1 in sorted(sample_dir.glob("*.ab1")):
            stem = ab1.stem
            key = (sample_dir.name, stem)
            ref_id = assignments.get(key)
            if ref_id is None:
                # Reverse reads are not BLAST-assigned directly; fall back to mate.
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
                continue
            # Pre-materialize the per-ref FASTA in the parent so workers
            # never race writing the same file.
            _ref_fa(ref_id)
            tasks.append((sample_dir.name, ab1, ref_id))

    threads = max(1, int(args.threads))
    common = dict(
        out_dir=args.out_dir,
        ref_dir=args.out_dir / "_ref",
        mode=args.mode,
        pratio=args.pratio,
        trim_left=args.trim_left,
        trim_right=args.trim_right,
        trim_quality=args.trim_quality,
        call_variants=_bool(args.call_variants),
        kmer=args.kmer,
        auto_decompose_q=args.auto_decompose_q,
        max_low_q_frac=args.max_low_q_frac,
    )

    if threads == 1 or len(tasks) <= 1:
        for sample, ab1, ref_id in tasks:
            print(_process_one(sample, ab1, ref_id, common), flush=True)
    else:
        print(f"[tracy_decompose] dispatching {len(tasks)} reads across {threads} workers",
              flush=True)
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futs = [pool.submit(_process_one, s, a, r, common) for s, a, r in tasks]
            for fut in as_completed(futs):
                print(fut.result(), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
