#!/usr/bin/env python3
"""Build one final consensus from overlapping inside/outside primer products.

The low-quality saver pipeline can review and merge normal two-read pairs.
For I_HC01109784, those pairs are intermediate evidence only:

    PH_2_DMP8-OF + PH_2_DMP8-IR -> PH_2_DMP8 consensus
    PH_3_DMP8-IF + PH_3_DMP8-OR -> PH_3_DMP8 consensus

This script merges those overlapping consensuses into one final DNA sequence.
It writes both a reference-length consensus and a covered-span consensus.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore

from Bio import SeqIO
from Bio.Align import PairwiseAligner

IUPAC_EXPAND = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}, "U": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": set(), "-": set(),
}
IUPAC_COLLAPSE = {
    frozenset({"A"}): "A", frozenset({"C"}): "C",
    frozenset({"G"}): "G", frozenset({"T"}): "T",
    frozenset({"A", "G"}): "R", frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S", frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K", frozenset({"A", "C"}): "M",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "G", "T"}): "N",
}


def _read_one_fasta(path: Path) -> tuple[str, str, str]:
    rec = next(SeqIO.parse(str(path), "fasta"), None)
    if rec is None:
        raise SystemExit(f"ERROR: no FASTA record in {path}")
    return rec.id, str(rec.seq).upper(), rec.description


def _ref_from_description(description: str) -> str | None:
    m = re.search(r"\bref=([^\s]+)", description)
    return m.group(1) if m else None


def _load_refs(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(path), "fasta")}


def _project_to_ref(seq: str, ref: str) -> list[str]:
    """Return one base per reference position, using N for uncovered gaps."""
    if len(seq) == len(ref):
        return list(seq)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    aln = aligner.align(ref, seq)[0]
    aref, aseq = str(aln[0]), str(aln[1])
    out: list[str] = []
    for r_base, q_base in zip(aref, aseq):
        if r_base == "-":
            continue
        out.append("N" if q_base == "-" else q_base.upper())
    if len(out) != len(ref):
        raise SystemExit(
            f"ERROR: projected sequence length {len(out)} does not match "
            f"reference length {len(ref)}"
        )
    return out


def _collapse(bases: set[str]) -> str:
    return IUPAC_COLLAPSE.get(frozenset(bases), "N")


def _merge_columns(projected: list[tuple[str, list[str]]]) -> tuple[str, list[dict]]:
    if not projected:
        raise SystemExit("ERROR: no input consensuses provided")
    length = len(projected[0][1])
    for label, seq in projected:
        if len(seq) != length:
            raise SystemExit(f"ERROR: projected length mismatch for {label}")

    out: list[str] = []
    rows: list[dict] = []
    for idx in range(length):
        evid = []
        for label, seq in projected:
            base = seq[idx].upper()
            expanded = IUPAC_EXPAND.get(base, set())
            if expanded:
                evid.append((label, base, expanded))

        if not evid:
            out.append("N")
            rows.append({
                "pos_1based": idx + 1,
                "consensus": "N",
                "coverage": 0,
                "status": "uncovered",
                "evidence": "",
            })
            continue

        sets = [e[2] for e in evid]
        shared = set.intersection(*sets) if len(sets) > 1 else set(sets[0])
        status = "single_source"
        if len(evid) > 1 and shared:
            status = "agree"
            merged = _collapse(shared)
        elif len(evid) > 1:
            status = "conflict_iupac"
            merged = _collapse(set.union(*sets))
        else:
            merged = _collapse(sets[0])

        out.append(merged)
        rows.append({
            "pos_1based": idx + 1,
            "consensus": merged,
            "coverage": len(evid),
            "status": status,
            "evidence": ";".join(f"{label}:{base}" for label, base, _ in evid),
        })
    return "".join(out), rows


def _write_fasta(path: Path, seq_id: str, description: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f">{seq_id} {description}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def _load_config(workspace: Path, config: Path) -> tuple[Path, Path, str, list[Path]]:
    cfg_path = config if config.is_absolute() else workspace / config
    data = tomllib.loads(cfg_path.read_text(encoding="utf-8"))
    project = data.get("project", {}) or {}
    overlap = data.get("overlap_consensus", {}) or {}

    reference = Path(project.get("reference", ""))
    if not reference:
        raise SystemExit("ERROR: [project].reference is required")
    reference = reference if reference.is_absolute() else workspace / reference

    output_dir = Path(overlap.get("output_dir", ""))
    if not output_dir:
        raise SystemExit("ERROR: [overlap_consensus].output_dir is required")
    output_dir = output_dir if output_dir.is_absolute() else workspace / output_dir

    seq_id = str(overlap.get("seq_id", "overlap_consensus"))
    input_fastas = overlap.get("input_fastas") or []
    if not isinstance(input_fastas, list) or not input_fastas:
        raise SystemExit("ERROR: [overlap_consensus].input_fastas must be a non-empty array")
    inputs = [Path(p) if Path(p).is_absolute() else workspace / str(p) for p in input_fastas]
    return reference, output_dir, seq_id, inputs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", type=Path, default=Path.cwd())
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args()

    workspace = args.workspace.resolve()
    reference, output_dir, seq_id, input_fastas = _load_config(workspace, args.config)
    refs = _load_refs(reference)
    if not refs:
        raise SystemExit(f"ERROR: no records in reference FASTA: {reference}")

    loaded = []
    ref_ids = []
    for fasta in input_fastas:
        if not fasta.is_file():
            raise SystemExit(f"ERROR: missing overlap input FASTA: {fasta}")
        rec_id, seq, desc = _read_one_fasta(fasta)
        ref_id = _ref_from_description(desc)
        if ref_id is None:
            raise SystemExit(f"ERROR: no ref= field in FASTA header: {fasta}")
        if ref_id not in refs:
            raise SystemExit(f"ERROR: ref '{ref_id}' from {fasta} is not in {reference}")
        loaded.append((fasta, rec_id, seq, ref_id))
        ref_ids.append(ref_id)

    unique_refs = sorted(set(ref_ids))
    if len(unique_refs) != 1:
        raise SystemExit(
            "ERROR: overlap inputs were assigned to different references: "
            + ", ".join(unique_refs)
        )
    ref_id = unique_refs[0]
    ref_seq = refs[ref_id]

    projected = []
    source_rows = []
    for fasta, rec_id, seq, _ref_id in loaded:
        ref_bases = _project_to_ref(seq, ref_seq)
        covered = [i for i, b in enumerate(ref_bases) if IUPAC_EXPAND.get(b, set())]
        start = covered[0] + 1 if covered else 0
        end = covered[-1] + 1 if covered else 0
        projected.append((rec_id, ref_bases))
        source_rows.append({
            "source_id": rec_id,
            "source_fasta": str(fasta.relative_to(workspace) if fasta.is_relative_to(workspace) else fasta),
            "ref_id": ref_id,
            "length": len(seq),
            "covered_bases": len(covered),
            "covered_start_1based": start,
            "covered_end_1based": end,
        })

    merged, per_pos = _merge_columns(projected)
    covered_final = [i for i, b in enumerate(merged) if IUPAC_EXPAND.get(b, set())]
    start = covered_final[0] + 1 if covered_final else 0
    end = covered_final[-1] + 1 if covered_final else 0
    covered_span = merged[start - 1:end] if covered_final else ""

    output_dir.mkdir(parents=True, exist_ok=True)
    ref_len_fa = output_dir / f"{seq_id}_reference_length.fa"
    span_fa = output_dir / f"{seq_id}_covered_span.fa"
    report_tsv = output_dir / f"{seq_id}_overlap_report.tsv"
    per_pos_tsv = output_dir / f"{seq_id}_per_position.tsv"

    conflict_count = sum(1 for row in per_pos if row["status"] == "conflict_iupac")
    overlap_count = sum(1 for row in per_pos if row["coverage"] >= 2)
    agree_count = sum(1 for row in per_pos if row["status"] == "agree")
    n_count = merged.count("N")

    desc = (
        f"ref={ref_id} reference_length={len(merged)} covered_start={start} "
        f"covered_end={end} overlap_positions={overlap_count} "
        f"agree_positions={agree_count} conflicts={conflict_count} N={n_count}"
    )
    _write_fasta(ref_len_fa, seq_id, desc, merged)
    _write_fasta(
        span_fa,
        f"{seq_id}_covered_span",
        f"ref={ref_id} ref_start={start} ref_end={end} length={len(covered_span)}",
        covered_span,
    )

    with report_tsv.open("w", encoding="utf-8") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"seq_id\t{seq_id}\n")
        fh.write(f"ref_id\t{ref_id}\n")
        fh.write(f"reference_length\t{len(merged)}\n")
        fh.write(f"covered_start_1based\t{start}\n")
        fh.write(f"covered_end_1based\t{end}\n")
        fh.write(f"covered_span_length\t{len(covered_span)}\n")
        fh.write(f"N_count_reference_length\t{n_count}\n")
        fh.write(f"overlap_positions\t{overlap_count}\n")
        fh.write(f"agree_positions\t{agree_count}\n")
        fh.write(f"conflict_iupac_positions\t{conflict_count}\n")
        fh.write("\nsource_id\tsource_fasta\tref_id\tlength\tcovered_bases\tcovered_start_1based\tcovered_end_1based\n")
        for row in source_rows:
            fh.write("\t".join(str(row[k]) for k in row) + "\n")

    with per_pos_tsv.open("w", encoding="utf-8") as fh:
        cols = ["pos_1based", "consensus", "coverage", "status", "evidence"]
        fh.write("\t".join(cols) + "\n")
        for row in per_pos:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    print(f"[overlap] wrote {ref_len_fa}")
    print(f"[overlap] wrote {span_fa}")
    print(f"[overlap] wrote {report_tsv}")
    print(f"[overlap] wrote {per_pos_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
