"""
Classify Sanger sequencing reads by Phred quality.

Reads are paired by gene (forward + reverse of the same amplicon belong
together). The pair is classified by the LOWER mean Phred score of its two
reads, so that both forward and reverse end up in the same tier folder.

Raw Sanger deliveries are read from I_Raw_Sequence_Data/ under the
PCR_Sequence root. Tier folders are created under II_QualityClassified/:

    1_HighQuality/    pair min mean Q >= 40   (>= 99.99% per-base accuracy)
    2_MediumQuality/  20 <= pair min mean Q < 40
    3_LowQuality/     pair min mean Q < 20    (< 99% per-base accuracy)

A summary CSV (quality_report.csv) is written to the PCR_Sequence root with
one row per individual read plus the pair tier it was assigned to.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_DATA_ROOT = ROOT / "I_Raw_Sequence_Data"
OUTPUT_ROOT = ROOT / "II_QualityClassified"

HIGH_THRESHOLD = 40.0
MEDIUM_THRESHOLD = 20.0

TIER_DIRS = {
    "HighQuality": "1_HighQuality",
    "MediumQuality": "2_MediumQuality",
    "LowQuality": "3_LowQuality",
}


def parse_phd_qualities(phd_path: Path) -> list[int]:
    """Return the list of per-base Phred quality scores from a .phd.1 file."""
    qualities: list[int] = []
    in_dna = False
    with phd_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN_DNA":
                in_dna = True
                continue
            if line == "END_DNA":
                break
            if not in_dna or not line:
                continue
            parts = line.split()
            # Format: <base> <quality> <peak_position>
            if len(parts) >= 2:
                try:
                    qualities.append(int(parts[1]))
                except ValueError:
                    pass
    return qualities


def classify(mean_q: float) -> str:
    if mean_q >= HIGH_THRESHOLD:
        return "HighQuality"
    if mean_q >= MEDIUM_THRESHOLD:
        return "MediumQuality"
    return "LowQuality"


def gene_key(stem: str) -> tuple[str, str | None]:
    """Split a read stem into (gene_id, direction).

    A stem like ``SmelDMP1_730_Primer1_730_forward`` becomes
    ``("SmelDMP1_730_Primer1_730", "forward")``. If neither suffix is
    present, direction is ``None`` and the stem itself is treated as the
    gene id.
    """
    for direction in ("forward", "reverse"):
        suffix = f"_{direction}"
        if stem.endswith(suffix):
            return stem[: -len(suffix)], direction
    return stem, None


def main() -> None:
    report_rows: list[dict[str, str]] = []
    if not RAW_DATA_ROOT.is_dir():
        raise SystemExit(f"ERROR: raw sequence-data directory not found: {RAW_DATA_ROOT}")

    for sample_dir in sorted(p for p in RAW_DATA_ROOT.iterdir() if p.is_dir()):

        # Collect per-read metrics first, grouped by gene id
        groups: dict[str, list[dict]] = {}
        for phd_file in sorted(sample_dir.glob("*.phd.1")):
            stem = phd_file.name[: -len(".phd.1")]
            qualities = parse_phd_qualities(phd_file)

            if not qualities:
                mean_q = 0.0
                n_bases = 0
                pct_q20 = 0.0
                pct_q30 = 0.0
            else:
                n_bases = len(qualities)
                mean_q = sum(qualities) / n_bases
                pct_q20 = 100.0 * sum(1 for q in qualities if q >= 20) / n_bases
                pct_q30 = 100.0 * sum(1 for q in qualities if q >= 30) / n_bases

            gene, direction = gene_key(stem)
            groups.setdefault(gene, []).append(
                {
                    "direction": direction,
                    "stem": stem,
                    "mean_q": mean_q,
                    "n_bases": n_bases,
                    "pct_q20": pct_q20,
                    "pct_q30": pct_q30,
                }
            )

        if not groups:
            continue

        # Classify each gene by the lower mean Q of its reads, then copy
        # both forward and reverse into the same tier folder.
        for gene, reads in sorted(groups.items()):
            pair_min_q = min(r["mean_q"] for r in reads)
            tier = classify(pair_min_q)
            tier_folder = TIER_DIRS[tier]

            dest_dir = OUTPUT_ROOT / tier_folder / sample_dir.name
            dest_dir.mkdir(parents=True, exist_ok=True)

            for r in reads:
                stem = r["stem"]
                for companion in sample_dir.glob(f"{stem}.*"):
                    shutil.copy2(companion, dest_dir / companion.name)

                report_rows.append(
                    {
                        "sample": sample_dir.name,
                        "gene": gene,
                        "direction": r["direction"] or "",
                        "read": stem,
                        "bases": str(r["n_bases"]),
                        "mean_phred": f"{r['mean_q']:.2f}",
                        "pair_min_phred": f"{pair_min_q:.2f}",
                        "pct_q20": f"{r['pct_q20']:.1f}",
                        "pct_q30": f"{r['pct_q30']:.1f}",
                        "tier": tier_folder,
                    }
                )

                print(
                    f"{sample_dir.name}/{stem}: "
                    f"n={r['n_bases']}, mean Q={r['mean_q']:.2f}, "
                    f"pair min Q={pair_min_q:.2f} -> {tier_folder}"
                )

    report_path = ROOT / "quality_report.csv"
    with report_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample",
                "gene",
                "direction",
                "read",
                "bases",
                "mean_phred",
                "pair_min_phred",
                "pct_q20",
                "pct_q30",
                "tier",
            ],
        )
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"\nWrote {len(report_rows)} rows to {report_path}")


if __name__ == "__main__":
    main()
