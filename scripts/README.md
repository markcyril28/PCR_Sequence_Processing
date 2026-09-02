# Script guide

Scripts are grouped by workflow. Number prefixes show the normal order within
each folder; numbering restarts for each workflow.

Run complete workflows through the shell drivers in the workspace root. The
drivers already resolve the numbered script paths, so they can be launched
from any working directory.

## Layout

```text
scripts/
├── common/
│   ├── 01_read_toml.py
│   ├── 02_alignment_viz.py
│   └── 03_trim_reads.py
├── merge_consensus/
│   └── 01_pairwise_consensus.py
├── low_quality_saver/
│   ├── 01_step_classify.py
│   ├── 02_ab1_to_fasta.py
│   ├── 03_step_trace_review.py
│   ├── 04_consensus_reviewed.py
│   ├── 05_step_consensus_assembly.py
│   ├── 06_step_tracy_decompose.py
│   ├── 07_step_merge_consensus.py
│   ├── 08_blast_consensus.sh
│   ├── 09_step_report.py
│   └── 10_curator_apply_overrides.py
├── reporting/
│   └── 01_build_full_report.py
├── tri_consensus/
│   └── 01_tri_consensus.py
├── inside_primer/
│   ├── 01_overlap_consensus.py
│   └── 02_overlap_consensus_report.py
└── standalone/
    ├── 01_classify_by_quality.py
    ├── 02_blast_classified.py
    ├── 03_merge_consensus.py
    └── 04_tracy_json_to_fa.py
```

The `__init__.py` files are package markers, not executable workflow steps,
so they do not receive number prefixes.

## Workflow entry points

| Root driver | Script folder | Purpose |
| --- | --- | --- |
| `a_merge_consensus_pipeline.sh` | `merge_consensus/` | Merge paired forward/reverse ABI reads without a reference |
| `b_low_quality_saver_pipeline.sh` | `low_quality_saver/` | Classify, assign, review, assemble, merge, BLAST, and summarize lower-quality reads |
| `c_build_report.sh` | `reporting/` | Build per-batch Markdown/PDF reports |
| `d_tri_consensus_pipeline.sh` | `tri_consensus/` | Combine two Sanger results with the reference |
| `f_inside_primer_pipeline.sh` | `inside_primer/` | Stitch inside/outside-primer products and report the result |

Every driver that reads TOML uses `common/01_read_toml.py`. Shared alignment
formatting and read trimming live in `common/02_alignment_viz.py` and
`common/03_trim_reads.py`.

## Low-quality-saver ordering

The method scripts numbered `03`, `05`, and `06` are alternative calling
methods selected by configuration, rather than three mandatory sequential
steps. `04_consensus_reviewed.py` is both the helper used by method `03` and a
standalone CLI. `10_curator_apply_overrides.py` is an optional manual curation
step after report review.

## Standalone tools

Files under `standalone/` are not called by the current root workflow drivers.
They remain available for manual runs and compatibility with earlier analyses:

- `01_classify_by_quality.py`: earlier workspace-wide quality classifier.
- `02_blast_classified.py`: BLAST reads already arranged by quality tier.
- `03_merge_consensus.py`: merge Tracy strand FASTA/QUAL files against a
  reference.
- `04_tracy_json_to_fa.py`: convert Tracy JSON basecalls to FASTA and QUAL.

## Python imports

Python identifiers cannot begin with a digit. Consumers therefore load
numbered helper modules with `importlib` and their canonical package names,
such as `scripts.common.03_trim_reads`. This also keeps worker imports stable
when multiprocessing starts a new Python interpreter.

| Consumer | Numbered helper modules |
| --- | --- |
| `merge_consensus/01_pairwise_consensus.py` | `common/02_alignment_viz.py`, `common/03_trim_reads.py` |
| `low_quality_saver/02_ab1_to_fasta.py` | `common/03_trim_reads.py` |
| `low_quality_saver/03_step_trace_review.py` | `common/02_alignment_viz.py`, `common/03_trim_reads.py`, `low_quality_saver/04_consensus_reviewed.py` |
| `low_quality_saver/04_consensus_reviewed.py` | `common/02_alignment_viz.py` |
| `low_quality_saver/05_step_consensus_assembly.py` | `common/02_alignment_viz.py`, `common/03_trim_reads.py` |
| `reporting/01_build_full_report.py` | `common/03_trim_reads.py` |
| `tri_consensus/01_tri_consensus.py` | `common/02_alignment_viz.py` |

## Maintenance conventions

- Keep number prefixes two digits and order new files within their workflow.
- Use canonical `scripts.<folder>.<numbered_module>` names for Python imports.
- Update the root driver and this guide whenever a script is renamed or moved.
- Keep shared modules free of work at import time.
- Do not commit `__pycache__` directories or `*.pyc` files.
