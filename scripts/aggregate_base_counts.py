#!/usr/bin/env python3
"""
aggregate_base_counts.py — Phase A of the gigabase-downsampling analysis.

Reads all per-sample base-count files produced by count_bases_in_fastq /
count_bases_in_fastq_sr, builds a summary table (rows = cell lines,
columns = platforms/samples), computes the global minimum total bases across
all platform × cell-line combinations, and writes:

  <summary_tsv>  — human-readable TSV with one row per cell line × sample
  <target_txt>   — single integer: the suggested target_bases for downsampling

Usage
-----
python3 aggregate_base_counts.py <summary_tsv> <target_txt> <count_file> [<count_file> ...]

The count files are named  {results_dir}/qc/base_counts/{sample}_{cell_line}.total_bases
and each contains a single integer (total bases in that FASTQ).
"""

import sys
import os
import re
import pandas as pd


def parse_count_file(path: str) -> tuple[str, str, int]:
    """
    Return (sample, cell_line, total_bases) inferred from the file path.

    Expected tail of path:  .../qc/base_counts/{sample}_{cell_line}.total_bases
    The cell_line is one of the 8 known LongBench lines; everything before the
    last underscore-delimited cell-line token is the sample name.
    """
    CELL_LINES = {"H146", "H69", "H526", "H211", "SHP77", "H1975", "H2228", "HCC827"}

    basename = os.path.basename(path)
    stem = basename.replace(".total_bases", "")

    # Split from the right: last token must be a known cell line
    parts = stem.rsplit("_", 1)
    if len(parts) != 2 or parts[1] not in CELL_LINES:
        raise ValueError(
            f"Cannot parse sample/cell_line from filename '{basename}'. "
            f"Expected pattern: {{sample}}_{{cell_line}}.total_bases"
        )
    sample, cell_line = parts

    with open(path) as fh:
        content = fh.read().strip()
    total_bases = int(content)

    return sample, cell_line, total_bases


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    summary_tsv = sys.argv[1]
    target_txt  = sys.argv[2]
    count_files = sys.argv[3:]

    records = []
    for path in count_files:
        try:
            sample, cell_line, bases = parse_count_file(path)
            records.append({"sample": sample, "cell_line": cell_line, "total_bases": bases})
        except Exception as exc:
            print(f"WARNING: skipping {path} — {exc}", file=sys.stderr)

    if not records:
        print("ERROR: no valid count files parsed", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(records)

    # Pivot: rows = cell_line, columns = sample
    pivot = df.pivot_table(index="cell_line", columns="sample", values="total_bases", aggfunc="first")

    # Compute per-row (cell line) min across platforms, then global min
    pivot["row_min"] = pivot.min(axis=1)
    global_min = int(pivot["row_min"].min())

    # Write summary TSV (drop the helper column)
    os.makedirs(os.path.dirname(summary_tsv) or ".", exist_ok=True)
    pivot.to_csv(summary_tsv, sep="\t")
    print(f"Summary written to: {summary_tsv}")

    # Write suggested target
    os.makedirs(os.path.dirname(target_txt) or ".", exist_ok=True)
    with open(target_txt, "w") as fh:
        fh.write(str(global_min) + "\n")
    print(f"Suggested target bases: {global_min:,}  →  {target_txt}")


if __name__ == "__main__":
    main()
