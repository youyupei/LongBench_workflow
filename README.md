# LongBench Analysis Workflow

> **Note — code reference, not a plug-and-run pipeline.**
> This repository documents the analysis code used in the LongBench study.
> It is not intended to be executed directly on new data without modification:
> input paths, sample identifiers, and several software locations are
> hard-coded to the original project environment. The code is made available
> so that methods can be inspected or adapted for
> related projects.
>
> **Reuse:** Individual sub-workflows are the most self-contained units and
> can be adapted with targeted changes to the config files (see
> [Configuration system](#configuration-system) below for what to modify).
> Examples of reusable components are given in the
> [Sub-workflows](#sub-workflows) section.

Snakemake workflow for the LongBench benchmarking project. Covers preprocessing,
QC, quantification, variant calling, and downstream R analysis across four sequencing
modalities: long-read bulk, long-read single-cell/nuclei, short-read bulk, and
short-read single-cell/nuclei.

---

## Repository layout

```
workflow/
├── Snakefile                   # Main entry point — imports all sub-workflows and rules
├── config/
│   ├── config.yaml             # Master config (paths, software, email)
│   ├── config_lr_bulk.yaml     # Per-sub-workflow sample/parameter configs
│   ├── config_lr_sc_sn.yaml
│   ├── config_sr_bulk.yaml
│   ├── config_sr_sc_sn.yaml
│   ├── config_hybrid_bulk.yaml
│   ├── config_main_wf.yaml
│   └── config_parser.py        # Loads + resolves relative paths in configs
├── sub_workflows/              # Self-contained per-modality pipelines (see below)
├── rules/                      # Cross-modality rules run from the main workflow
├── modules/
│   └── utilities.smk           # Shared rules (BAM sort/index, subsample, git helpers)
├── scripts/                    # Standalone analysis scripts (see below)
├── rmarkdown/                  # R Markdown reports and their preprocessing R scripts
├── conda/config/               # Conda environment definitions keyed by tool name
├── singularity/                # Singularity container definitions
├── profile/                    # Snakemake SLURM profile (cluster submission settings)
└── custom_analysis/            # One-off analyses not part of the main DAG
```

---

## Sub-workflows

Each sub-workflow under `sub_workflows/` is a self-contained Snakemake module with
its own `Snakefile`, `rules/`, `config/`, and `scripts/`. They are imported into the
main `Snakefile` via `module` declarations and all their rules are namespaced
(e.g. `lr_bulk_*`).

| Sub-workflow | Modality | Key tools | Key rules |
|---|---|---|---|
| `lr_bulk` | Long-read bulk (ONT cDNA, PacBio, dRNA) | minimap2, oarfish, IsoQuant, NanoPlot, clair3 | minimap2, qc, quantification, internal_priming, mutation, subsample |
| `lr_sc_sn` | Long-read single-cell / single-nuclei | FLAMES, minimap2, demuxlet, oarfish, clair3, whatshap | flames, qc, demuxlet, pseudo_bulk_map_n_quant, pseudo_bulk_qc, mutation, sc_clustering |
| `sr_bulk` | Short-read bulk (Illumina) | STAR, salmon, fastp, clair3 | mapping, qc, quantification, mutation |
| `sr_sc_sn` | Short-read single-cell / single-nuclei | CellRanger, demuxlet | cellranger, demuxlet, qc |
| `hybrid_bulk` | Hybrid long+short read quantification | MiniQuant | miniquant |

Each sub-workflow's `Snakefile` is self-contained and loads its own config, making
them the most practical units for reuse. The table below gives concrete examples of
what can be adapted with config-only changes vs. what requires code edits.

| Goal | Sub-workflow to start from | What to change |
|---|---|---|
| Run the ONT/PacBio bulk preprocessing pipeline on new samples | `lr_bulk` | `config_lr_bulk.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`, and `barcode_list`; update `minimap2` path in `config.yaml` if not using the same binary |
| Apply the FLAMES single-cell pipeline to a new ONT SC dataset | `lr_sc_sn` | `config_lr_sc_sn.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`; provide a new genome/GTF under the `reference:` block |
| Re-run short-read bulk quantification (STAR + salmon) on new samples | `sr_bulk` | `config_sr_bulk.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`; point `star_index` and `salmon_index` to pre-built indices for the new genome |
| Reuse the variant-calling rules (clair3 + whatshap) alone | `lr_bulk` — `rules/mutation.smk` | Include only `mutation.smk` in a new Snakefile; the rule inputs are standard sorted BAMs so no other `lr_bulk` rules are required |
| Adapt the junction saturation QC plots for a different project | `scripts/junctionSaturation_plot_known.R` | The script reads paths from `snakemake@input`; replace the input block or pass paths as arguments to run standalone |

---

## Main workflow rules (`rules/`)

Rules that aggregate or compare across modalities, run after the sub-workflows complete.

| File | What it does |
|---|---|
| `rmarkdown.smk` | Knits all R Markdown reports via `versioned_knit.R`; defines the preprocessing R script rules (`r_tx2gene_map`, `r_get_bulk_DGE_objects`, etc.) |
| `qc_plot.smk` | Combined QC plots across all modalities |
| `variant_analysis.smk` | Variant calling evaluation (clair3, longcallR, CCLE truth sets) |
| `rarefraction_analysis.smk` | Rarefaction / downsampling analysis |
| `base_count_analysis.smk` | Per-sample base count aggregation |
| `sc_cell_line_anno.smk` | Single-cell cell-line annotation pipeline |
| `rule_recycle_bin.smk` | Inactive rules kept for reference — not included in the main workflow |

---

## Shared utilities (`modules/utilities.smk`)

Importable rules available to any sub-workflow via `include:`:

- `utilities_sort_and_index_bam` — samtools sort + index
- `utilities_subsample_bam` — samtools fractional subsample
- `utilities_gtf_to_bed` — GTF → BED conversion via gxf2bed
- `utilities_git_check_commit` / `utilities__git_clone` — pin and auto-update external git repos

---

## Scripts (`scripts/`)

Standalone scripts not directly called by Snakemake rules, organised by topic.

### QC and read-level summaries
| Script | Called by | Purpose |
|---|---|---|
| `read_length_and_quality_plot.R` | `qc_plot.smk` | Per-sample read length and quality distributions |
| `read_count_plot.R` | `qc_plot.smk` | Read count summary across bulk and SC samples |
| `qc_read_length_quality_summary.R` | sub-workflow qc rules | Tabulate length/quality from NanoPlot/FastQC output |
| `qc_read_count_summary.R` | sub-workflow qc rules | Parse and aggregate read count JSON/log files |
| `junctionSaturation_plot_known.R` | `qc_plot.smk` | Junction saturation curves (RSeQC output) |

### Rarefaction analysis
| Script | Called by | Purpose |
|---|---|---|
| `rarefraction_analysis.R` | `rarefraction_analysis.smk` | Compute rarefaction curves from count matrices |
| `rarefraction_plot_only.R` | standalone | Re-plot rarefaction results without re-running the analysis |

### Utility scripts
| Script | Called by | Purpose |
|---|---|---|
| `aggregate_base_counts.py` | `base_count_analysis.smk` | Aggregate per-sample base-count files from gigabase-downsampling runs |
| `find_longest_reads_in_bam.py` | standalone / manual | Extract top-N longest mapped reads from a BAM; outputs TSV of read name, length, aligned length, mean quality |
| `Rfunctions.R` | sourced by Rmd files | Shared R helper functions and metadata loading |
| `FLAMES_transcript_coverage_plot.R` | standalone | Transcript coverage visualisation from FLAMES output (requires FLAMES ≥ 2.3.4) |

### Subfolders
| Folder | Contents |
|---|---|
| `DTU_analysis/` | Differential transcript usage analysis (`bulkDTU_analysis.Rmd`, `dtu_Rfunction.R`, `major_isoform_analysis.Rmd`) — standalone, not in the main DAG |
| `fusion_analysis/` | Fusion gene detection with JAFFAL across ONT cDNA, dRNA, and PacBio; shell scripts for adapter trimming, subsampling, and JAFFAL runs, plus a figure-generation Rmd |
| `RNAmod_analysis/` | RNA modification analysis via modkit pileup; includes a Nextflow workflow (`RNAmod_analysis.nf`) and comparison R script |
| `variant_analysis/` | Manual / exploratory shell scripts for clair3 variant calling on spike-in samples and CCLE/FP variant detection (the automated versions are in `rules/variant_analysis.smk`) |

---

## Configuration system

All paths are centralised in `config/config.yaml`. The helper `config/config_parser.py`
resolves relative paths to absolute and merges per-sub-workflow configs so each
sub-workflow only sees its own keys.

To adapt the workflow to new data, edit the relevant `config/config_<sub_wf>.yaml`:
- `samples_fastq_dir` — per-sample input paths
- `sample_id` — sample list
- `output_path` — where results are written
- `subsample_read_n` — optional; triggers the subsampling rules if set

Software not managed by conda (SQANTI3, minimap2, kallisto, bustools, demuxlet) is
pointed to by absolute paths under the `software:` block in `config.yaml`.

---

## Execution context (how the workflow was run)

The workflow was executed on a SLURM cluster using the profile in `profile/config.yaml`.
The commands below are provided as a reference for anyone adapting the code, not as
ready-to-run instructions on new data.

```bash
# Dry-run to inspect the DAG without executing anything
snakemake -n --profile profile/

# Full run (SLURM submission via profile)
snakemake --profile profile/

# Run a single sub-workflow in isolation
cd sub_workflows/lr_bulk
snakemake --profile ../../profile/

# Force re-run of a specific rule
snakemake --profile profile/ -R lr_bulk_minimap2_align
```

Conda environments are stored under `conda/envs/` (referenced by hash, not name).
Singularity containers are under `singularity/`. The two are used together: the
container provides the base OS, and conda sets the tool environment inside it.

---

## R Markdown reports (`rmarkdown/`)

Knitting is managed by `rmarkdown/versioned_knit.R` which passes cache and figure
paths so outputs are versioned. To re-knit a single report:

```bash
module load pandoc
Rscript rmarkdown/versioned_knit.R rmarkdown/<Report>.Rmd
```

The dependency order and intermediate RDS files are documented in
[`rmarkdown/DEPENDENCIES.md`](rmarkdown/DEPENDENCIES.md).
