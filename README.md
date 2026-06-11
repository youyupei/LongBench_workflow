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
>
> **HPC dependency — `module load`:** Several rules call `module load` to
> activate tools via the WEHI SLURM environment module system. These calls
> are not portable and will fail on systems without those modules. See the
> [HPC module dependencies](#hpc-module-dependencies) section for a full list
> and suggested replacements.

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

## Figure source map

The table below maps each paper figure to the Rmd or script that generates it and
the Snakemake rule chain that produces its inputs.

#### Figure 1

| Figure | Source | Upstream Snakemake rules |
|---|---|---|
| Fig 1B, 1C | `00_Bulk_BCV.Rmd` | quantify with Salmon (`sr_bulk_salmon_quant`); Rmd reads `quant.sf` files directly |
| Fig 1E, 1F | `01_QC_plot.Rmd` | aggregate read and base counts across all modalities (`qc_read_count_summary`) — pulls from per-FASTQ read counters in lr_bulk (`count_reads_in_fastq`, `count_bases_in_fastq`), FLAMES `summary.txt` in lr_sc_sn, Salmon output in sr_bulk, and CellRanger metrics in sr_sc; also collects NanoPlot length/quality TSVs from all sub-workflows (`qc_read_length_quality_summary`) |

#### Figure 2

| Figure | Source | Upstream Snakemake rules |
|---|---|---|
| Fig 2A | `01_QC_plot.Rmd` | subsample BAM and run AlignQC to produce an alignment-feature HTML report for long reads (`lr_bulk_alignQC_analysis_subsample`) and short reads (`sr_bulk_alignQC_analysis`) |
| Fig 2B | `01_QC_plot.Rmd` | detect internal priming sites on transcript-aligned BAMs produced by minimap2 (`lr_bulk_minimap2_transcript`), identify internal-priming reads with PrimeSpotter (`rules/internal_priming.smk`) |
| Fig 2C | `scripts/FLAMES_transcript_coverage_plot.R` (standalone) | run FLAMES pipeline (barcode demux + transcript assembly + alignment) (`lr_sc_sn_flames`); standalone script reads the FLAMES output BAM directly (requires FLAMES ≥ 2.3.4) |
| Fig 2D | `03_Bulk_quantification.Rmd` | quantify transcript abundance from transcript-aligned BAMs with oarfish for long reads (`lr_bulk_oarfish_cov`) and with Salmon for short reads (`sr_bulk_salmon_quant`); results loaded via shared bulk DGE RDS (`r_get_bulk_DGE_objects`) |
| Fig 2E, 2G | `02_Bulk_Identification.Rmd` | same oarfish (`lr_bulk_oarfish_cov`) and Salmon (`sr_bulk_salmon_quant`) quantification as Fig 2D; results loaded via shared bulk DGE RDS (`r_get_bulk_DGE_objects`) |
| Fig 2F, 2H | `scripts/rarefraction_plot_only.R` via `rarefraction_plot` rule | subsample long-read BAMs at multiple depths and re-run oarfish (`oarfish_cov_rare_fraction_analysis`); subsample short-read FASTQs and re-run Salmon (`salmon_quant_downsample`); compute gene/transcript detection curves (`rarefraction_analysis` → `rarefraction_table`); plot with `rarefraction_plot_only.R` — all in `rules/rarefraction_analysis.smk` |

#### Figure 3

| Figure | Source | Upstream Snakemake rules |
|---|---|---|
| Fig 3A–F | `02_Bulk_Identification.Rmd` | shared bulk DGE RDS — see note above |
| Fig 3G–J | `02_Bulk_Identification.Rmd` (extract commonly detected genes and transcripts) → `04_Bulk_DE_Summary.Rmd` | shared bulk DGE RDS fed into both reports sequentially |

#### Figure 4

| Figure | Source | Upstream Snakemake rules |
|---|---|---|
| Fig 4A | `01_QC_plot.Rmd` | detect cell barcodes from long-read FASTQs with BLAZE (`lr_sc_sn_blaze`); genotype cells at known SNPs with cellsnp-lite and assign each barcode to a donor with Vireo (`Cellsnp_lite_rule` → `vireo_rule` for LR; `Cellsnp_lite_rule_short_read` → `vireo_rule_short_read` for SR) |
| Fig 4B | `05_sc_clustering_annotation.Rmd` | BLAZE barcode detection (`lr_sc_sn_blaze`) → FLAMES alignment and transcript quantification (`lr_sc_sn_flames`) for LR; CellRanger demux + count for SR (`sr_sc_sn_cellranger`); Vireo cell-line assignment for all three platforms (`sc_cell_line_anno.smk`); Rmd produces alluvial barcode-overlap diagram |
| Fig 4C | `07_sc_sn_umap.Rmd` | BLAZE → FLAMES pipeline for long-read gene counts (`lr_sc_sn_blaze` → `lr_sc_sn_flames`); CellRanger feature–barcode matrix for short reads (`sr_sc_sn_cellranger`) |
| Fig 4D | `08_SC_identification_DE_analysis.Rmd` | split BAMs by cell-line barcode list, re-align and quantify as pseudo-bulk samples (`pseudo_bulk_map_n_quant.smk` — minimap2 + oarfish per cell line) |
| Fig 4E–J | `08_SC_identification_DE_analysis.Rmd` | pseudo-bulk inputs (same as 4D) + shared bulk DGE RDS as bulk reference; Fig 4E pseudobulk–bulk correlation is computed at gene and transcript level (log CPM) across all platform combinations |

#### Figure 5

| Figure | Source | Upstream Snakemake rules |
|---|---|---|
| Fig 5A, 5B, 5F | `09_Mutation_analysis.Rmd` | call SNVs/indels from genome-aligned BAMs using clair3-rna (`lr_bulk_clair3_rna`), phase variants against the BAM with WhatsHap (`lr_bulk_whatshap`) — both in `sub_workflows/lr_bulk/rules/mutation.smk`; sr_bulk equivalent clair3 rules; `rules/variant_analysis.smk` aggregates all callsets |
| Fig 5C or 5D | `10_Mutation_analysis_downsample.Rmd` | repeat clair3 variant calling on progressively downsampled BAMs (`custom_analysis/lr_bulk_ds/variant_analysis_ds.smk`) |
| Fig 5D, 5E | `scripts/fusion_analysis/` (standalone) | run JAFFAL v2.5 fusion detection on ONT cDNA, dRNA, and PacBio FASTQs via numbered SLURM shell scripts; no Snakemake integration |
| Fig 5G–J | `09_Mutation_analysis.Rmd` | call variants and phase a BAM with LongCallR (`longcallR`); run allele-specific junction analysis on the phased BAM (`longcallR_analsyis_asj`); run allele-specific expression (`longcallR_analsyis_ase`) — all in `sub_workflows/lr_bulk/rules/mutation.smk` |

---

## Sub-workflows

Each sub-workflow under `sub_workflows/` is a Snakemake module with its own
`Snakefile`, `rules/`, `config/`, and `scripts/`. They are imported into the main
`Snakefile` via `module` declarations and all their rules are namespaced
(e.g. `lr_bulk_*`). Each sub-workflow is self-contained and loads its own config,
making them the most practical units for reuse — though adapting one to a new dataset
requires updating paths and sample identifiers in the relevant config file.

| Sub-workflow | Modality | Key tools | Key rules |
|---|---|---|---|
| `lr_bulk` | Long-read bulk (ONT cDNA, PacBio, dRNA) | minimap2, oarfish, IsoQuant, NanoPlot, clair3, LongCallR | minimap2, qc, quantification, internal_priming, mutation, subsample |
| `lr_sc_sn` | Long-read single-cell / single-nuclei | FLAMES, minimap2, oarfish, cellsnp-lite, Vireo, clair3, whatshap | flames, qc, pseudo_bulk_map_n_quant, pseudo_bulk_qc, mutation, sc_clustering |
| `sr_bulk` | Short-read bulk (Illumina) | Rsubread (Subjunc), salmon, fastp, clair3 | mapping, qc, quantification, mutation |
| `sr_sc_sn` | Short-read single-cell / single-nuclei | CellRanger, cellsnp-lite, Vireo | cellranger, qc |
| `hybrid_bulk` | Hybrid long+short read quantification | MiniQuant | miniquant |

---

## Reusable rules

The rules below are general enough to lift into another project with minimal changes.
The "Minimum config keys" column lists exactly what must be present in the config for
the rule to work; everything else in the original configs is LongBench-specific.

### Long-read alignment — `sub_workflows/lr_bulk/rules/minimap2.smk`

Rules: `lr_bulk_minimap2_transcript`, `lr_bulk_minimap2_Genome`

Aligns per-sample FASTQ files to a transcript or genome reference using minimap2,
producing sorted BAMs. The preset flags (e.g. `-ax splice`) are passed through
`minimap2_trans_options` / `minimap2_genome_options` so they are easy to swap for
a different chemistry.

**Minimum config keys:**
```yaml
samples_fastq_dir:
  my_sample: "/path/to/fastqs"       # directory containing per-cell-line FASTQ files
sample_id:
  - my_sample
barcode_list:                        # maps barcode number to cell-line label
  - 1: SampleA
cell_lines:
  - SampleA
reference:
  transcript: "/path/to/transcriptome.fa"
  genome:     "/path/to/genome.fa"
software:
  minimap2: "/path/to/minimap2"
minimap2_trans_options:
  my_sample: "-ax map-ont --secondary=no"
minimap2_genome_options:
  my_sample: "-ax splice --secondary=no"
output_path: "/path/to/results"
```

---

### Long-read QC suite — `sub_workflows/lr_bulk/rules/qc.smk`

Rules: `NanoPlot`, `RSeQC_junction_saturation`, `RSeQC_junction_annotation`,
`RSeQC_gene_body_coverage`, `alignQC_analysis_subsample`, `count_reads_in_fastq`,
`count_bases_in_fastq`

A self-contained QC battery for long-read RNA-seq. `NanoPlot` runs on subsampled
FASTQs; the RSeQC rules run on genome-aligned BAMs; `alignQC` runs inside a
Singularity container (`docker://vacation/alignqc`) and produces a full HTML report.
Each rule can be included independently — they do not form a strict chain.

**Minimum config keys:**
```yaml
reference:
  bed_human: "/path/to/annotation.bed"   # for RSeQC; convert GTF with utilities_gtf_to_bed
  genome:    "/path/to/genome.fa"        # for AlignQC
  gtf_gz:   "/path/to/annotation.gtf.gz" # for AlignQC
conda:
  NanoPlot: "/path/to/conda/NanoPlot.yaml"
  RSeQC:    "/path/to/conda/RSeQC.yaml"
output_path: "/path/to/results"
scratch_dir: "/path/to/scratch"          # intermediate subsampled BAMs written here
random_seed: 2024
```

---

### Transcript quantification — `sub_workflows/lr_bulk/rules/quantification.smk`

Rules: `oarfish_cov`, `salmon`

`oarfish_cov` quantify transcript abundance from
transcript-aligned BAMs using EM (with optional coverage correction).
`salmon` runs in alignment-based mode on the same BAMs.
All three take a transcript-sorted BAM as input, so they depend on the
minimap2 transcript alignment rule above but nothing else.

**Minimum config keys:**
```yaml
reference:
  transcript: "/path/to/transcriptome.fa"
  gtf:        "/path/to/annotation.gtf"   # featureCounts only
conda:
  oarfish: "/path/to/conda/oarfish.yaml"
  main:    "/path/to/conda/main.yaml"     # salmon / featureCounts
output_path: "/path/to/results"
```
---

### Variant calling — `sub_workflows/lr_bulk/rules/mutation.smk`

Rules: `clair3_rna`, `whatshap`

`clair3_rna` calls SNVs/indels from a genome-aligned BAM using the clair3-rna
Singularity container. `whatshap` phases the resulting VCF back against the BAM.
The platform model is selected per sample via the `_clair3_rna_platform` dict at
the top of the file — edit that dict to add new platform presets.

Both rules take a sorted, indexed genome BAM as their only upstream dependency,
so they can be dropped into any workflow that produces one.

**Minimum config keys:**
```yaml
reference:
  genome: "/path/to/genome.fa"
conda:
  whatshap: "/path/to/conda/whatshap.yaml"
output_path: "/path/to/results"
scratch_dir: "/path/to/scratch"    # clair3 VCFs written here first to save project storage
```

Also update the `_clair3_rna_platform` dict in `mutation.smk` to map your sample
names to a clair3-rna platform string (e.g. `ont_r10_dorado_cdna`, `hifi_mas_minimap2`).
---

## HPC module dependencies

Several rules use `module load` to activate software through the WEHI HPC
environment module system. This is **not portable** — the commands will fail on
any system without those specific modules installed. When adapting these rules,
replace each `module load` call with the appropriate `conda:` directive or a
Singularity container.

**How to add a `conda:` directive to an existing rule:**
```python
rule my_rule:
    input: ...
    output: ...
    conda: "path/to/env.yaml"   # relative to the Snakefile, or absolute
    shell: "my-tool ..."
```

Conda env YAML files for tools already used in this workflow are in `conda/config/`.

---

## R Markdown reports (`rmarkdown/`)
 see [`rmarkdown/README.md`](rmarkdown/README.md) for conventions and the full report index.
