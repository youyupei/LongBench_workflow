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

## Sub-workflows

Each sub-workflow under `sub_workflows/` is a self-contained Snakemake module with
its own `Snakefile`, `rules/`, `config/`, and `scripts/`. They are imported into the
main `Snakefile` via `module` declarations and all their rules are namespaced
(e.g. `lr_bulk_*`).

| Sub-workflow | Modality | Key tools | Key rules |
|---|---|---|---|
| `lr_bulk` | Long-read bulk (ONT cDNA, PacBio, dRNA) | minimap2, oarfish, IsoQuant, NanoPlot, clair3, LongCallR | minimap2, qc, quantification, internal_priming, mutation, subsample |
| `lr_sc_sn` | Long-read single-cell / single-nuclei | FLAMES, minimap2, oarfish, cellsnp-lite, Vireo, clair3, whatshap | flames, qc, pseudo_bulk_map_n_quant, pseudo_bulk_qc, mutation, sc_clustering |
| `sr_bulk` | Short-read bulk (Illumina) | Rsubread (Subjunc), salmon, fastp, clair3 | mapping, qc, quantification, mutation |
| `sr_sc_sn` | Short-read single-cell / single-nuclei | CellRanger, cellsnp-lite, Vireo | cellranger, qc |
| `hybrid_bulk` | Hybrid long+short read quantification | MiniQuant | miniquant |

Each sub-workflow's `Snakefile` is self-contained and loads its own config, making
them the most practical units for reuse. The table below gives concrete examples of
what can be adapted with config-only changes vs. what requires code edits.

| Goal | Sub-workflow to start from | What to change |
|---|---|---|
| Run the ONT/PacBio bulk preprocessing pipeline on new samples | `lr_bulk` | `config_lr_bulk.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`, and `barcode_list`; update `minimap2` path in `config.yaml` if not using the same binary |
| Apply the FLAMES single-cell pipeline to SC/SN dataset | `lr_sc_sn` | `config_lr_sc_sn.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`; provide a new genome/GTF under the `reference:` block |
| Run short-read bulk quantification (Salmon) on new samples | `sr_bulk` | `config_sr_bulk.yaml`: update `samples_fastq_dir`, `sample_id`, `output_path`; point `star_index` and `salmon_index` to pre-built indices for the new genome |
| Reuse the variant-calling rules (clair3 + whatshap) alone | `lr_bulk` — `rules/mutation.smk` | Include only `mutation.smk` in a new Snakefile; the rule inputs are standard sorted BAMs so no other `lr_bulk` rules are required |


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
| `rarefraction_analysis.R` | `rarefraction_analysis.smk` | Compute rarefaction curves from count matrices (note: filename contains a typo — "rarefraction" should be "rarefaction") |
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
| `fusion_analysis/` | Complete standalone fusion gene detection sub-study, not part of the main Snakemake DAG. Contains: (1) numbered SLURM shell scripts for adapter trimming, quality filtering, subsampling, and JAFFAL v2.5 runs for ONT cDNA, dRNA, and PacBio; (2) a figure-generation Rmd (`*_JW.Rmd`) that reads pre-computed JAFFAL results and CCLE translocation ground truth to produce benchmarking figures. To re-run, execute the shell scripts in order within each platform subdirectory, then knit the Rmd. |
| `RNAmod_analysis/` | Standalone RNA modification analysis via modkit pileup. Includes a Nextflow workflow (`RNAmod_analysis.nf`) for running modkit across samples, and an R script (`analysis_modkit_pileup_two_runs.R`) for comparing modification calls between two dRNA sequencing runs. Not integrated into the main Snakemake DAG. |


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

## Configuration system

All paths are centralised in `config/config.yaml`. The helper `config/config_parser.py`
resolves relative paths to absolute and merges per-sub-workflow configs so each
sub-workflow only sees its own keys.


---

## R Markdown reports (`rmarkdown/`)
Script generated most of the paper figures.
The dependency order are documented in
[`rmarkdown/DEPENDENCIES.md`](rmarkdown/DEPENDENCIES.md).
