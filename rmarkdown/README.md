# R Markdown Reports

This directory contains all R Markdown reports for the LongBench analysis, along
with the knitting infrastructure and preprocessing R scripts that produce the
intermediate RDS objects they depend on.

For the full dependency graph between reports and intermediate files, see
[DEPENDENCIES.md](DEPENDENCIES.md).

---

## Knitting reports

All reports are knitted via `versioned_knit.R`, which handles output naming,
cache directories, and figure paths automatically.

```bash
module load pandoc
Rscript versioned_knit.R <Report>.Rmd
```

Output is written to `<Report>_outdir/` in the same directory. If a `version:`
field is set in the Rmd YAML front matter, it is included in the output filename:

```
<Report>_outdir/
├── <Report>_<version>_<date>.html
├── <Report>_<version>_cache/        # knitr chunk cache
└── <Report>_<version>_<date>_figures/   # saved PNG + SVG figures
```

If no `version:` is set, the version component is omitted from the filename.

---

## Rmd conventions

All reports in this directory follow these conventions so they work correctly
with `versioned_knit.R` and can be re-knitted independently.

**Required YAML front matter:**
```yaml
---
title: "..."
version: "1.0"   # bump to invalidate the cache and force a clean run
params:
  cache_dir: NULL
  fig.path: NULL
---
```

**Required setup chunk:**
````r
```{r setup, include=FALSE}
knitr::opts_chunk$set(
  cache    = !is.null(params$cache_dir),
  cache.path = params$cache_dir,
  fig.path = params$fig.path,
  fig.width  = 15,
  fig.height = 10,
  dev = if (!is.null(params$fig.path)) c("png", "svg") else "png"
)
```
````

**Chunk discipline:**
- Separate data-generation chunks (heavy computation, `cache = TRUE`) from
  plotting chunks (`cache = FALSE`). This lets plots be tweaked without
  re-running expensive steps.
- Name all chunks that produce figures — the chunk name becomes the figure
  filename on disk.
- Do not call `ggsave()` or write files directly from chunks; let knitr save
  figures via `fig.path`.

**Cache invalidation:**
Bump the `version:` field in the YAML front matter to force a clean run.
The cache directory name includes the version, so old caches are left intact
and the new run starts fresh.

---

## Preprocessing R scripts (`Rscript/`)

Scripts under `Rscript/` are called by Snakemake rules (defined in
`rules/rmarkdown.smk`) to produce intermediate RDS objects that multiple
reports share. They should not be run manually unless re-running in isolation.

| Script | Output RDS | Used by |
|---|---|---|
| `Tx2Gene.map.R` | `RDS/Tx2Gene.map.rds` | Bulk DGE preprocessing, GC analysis |
| `Bulk.DGElist.preprocessing.R` | `RDS/bulk_DGE.obj.rds` | Most bulk Rmd reports |
| `Bulk.DGElist.preprocessing.IP_filtered.R` | `RDS/bulk_DGE.IP_filtered.obj.rds` | `Bulk_identification_IP_filtered.Rmd` |
| `Sc.DGElist.preprocessing.R` | `RDS/sc_DGE.obj.rds` | `SC_identification_DE_analysis.Rmd` |
| `get_intronic_gene.R` | `RDS/intronic_gene_and_exon_count.rds` | Bulk identification reports |

---

## Report index

### QC
| Report | Description |
|---|---|
| `QC_plot.Rmd` | Combined QC across all modalities (read length, quality, alignment, junction saturation, internal priming) |

### Bulk — transcript identification
| Report | Description |
|---|---|
| `Bulk_identification.Rmd` | Main bulk transcript identification benchmarking (oarfish/FLAMES vs reference) |
| `Bulk_identification_IsoQuant.Rmd` | Same analysis using IsoQuant quantification |
| `Bulk_identification_IP_filtered.Rmd` | Identification after internal-priming filtering |
| `Bulk_identification_3way_comparison.Rmd` | Three-tool comparison of identification results |
| `Bulk_identification_tool_comparison.Rmd` | Broader tool comparison |
| `02_Bulk_Identification.Rmd` | Revised/numbered version of the identification report |

### Bulk — quantification
| Report | Description |
|---|---|
| `Bulk_quantification_analysis.Rmd` | Quantification accuracy benchmarking (oarfish) |
| `Bulk_quantification_analysis_IsoQuant.Rmd` | Quantification with IsoQuant |
| `Bulk_quantification_analysis_MiniQuant.Rmd` | Quantification with MiniQuant |
| `Bulk_quantification_cross_tool_comparison.Rmd` | Cross-tool quantification comparison |
| `04_Bulk_quantification.Rmd` | Revised/numbered quantification report |

### Bulk — differential expression
| Report | Description |
|---|---|
| `Bulk_DE_Summary.Rmd` | Main bulk DE benchmarking summary |
| `Bulk_DE_Summary_20M.Rmd` | DE summary at 20M read depth |
| `Bulk_DE_benchmarking.Rmd` | Detailed DE benchmarking |
| `Bulk_DE_benchmarkin_20M.Rmd` | DE benchmarking at 20M depth |
| `Bulk_DE_benchmarking_IsoQuant.Rmd` | DE benchmarking with IsoQuant |
| `Bulk_DE_IsoQuant.Rmd` | DE analysis using IsoQuant counts |
| `Bulk_DE_analysis_continue.Rmd` | Extended DE analysis (reads `bulk_DE.rds`) |
| `Bulk_DE_analysis_continue_IsoQuant.Rmd` | Extended DE with IsoQuant |
| `Bulk_DE.spikeins.Rmd` | DE analysis on spike-in transcripts |
| `Bulk_DE.spikeins_IsoQuant.Rmd` | Spike-in DE with IsoQuant |
| `Bulk_DE_platform_specificity.Rmd` | Platform-specific DE patterns |
| `Bulk_DE_rarefaction.Rmd` | DE performance vs. sequencing depth |
| `05_Bulk_DE_logFC_analysis.Rmd` | logFC concordance analysis |
| `Bulk_BCV.Rmd` | Biological coefficient of variation analysis |

### Bulk — other
| Report | Description |
|---|---|
| `Bulk_GC_analysis_20M.Rmd` | GC content bias analysis at 20M reads |
| `dRNA_run_comparison.Rmd` | Comparison between two dRNA sequencing runs |

### Single-cell / single-nuclei
| Report | Description |
|---|---|
| `sc_sn_umap.Rmd` | UMAP visualisation for SC and SN data; produces `sc_sn_filtered_so.rds` |
| `sc_clustering_annotation.Rmd` | SC cluster annotation |
| `sn_clustering_annotation.Rmd` | SN cluster annotation |
| `SC_identification_DE_analysis.Rmd` | SC pseudo-bulk identification and DE analysis |

### Variant / mutation
| Report | Description |
|---|---|
| `Mutation_analysis.Rmd` | Variant calling benchmarking (clair3, longcallR) |
| `Mutation_analysis_downsample.Rmd` | Variant calling at downsampled read depths |
