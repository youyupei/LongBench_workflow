# R Markdown Scripts

This directory contains all R Markdown reports for the LongBench analysis

## Preprocessing R scripts (`Rscript/`)

Scripts under `Rscript/` are called by Snakemake rules (defined in
`rules/rmarkdown.smk`) to produce intermediate RDS objects that multiple
reports share. They should not be run manually unless re-running in isolation.

| Script | Output RDS | Used by |
|---|---|---|
| `Tx2Gene.map.R` | `RDS/Tx2Gene.map.rds` | Bulk DGE preprocessing, GC analysis |
| `Bulk.DGElist.preprocessing.R` | `RDS/bulk_DGE.obj.rds` | Most bulk Rmd reports |
| `Bulk.DGElist.preprocessing.IP_filtered.R` | `RDS/bulk_DGE.IP_filtered.obj.rds` | *(report removed; RDS retained for reference)* |
| `Sc.DGElist.preprocessing.R` | `RDS/sc_DGE.obj.rds` | `SC_identification_DE_analysis.Rmd` |
| `get_intronic_gene.R` | `RDS/intronic_gene_and_exon_count.rds` | Bulk identification reports |

---

## Rmd index

### Main paper figures (numbered Rmds — run manually via `versioned_knit.R`)

| Report | Paper figures | Description |
|---|---|---|
| `00_Bulk_BCV.Rmd` | Fig 1B, 1C | BCV / dispersion plots using short-read Salmon quantification |
| `01_QC_plot.Rmd` | Fig 1E, 1F; 2A, 2B; 4A | Combined QC across all modalities: read length/quality (NanoPlot), alignment (AlignQC), internal priming, SC/SN Vireo doublet/ambient stats |
| `02_Bulk_Identification.Rmd` | Fig 2E, 2G; 3A–F; 3G–J (input) | Bulk transcript identification benchmarking (long-read oarfish + short-read Salmon) |
| `03_Bulk_quantification.Rmd` | Fig 2D | Over-dispersion and raw count distribution benchmarking |
| `04_Bulk_DE_Summary.Rmd` | Fig 3G–J | Bulk differential expression summary |

### Single-cell / single-nuclei

| Report | Paper figures | Description |
|---|---|---|
| `05_sc_clustering_annotation.Rmd` | Fig 4B | Cell and cell-line annotation; alluvial barcode-overlap plot (BLAZE + CellRanger + Vireo across ONT, PacBio, Illumina) |
| `06_sn_clustering_annotation.Rmd` | — | SN equivalent of `05_sc_clustering_annotation.Rmd` |
| `07_sc_sn_umap.Rmd` | Fig 4C | UMAP visualisation for SC and SN; reads FLAMES gene counts (LR) and CellRanger feature matrix (SR) |
| `08_SC_identification_DE_analysis.Rmd` | Fig 4D, 4E–J | SC pseudo-bulk identification and DE; pseudobulk–bulk correlation at gene and transcript level |

### Variant / mutation

| Report | Paper figures | Description |
|---|---|---|
| `09_Mutation_analysis.Rmd` | Fig 5A, 5B, 5F, 5G–J | Variant calling benchmarking: clair3-RNA, clair3-Illumina, whatshap phasing, LongCallR |
| `10_Mutation_analysis_downsample.Rmd` | Fig 5C or 5D | Variant calling at downsampled read depths (uses `custom_analysis/lr_bulk_ds/variant_analysis_ds.smk`) |

### Supplementary / exploratory

| Report | Description |
|---|---|
| `Sup_Bulk_quantification_analysis_MiniQuant.Rmd` | Quantification benchmarking with MiniQuant |
| `Sup_Bulk_identification_3way_comparison.Rmd` | Three-tool transcript identification comparison (Oarfish vs )|
| `Sup_Bulk_quantification_analysis_IsoQuant.Rmd` | Quantification benchmarking with IsoQuant |
| `Sup_Bulk_DE_Summary_20M.Rmd` | DE summary at 20M read depth |
| `Sup_Bulk_DE_IsoQuant.Rmd` | DE analysis using IsoQuant counts |
| `Sup_Bulk_DE.spikeins.Rmd` | DE analysis on spike-in transcripts |
| `Sup_Bulk_DE.spikeins_IsoQuant.Rmd` | Spike-in DE with IsoQuant |
| `Sup_Bulk_DE_logFC_analysis.Rmd` | logFC concordance analysis |
| `Sup_Bulk_GC_analysis_20M.Rmd` | GC content bias analysis at 20M reads |
| `Sup_dRNA_run_comparison.Rmd` | Comparison between two dRNA sequencing runs |
