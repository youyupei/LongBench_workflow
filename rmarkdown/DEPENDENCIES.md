# Rmarkdown Dependency Map

All active Rmd files in `rmarkdown/`. Numbered Rmds (`00_`–`10_`) are knitted
manually via `versioned_knit.R`; `Sup_` Rmds are supplementary analyses.

---

## Dependency DAG

```
Reference files (GTF, SIRV, Sequins, metadata)
│
├──> Rscript/Tx2Gene.map.R  (r_tx2gene_map)  ──────────────────> RDS/Tx2Gene.map.rds
│                                                                          │
│   Workflow outputs (oarfish, salmon, pseudobulk)                        │
│   ├──> Rscript/Bulk.DGElist.preprocessing.R  (r_get_bulk_DGE_objects) ──┤──> RDS/bulk_DGE.obj.rds
│   └──> Rscript/Sc.DGElist.preprocessing.R    (r_get_sc_DGE_objects)  ──┘──> RDS/sc_DGE.obj.rds
│
├──> Rscript/get_intronic_gene.R  (rmd_intronic_gene_and_exon_count)
│         └──> RDS/intronic_gene_and_exon_count.rds
│
│   ┌──────────── RDS/bulk_DGE.obj.rds ──────────────────────────────────────────────┐
│   │                                                                                  │
│   ├──> 02_Bulk_Identification.Rmd ──────────────────> RDS/bulk_identification.rds   │
│   │                                                              │                   │
│   ├──> 04_Bulk_DE_Summary.Rmd ─────────────────────> RDS/bulk_DE.rds               │
│   │                                                              │                   │
│   │                           RDS/bulk_identification.rds ───────┤                   │
│   │                           RDS/bulk_DE.rds ────────────────────────────────────► 08_SC_identification_DE_analysis.Rmd
│   │                                                                                  │  (also needs sc_DGE.obj.rds,
│   ├──> 00_Bulk_BCV.Rmd        (reads sr_bulk salmon directly)                       │   pseudo_bulk CSV)
│   ├──> 03_Bulk_quantification.Rmd  (reads bulk_DGE.obj.rds)                         │
│   └──> Sup_* bulk Rmds             (read bulk_DGE.obj.rds / Tx2Gene / IsoQuant)     │
│                                                                                       │
│   RDS/sc_DGE.obj.rds ────────────────────────────────────────────────────────────► 08_SC_identification_DE_analysis.Rmd
│
│   lr_sc_sn + sr_sc_sn outputs (FLAMES, CellRanger, Vireo) ────┬──> 05_sc_clustering_annotation.Rmd
│                                                                 ├──> 06_sn_clustering_annotation.Rmd
│                                                                 └──> 07_sc_sn_umap.Rmd ──> RDS/sc_sn_filtered_so.rds
│
│   Standalone (no intermediate RDS produced):
│   ├──> 01_QC_plot.Rmd                (AlignQC, NanoPlot, Vireo, IP outputs)
│   ├──> 09_Mutation_analysis.Rmd      (clair3, LongCallR, WhatsHap output dirs)
│   ├──> 10_Mutation_analysis_downsample.Rmd  (custom_analysis/lr_bulk_ds outputs)
│   └──> Sup_dRNA_run_comparison.Rmd   (custom_analysis/dRNA_comparison outputs)
```

---

## Preprocessing R Scripts (`Rscript/`)

### `Rscript/Tx2Gene.map.R`
| | |
|---|---|
| **Snakemake rule** | `r_tx2gene_map` |
| **Reads** | `gencode.v44.annotation.gtf`, SIRV GTF, Sequins GTF |
| **Writes** | `RDS/Tx2Gene.map.rds` |

### `Rscript/Bulk.DGElist.preprocessing.R`
| | |
|---|---|
| **Snakemake rule** | `r_get_bulk_DGE_objects` |
| **Reads** | `RDS/Tx2Gene.map.rds`, oarfish outputs (`lr_bulk/result/oarfish_cov_output/{ont,pb,dRNA}_bulk`), salmon outputs (`sr_bulk/result/salmon/salmon_quant`) |
| **Writes** | `RDS/bulk_DGE.obj.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `Rscript/Sc.DGElist.preprocessing.R`
| | |
|---|---|
| **Snakemake rule** | `r_get_sc_DGE_objects` |
| **Reads** | `RDS/Tx2Gene.map.rds`, pseudobulk oarfish outputs (`lr_sc_sn/result/PseudoBulkOarfishCov/{ont,pb}_{sc,sn}`) |
| **Writes** | `RDS/sc_DGE.obj.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `Rscript/get_intronic_gene.R`
| | |
|---|---|
| **Snakemake rule** | `rmd_intronic_gene_and_exon_count` |
| **Reads** | `gencode.v44.annotation.gtf` |
| **Writes** | `RDS/intronic_gene_and_exon_count.rds` |

### `versioned_knit.R`
| | |
|---|---|
| **Role** | Knitting helper — called manually to render Rmd with versioned output filenames |
| **Usage** | `module load pandoc && Rscript versioned_knit.R <Report>.Rmd` |

---

## Rmd Files

### `00_Bulk_BCV.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 1B, 1C |
| **Reads** | sr_bulk salmon `quant.sf` files (hard-coded path) |
| **Writes** | HTML report only |

### `01_QC_plot.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 1E, 1F; 2A, 2B; 4A |
| **Reads** | AlignQC HTML dirs (`lr_bulk/result/qc/AlignQC/`, `sr_bulk/result/qc/AlignQC/`), internal priming summaries (`lr_bulk/result/int_prim_analysis/`), NanoPlot TSVs (all sub-workflows), Vireo outputs (`lr_sc_sn/result/vireo/`), barcode lists, `figures/qc/read_count_summary.tsv` (from `qc_read_count_summary` rule) |
| **Writes** | HTML report; optional CSV raw-data exports |

### `02_Bulk_Identification.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 2E, 2G; 3A–F; 3G–J (input to `04_Bulk_DE_Summary.Rmd`) |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/intronic_gene_and_exon_count.rds` |
| **Writes** | `RDS/bulk_identification.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `03_Bulk_quantification.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 2D |
| **Reads** | `RDS/bulk_DGE.obj.rds` |
| **Writes** | HTML report only |

### `04_Bulk_DE_Summary.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 3G–J |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `reference_files/Expected_DE.csv` |
| **Writes** | `RDS/bulk_DE.rds` |

### `05_sc_clustering_annotation.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 4B |
| **Reads** | FLAMES gene counts (`lr_sc_sn/result/flames_out/ont_sc/gene_count.csv`, `pb_sc/gene_count.csv`), CellRanger output (`illumina_sc_sn/result/cellranger/ill_sc/`), Vireo TSVs (`lr_sc_sn/result/vireo/{ont,pb,ill}_sc_vireo/`) |
| **Writes** | `sc_cell_line_annotation.csv`, `sc_cell_line_annotation_post_qc.csv`; optionally an annotated Seurat object RDS (`out_annotated_so_rds`, default NULL) |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |

### `06_sn_clustering_annotation.Rmd`
| | |
|---|---|
| **Paper figures** | — (SN equivalent of `05_sc_clustering_annotation.Rmd`) |
| **Reads** | FLAMES gene counts (`ont_sn`, `pb_sn`), CellRanger SN output, Vireo SN TSVs |
| **Writes** | `sn_cell_line_annotation.csv`, `sn_cell_line_annotation_post_qc.csv` |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |

### `07_sc_sn_umap.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 4C |
| **Reads** | FLAMES gene counts (all 4: ont/pb × sc/sn), CellRanger outputs (sc + sn), Vireo TSVs (all 6 combinations) |
| **Writes** | `RDS/sc_sn_filtered_so.rds` (list of filtered Seurat objects) |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |
| **Note** | Reads same raw inputs as `05/06_*_clustering_annotation.Rmd` independently — does not consume their output |

### `08_SC_identification_DE_analysis.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 4D, 4E–J |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/sc_DGE.obj.rds`, `RDS/bulk_identification.rds`, `RDS/bulk_DE.rds`, `lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv` |
| **Writes** | HTML report only |
| **Sources** | `scripts/Rfunctions.R` |

### `09_Mutation_analysis.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 5A, 5B, 5F, 5G–J |
| **Reads** | Variant call outputs (`lr_bulk/result/Mutation/`, `sr_bulk/result/Mutation/`, `lr_bulk/result/LongcallR/`) |
| **Writes** | HTML report only |

### `10_Mutation_analysis_downsample.Rmd`
| | |
|---|---|
| **Paper figures** | Fig 5C or 5D |
| **Reads** | Downsampled variant call outputs (`lr_bulk_ds/result/Mutation/`) |
| **Writes** | HTML report only |

---

## Intermediate RDS Files Summary

| RDS File | Produced by | Consumed by |
|---|---|---|
| `Tx2Gene.map.rds` | `Rscript/Tx2Gene.map.R` | `Bulk.DGElist.preprocessing.R`, `Sc.DGElist.preprocessing.R` |
| `bulk_DGE.obj.rds` | `Rscript/Bulk.DGElist.preprocessing.R` | `02_Bulk_Identification.Rmd`, `03_Bulk_quantification.Rmd`, `04_Bulk_DE_Summary.Rmd`, `08_SC_identification_DE_analysis.Rmd`, all `Sup_Bulk_*` Rmds |
| `sc_DGE.obj.rds` | `Rscript/Sc.DGElist.preprocessing.R` | `08_SC_identification_DE_analysis.Rmd` |
| `intronic_gene_and_exon_count.rds` | `Rscript/get_intronic_gene.R` | `02_Bulk_Identification.Rmd` |
| `bulk_identification.rds` | `02_Bulk_Identification.Rmd` | `08_SC_identification_DE_analysis.Rmd` |
| `bulk_DE.rds` | `04_Bulk_DE_Summary.Rmd` | `08_SC_identification_DE_analysis.Rmd` |
| `sc_sn_filtered_so.rds` | `07_sc_sn_umap.Rmd` | *(not consumed by any current Rmd)* |
