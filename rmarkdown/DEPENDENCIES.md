# Rmarkdown Dependency Map

All active Rmd files and standalone R scripts in this directory. Excludes `old_analysis_backup/`.

---

## Dependency DAG

```
Reference files (GTF, SIRV, Sequins, metadata)
│
├──> Tx2Gene.map.R ──────────────────────────────> Tx2Gene.map.rds
│                                                         │
│   Workflow outputs (oarfish, salmon, pseudobulk) ───────┤
│                                                         │
│   ├──> Bulk.DGElist.preprocessing.R ──────────> bulk_DGE.obj.rds
│   ├──> Bulk.DGElist.preprocessing.IP_filtered.R > bulk_DGE.IP_filtered.obj.rds
│   └──> Sc.DGElist.preprocessing.R ────────────> sc_DGE.obj.rds
│
├──> get_intronic_gene.R ──────────────────────> intronic_gene_and_exon_count.rds
│
│   ┌─────────── bulk_DGE.obj.rds ──────────────────────────────────────────┐
│   │                                                                        │
│   ├──> Bulk_identification.Rmd ──────────────────────> bulk_identification.rds
│   │                                                            │
│   ├──> Bulk_DE_Summary.Rmd ──────────────────────────> bulk_DE.rds
│   │                                                            │
│   │                         bulk_DE.rds ──────────────────────┤
│   │                                   ├──> Bulk_DE_analysis_continue.Rmd
│   │                                   ├──> Bulk_DE_Summary_20M.Rmd
│   │                                   └──> SC_identification_DE_analysis.Rmd
│   │                                              (also needs bulk_identification.rds,
│   │                                               sc_DGE.obj.rds, pseudo_bulk CSV)
│   │
│   ├──> Bulk_identification_analysis.Rmd   (reads bulk_DGE + intronic_gene_count)
│   ├──> Bulk_quantification_analysis.Rmd   (reads bulk_DGE + IsoQuant/MiniQuant RDS)
│   ├──> Bulk_quantification_analysis_IsoQuant.Rmd
│   ├──> Bulk_quantification_cross_tool_comparison.Rmd
│   └──> Bulk_DE.spikeins.Rmd
│
│   bulk_DGE.IP_filtered.obj.rds
│   └──> Bulk_identification_IP_filtered.Rmd ─────> bulk_identification_IP_filtered.rds
│
│   sc_sn_umap.Rmd ──────────────────────────────> sc_sn_filtered_so.rds
│       (also needs FLAMES counts, CellRanger,              │
│        Vireo outputs, sc_long_preprocessing.R)            │
│                                                           └──> sc_sn_marker_analysis.Rmd ⚠
│
│   Standalone (no intermediate RDS produced):
│   ├──> QC_plot.Rmd               (reads AlignQC, Picard, Vireo, IP summary outputs)
│   ├──> Mutation_analysis.Rmd     (reads Mutation/LongcallR output dirs)
│   ├──> Bulk_GC_analysis_20M.Rmd  (reads 20M rarefaction outputs + Tx2Gene)
│   ├──> Bulk_identification_20M.Rmd
│   ├──> Bulk_identification_tool_comparison.Rmd ⚠
│   ├──> sc_clustering_annotation.Rmd
│   ├──> sn_clustering_annotation.Rmd
│   └──> dRNA_run_comparison.Rmd   (reads custom_analysis/dRNA_comparison outputs)
```

`⚠` = No Snakemake rule defined yet.

---

## Standalone R Scripts

### `Rscript/Tx2Gene.map.R`
| | |
|---|---|
| **Snakemake rule** | `r_tx2gene_map` |
| **Reads** | `gencode.v44.annotation.gtf`, SIRV GTF, Sequins GTF |
| **Writes** | `RDS/Tx2Gene.map.rds` |
| **Sources** | — |

### `Rscript/Bulk.DGElist.preprocessing.R`
| | |
|---|---|
| **Snakemake rule** | `r_get_bulk_DGE_objects` |
| **Reads** | `RDS/Tx2Gene.map.rds`, oarfish outputs (`lr_bulk/result/oarfish_cov_output/{ont,pb,dRNA}_bulk`), salmon outputs (`sr_bulk/result/salmon/salmon_quant`), `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/bulk_DGE.obj.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `Rscript/Bulk.DGElist.preprocessing.IP_filtered.R`
| | |
|---|---|
| **Snakemake rule** | `r_get_bulk_DGE_objects_IP_filtered` |
| **Reads** | `RDS/Tx2Gene.map.rds`, IP-filtered oarfish outputs (`lr_bulk/result/IP_filtered/oarfish_cov_output/{ont,pb,dRNA}_bulk`), salmon outputs, `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/bulk_DGE.IP_filtered.obj.rds` |
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
| **Reads** | `gencode.v44.annotation.gtf` (via `snakemake@input`) |
| **Writes** | `RDS/intronic_gene_and_exon_count.rds` (via `snakemake@output`) |
| **Sources** | — |

### `Rscript/Annotation_overlap.R`
| | |
|---|---|
| **Snakemake rule** | None |
| **Reads** | `gencode.v44.annotation.gtf` (hardcoded) |
| **Writes** | — |
| **Sources** | — |

### `versioned_knit.R`
| | |
|---|---|
| **Role** | Knitting helper — called by most Snakemake rules to render Rmd with versioned output filenames |
| **Reads** | YAML front matter of the target Rmd |
| **Writes** | Versioned HTML output, cache dir, figure dir |

---

## Rmd Files

### `QC_plot.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_QC_plots` |
| **Reads** | AlignQC outputs (`lr_bulk/result/qc/AlignQC/`, `sr_bulk/result/qc/AlignQC/`), Picard RNA metrics (`lr_bulk/result/qc/coverage/`), internal priming summaries (`lr_bulk/result/int_prim_analysis/`), Vireo outputs (`lr_sc_sn/result/vireo/`), barcode lists (`lr_sc_sn/result/misc/cell_line_bc_list/`), SC/SN AlignQC + Picard + IP outputs |
| **Writes** | HTML report; optional CSV raw-data exports to `params$fig.path` |
| **Sources** | — |
| **Workflow inputs** | lr_bulk QC, sr_bulk QC, lr_sc_sn QC sub-workflow outputs |

---

### `Bulk_identification.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_identification_analysis` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/intronic_gene_and_exon_count.rds`, `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/bulk_identification.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `Bulk_identification_IP_filtered.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_identification_analysis_IP_filtered` |
| **Reads** | `RDS/bulk_DGE.IP_filtered.obj.rds`, `RDS/intronic_gene_and_exon_count.rds`, `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/bulk_identification_IP_filtered.rds` |
| **Sources** | `scripts/Rfunctions.R` |

### `Bulk_identification_20M.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_identification_analysis_20M` |
| **Reads** | `RDS/Tx2Gene.map.rds`, 20M rarefaction oarfish/salmon outputs (`main_workflow/result/rarefraction_analysis/`), `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/bulk_identification_20M.rds` (to reports dir) |
| **Sources** | — |

### `Bulk_identification_tool_comparison.Rmd` ⚠
| | |
|---|---|
| **Snakemake rule** | **None — not in workflow yet** |
| **Reads** | `RDS/bulk_DGE.obj.rds` (oarfish), IsoQuant DGE RDS files (`lr_bulk/result/IsoQuant_output_no_novel/isoquant_bulk_dge.{ont,pb,dRNA}_bulk.rds`) |
| **Writes** | `RDS/bulk_identification_tool_comparison.rds` |
| **Sources** | `scripts/Rfunctions.R` |

---

### `Bulk_DE_Summary.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_DE_analysis` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `reference_files/Expected_DE.csv` |
| **Writes** | `RDS/bulk_DE.rds` |
| **Sources** | — |

### `Bulk_DE_analysis_continue.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_DE_analysis_continue` |
| **Reads** | `RDS/bulk_DE.rds` |
| **Writes** | HTML report only |
| **Sources** | — |

### `Bulk_DE_Summary_20M.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_DE_analysis_subsample_20M` |
| **Reads** | `RDS/bulk_DE.rds`, 20M rarefaction oarfish/salmon outputs, Sequins/SIRV references, `illumina_bulk/metadata.txt` |
| **Writes** | `RDS/Bulk_identification_DE_M20.rds` |
| **Sources** | — |

### `Bulk_DE.spikeins.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_DE_spikeins` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/Tx2Gene.map.rds`, Sequins TSV, SIRV CSV, `illumina_bulk/metadata.txt` |
| **Writes** | HTML report only |
| **Sources** | — |

---

### `Bulk_quantification_analysis.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_quantification_analysis` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/Tx2Gene.map.rds`, IsoQuant DGE RDS (`lr_bulk/result/IsoQuant_output_no_novel_bc/`), MiniQuant DGE RDS (`hybrid_bulk/result/miniquant_output/`), Sequins TSV, SIRV CSV |
| **Writes** | HTML report only |
| **Sources** | — |
| **Workflow inputs** | hybrid_bulk sub-workflow outputs |

### `Bulk_quantification_analysis_IsoQuant.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_quantification_analysis_isoquant` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/Tx2Gene.map.rds`, IsoQuant DGE RDS (`lr_bulk/result/IsoQuant_output_no_novel_bc/`), Sequins TSV, SIRV CSV |
| **Writes** | HTML report only |
| **Sources** | — |
| **Workflow inputs** | lr_bulk IsoQuant sub-workflow outputs |

### `Bulk_quantification_cross_tool_comparison.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_quantification_cross_tool_comparison` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, IsoQuant DGE RDS, MiniQuant DGE RDS, Sequins TSV |
| **Writes** | HTML report only |
| **Sources** | — |
| **Workflow inputs** | lr_bulk (IsoQuant) and hybrid_bulk (MiniQuant) outputs |

---

### `Bulk_GC_analysis_20M.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_bulk_GC_analysis_20M` |
| **Reads** | `RDS/Tx2Gene.map.rds`, `RDS/intronic_gene_and_exon_count.rds`, 20M rarefaction oarfish/salmon outputs |
| **Writes** | HTML report only |
| **Sources** | — |

---

### `Mutation_analysis.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_mutation_analysis` |
| **Reads** | Mutation/variant call outputs (`lr_bulk/result/Mutation/`, `sr_bulk/result/Mutation/`, `lr_bulk/result/LongcallR/`), allele-specific summary dir (`variant_allele_specific_analysis/`) |
| **Writes** | HTML report only |
| **Sources** | — |

---

### `sc_clustering_annotation.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_sc_clustering_annotation` |
| **Reads** | FLAMES gene counts (`lr_sc_sn/result/flames_out/ont_sc/gene_count.csv`, `pb_sc/gene_count.csv`), CellRanger output (`illumina_sc_sn/result/cellranger/ill_sc/`), Vireo TSVs (`lr_sc_sn/result/vireo/{ont,pb,ill}_sc_vireo/donor_ids.tsv`) |
| **Writes** | `sc_cell_line_annotation.csv`, `sc_cell_line_annotation_post_qc.csv` (to figures dir) |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |
| **Workflow inputs** | lr_sc_sn sub-workflow outputs |

### `sn_clustering_annotation.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_sn_clustering_annotation` |
| **Reads** | FLAMES gene counts (`lr_sc_sn/result/flames_out/ont_sn/gene_count.csv`, `pb_sn/gene_count.csv`), CellRanger output (`illumina_sc_sn/result/cellranger/ill_sn/`), Vireo TSVs (`lr_sc_sn/result/vireo/{ont,pb,ill}_sn_vireo/donor_ids.tsv`) |
| **Writes** | `sn_cell_line_annotation.csv`, `sn_cell_line_annotation_post_qc.csv` (to figures dir) |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |
| **Workflow inputs** | lr_sc_sn sub-workflow outputs |

### `sc_sn_umap.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_sc_sn_umap` |
| **Reads** | FLAMES gene counts (all 4: ont/pb × sc/sn), CellRanger outputs (sc + sn), Vireo TSVs (all 6 donors), `sc_cell_line_annotation_post_qc.csv`, `sn_cell_line_annotation_post_qc.csv` |
| **Writes** | `RDS/sc_sn_filtered_so.rds` (list of 6 Seurat objects) |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |
| **Note** | Snakemake rule does NOT declare the annotation CSVs as inputs — implicit dependency on `sc/sn_clustering_annotation.Rmd` |

### `RDS/sc_sn_marker_analysis.Rmd` ⚠
| | |
|---|---|
| **Snakemake rule** | **None — not in workflow yet** |
| **Reads** | `RDS/sc_sn_filtered_so.rds` |
| **Writes** | HTML report only |
| **Sources** | `sub_workflows/lr_sc_sn/rmarkdown/sc_long_preprocessing.R` |

### `SC_identification_DE_analysis.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_sc_pseudo_bulk_analysis` |
| **Reads** | `RDS/bulk_DGE.obj.rds`, `RDS/sc_DGE.obj.rds`, `RDS/bulk_identification.rds`, `RDS/bulk_DE.rds`, `lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv`, `reference_files/GRCh38/gencode.v44.annotation.gtf`, `illumina_bulk/metadata.txt`, `reference_files/Expected_DE.csv` |
| **Writes** | HTML report only |
| **Sources** | `scripts/Rfunctions.R` |

---

### `dRNA_run_comparison.Rmd`
| | |
|---|---|
| **Snakemake rule** | `rmd_dRNA_run_comparison` |
| **Reads** | dRNA comparison oarfish outputs (`custom_analysis/dRNA_comparison` result dir: `oarfish_cov_output/{dRNA_bulk_run1,dRNA_bulk_run2}/`, read count files) |
| **Writes** | HTML report only |
| **Sources** | — |
| **Workflow inputs** | `custom_analysis/dRNA_comparison` sub-workflow outputs |

---

## Intermediate RDS Files Summary

| RDS File | Produced by | Consumed by |
|---|---|---|
| `Tx2Gene.map.rds` | `Tx2Gene.map.R` | `Bulk.DGElist.preprocessing.R`, `Bulk.DGElist.preprocessing.IP_filtered.R`, `Sc.DGElist.preprocessing.R`, `Bulk_DE.spikeins.Rmd`, `Bulk_quantification_analysis.Rmd`, `Bulk_quantification_analysis_IsoQuant.Rmd`, `Bulk_identification_20M.Rmd`, `Bulk_GC_analysis_20M.Rmd` |
| `bulk_DGE.obj.rds` | `Bulk.DGElist.preprocessing.R` | `Bulk_identification.Rmd`, `Bulk_DE_Summary.Rmd`, `Bulk_DE.spikeins.Rmd`, `Bulk_quantification_analysis.Rmd`, `Bulk_quantification_analysis_IsoQuant.Rmd`, `Bulk_quantification_cross_tool_comparison.Rmd`, `Bulk_identification_tool_comparison.Rmd`, `SC_identification_DE_analysis.Rmd` |
| `bulk_DGE.IP_filtered.obj.rds` | `Bulk.DGElist.preprocessing.IP_filtered.R` | `Bulk_identification_IP_filtered.Rmd` |
| `sc_DGE.obj.rds` | `Sc.DGElist.preprocessing.R` | `SC_identification_DE_analysis.Rmd` |
| `intronic_gene_and_exon_count.rds` | `get_intronic_gene.R` | `Bulk_identification.Rmd`, `Bulk_identification_IP_filtered.Rmd`, `Bulk_GC_analysis_20M.Rmd` |
| `bulk_identification.rds` | `Bulk_identification.Rmd` | `SC_identification_DE_analysis.Rmd` |
| `bulk_identification_IP_filtered.rds` | `Bulk_identification_IP_filtered.Rmd` | *(standalone, not consumed downstream)* |
| `bulk_DE.rds` | `Bulk_DE_Summary.Rmd` | `Bulk_DE_analysis_continue.Rmd`, `Bulk_DE_Summary_20M.Rmd`, `SC_identification_DE_analysis.Rmd` |
| `sc_sn_filtered_so.rds` | `sc_sn_umap.Rmd` | `sc_sn_marker_analysis.Rmd` |

## Missing Snakemake Rules

| File | What's needed |
|---|---|
| `RDS/sc_sn_marker_analysis.Rmd` | Rule with `RDS/sc_sn_filtered_so.rds` as declared input |
| `Bulk_identification_tool_comparison.Rmd` | Rule with `bulk_DGE.obj.rds` + IsoQuant RDS files as inputs, outputs `RDS/bulk_identification_tool_comparison.rds` |
