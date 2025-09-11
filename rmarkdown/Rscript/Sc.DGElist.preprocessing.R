# params:
#   random_seed: 2024
#   cache_dir: NULL
#   ont_sc_dir: "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sc"
#   ont_sn_dir: "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sn"
#   pb_sc_dir: "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sc"
#   pb_sn_dir: "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sn"
#   human_gtf: '/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf'
#   bulk_meta: "/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt"
#   bulk_ident_rds: "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification.rds"
#   bulk_de_rds: "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DE.rds"
#   expected_de_table: '/vast/projects/LongBench/reference_files/Expected_DE.csv'
#   fig.path: NULL

# This script reads the Quantification Pseudobulk result from Oarfish ( long Read)  saves the DGE objects in RDS format
library(tximport) # install latest version from github
library(limma)
library(edgeR)
library(tidyverse)
library(conflicted)
source("/vast/projects/LongBench/analysis/workflow/scripts/Rfunctions.R")

conflicts_prefer(
    dplyr::select,
    dplyr::rename,
    dplyr::filter,
    base::intersect,
    base::setdiff
)

# --------------------------------------------------------------
# SETUP
# --------------------------------------------------------------
params <- list(
  random_seed = 2024,
  cache_dir = NULL,
  ont_sc_dir = "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sc",
  ont_sn_dir = "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sn",
  pb_sc_dir = "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sc",
  pb_sn_dir = "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sn",
  input_tx2gene.rds = "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds",
  output.rds = "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/sc_DGE.obj.rds"
)

# SR salmon
# ill_bulk.dge <- get_dge_from_salmon(params$ill_bulk_salmon_dir, ".*/")

# --------------------------------------------------------------
# Load the transcript level DGE objects
# --------------------------------------------------------------
ont_sc.dge <- get_dge_from_oarfish(params$ont_sc_dir, ".*/")
ont_sn.dge <- get_dge_from_oarfish(params$ont_sn_dir, ".*/")
pb_sc.dge <- get_dge_from_oarfish(params$pb_sc_dir, ".*/")
pb_sn.dge <- get_dge_from_oarfish(params$pb_sn_dir, ".*/")

# --------------------------------------------------------------
# Load the gene level DGE objects
# --------------------------------------------------------------

## Extract tx_name and gene_id columns from human , sequins and SIRV
list2env(readRDS(params$input_tx2gene.rds), .GlobalEnv)
human_tx2gene <- human.G.Tx.map[, c("tx_name", "gene_id")]
sequins_tx2gene <- sequins.G.Tx.map[, c("tx_name", "gene_id")]
SIRV_tx2gene <- SIRV.G.Tx.map[, c("tx_name", "gene_id")]
combined_tx2gene <- rbind(human_tx2gene, sequins_tx2gene, SIRV_tx2gene)

## Get quant.sf file paths
get_quant_paths <- function(dirs, file_name) {
  # Generate file paths for quantification files
  quant_paths <- file.path(dirs, file_name)
  return(quant_paths)
}
## Define directories with their corresponding expected file names
oarfish_dirs <- list(
    ont_sc = list(dir = params$ont_sc_dir, file_name = ".quant"),
    ont_sn = list(dir = params$ont_sn_dir, file_name = ".quant"),
    pb_sc = list(dir = params$pb_sc_dir, file_name = ".quant"),
    pb_sn = list(dir = params$pb_sn_dir, file_name = ".quant")
)

# Apply function to each directory to get all quant file paths
quant_paths <- lapply(oarfish_dirs, function(dir_info) {
  # Get the directory and file name from dir_info
  dir <- dir_info$dir
  file_name <- dir_info$file_name
  
  # List all subdirectories and get the quant file paths
  sub_dirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
  get_quant_paths(sub_dirs, file_name)
})

# Flatten the list of quant_paths into a named list
quant_paths <- lapply(quant_paths, unlist)

# Function to import DGE from tximport
get_dge_from_txi <- function(quant_dir, sample_prefix_regex, type) {
  # Import data using tximport
  txi <- tximport(quant_dir,
                  type = type,
                  tx2gene = human_tx2gene,
                  ignoreAfterBar = TRUE,
                  countsFromAbundance = "no") # raw counts
  counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)
  sample_names <- gsub(sample_prefix_regex, '', basename(dirname(quant_dir)))
  colnames(counts) <- sample_names
  rst.dge <- DGEList(counts = counts)
  rst.dge$genes <- data.frame(Length=txi$length %>% rowMeans())
  rownames(rst.dge$samples) <- sub(sample_prefix_regex, '', rownames(rst.dge$samples))
  colnames(rst.dge$counts) <- sub(sample_prefix_regex, "", colnames(rst.dge$counts))
  # colnames(rst.dge$gene$Length) <- sub(sample_prefix_regex, "", colnames(rst.dge$gene$Length))
  
  rst.dge$samples <- merge(rst.dge$samples %>% select(-group), bulk.meta, by = "row.names", all.x = TRUE) %>% select(-Row.names)
  rst.dge$samples <- rst.dge$samples[match(colnames(rst.dge$counts), rst.dge$samples$sample),]
  return(rst.dge)
}

# Process each quant path for different bulk types
ont_sc.gene.dge <- get_dge_from_txi(quant_paths$ont_sc, ".*/", "oarfish")
ont_sn.gene.dge <- get_dge_from_txi(quant_paths$ont_sn, ".*/", "oarfish")
pb_sc.gene.dge <- get_dge_from_txi(quant_paths$pb_sc, ".*/", "oarfish")
pb_sn.gene.dge <- get_dge_from_txi(quant_paths$pb_sn, ".*/", "oarfish")

# --------------------------------------------------------------
# OUTPUT
# --------------------------------------------------------------

saveRDS(
    list(
        ont_sc.tx.dge = ont_sc.dge,
        ont_sn.tx.dge = ont_sn.dge,
        pb_sc.tx.dge = pb_sc.dge,
        pb_sn.tx.dge = pb_sn.dge,
        ont_sc.gene.dge = ont_sc.gene.dge,
        ont_sn.gene.dge = ont_sn.gene.dge,
        pb_sc.gene.dge = pb_sc.gene.dge,
        pb_sn.gene.dge = pb_sn.gene.dge
    ), 
    file = params$output.rds
)

