# This script reads the Quantification result from Oarfish (for long Read) and Salmon (for short read) and saves the DGE objects in RDS format
library(tximport) # install latest version from github
library(limma)
library(edgeR)
library(tidyverse)
library(conflicted)


# --------------------------------------------------------------
# SETUP
# --------------------------------------------------------------
params <- list(
    random_seed = 2024,
    ont_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/ont_bulk",
    pb_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/pb_bulk",
    ont_drnd_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/dRNA_bulk",
    ill_bulk_salmon_dir = "/vast/projects/LongBench/analysis/sr_bulk/result/salmon/salmon_quant",
    bulk_meta = "/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt",
    input_tx2gene.rds = "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds",
    output.rds = "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DGE.obj.rds"
)
source("/vast/projects/LongBench/analysis/workflow/scripts/Rfunctions.R")

conflicts_prefer(
    dplyr::select,
    dplyr::rename,
    dplyr::filter,
    base::intersect,
    base::setdiff
)


# --------------------------------------------------------------
# filtering functions for getting different types of gene and transcript
# --------------------------------------------------------------
is.sequins <- function(dge) {grepl("^R", rownames(dge))}
is.human <- function(dge) {grepl("^ENS", rownames(dge))}
is.sirv_or_ercc <- function(dge) {!is.human(dge) & !is.sequins(dge)} # this is the whole SIRV set4 (ERCC + SIRV E0 + Long SIRV)
is.sirv_E0 <- function(dge) {grepl("^SIRV\\d{3}$", rownames(dge))} 
is.sirv_all <- function(dge) {grepl("^SIRV", rownames(dge))} # including all SIRV E0 and Long SIRV
is.long_sirv <- function(dge) {grepl("^SIRV\\d{4}", rownames(dge))}
is.ercc <- function(dge) {!is.human(dge) & !is.sequins(dge) & !is.sirv_all(dge)}

# --------------------------------------------------------------
# Load the transcript level DGE objects
# --------------------------------------------------------------
bulk.meta <- read.csv(params$bulk_meta)
rownames(bulk.meta) <- bulk.meta$sample
## SR salmon
ill_bulk.dge <- get_dge_from_salmon(params$ill_bulk_salmon_dir, ".*/")
## LR oarfish
ont_bulk.dge.oarfish <- get_dge_from_oarfish(params$ont_bulk_oarfish_dir, ".*/")
pb_bulk.dge.oarfish <- get_dge_from_oarfish(params$pb_bulk_oarfish_dir, ".*/")
ont_drna.dge.oarfish <- get_dge_from_oarfish(params$ont_drnd_oarfish_dir, ".*/")


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
  ill_bulk = list(dir = params$ill_bulk_salmon_dir, file_name = "quant.sf"),
  ont_bulk = list(dir = params$ont_bulk_oarfish_dir, file_name = ".quant"),
  pb_bulk = list(dir = params$pb_bulk_oarfish_dir, file_name = ".quant"),
  drna_bulk = list(dir = params$ont_drnd_oarfish_dir, file_name = ".quant")
)

## Apply function to each directory to get all quant file paths
quant_paths <- lapply(oarfish_dirs, function(dir_info) {
  # Get the directory and file name from dir_info
  dir <- dir_info$dir
  file_name <- dir_info$file_name
  
  # List all subdirectories and get the quant file paths
  sub_dirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
  get_quant_paths(sub_dirs, file_name)
})

## Flatten the list of quant_paths into a named list
quant_paths <- lapply(quant_paths, unlist)

## Process each quant path for different bulk types
ill_bulk.gene.dge <- get_dge_from_txi(quant_paths$ill_bulk, ".*/", "salmon")
ont_bulk.gene.dge <- get_dge_from_txi(quant_paths$ont_bulk, ".*/", "oarfish")
pb_bulk.gene.dge <- get_dge_from_txi(quant_paths$pb_bulk, ".*/", "oarfish")
drna_bulk.gene.dge <- get_dge_from_txi(quant_paths$drna_bulk, ".*/", "oarfish")


# --------------------------------------------------------------
# OUTPUT
# --------------------------------------------------------------

saveRDS(
    list(
        ill_bulk.tx.dge = ill_bulk.dge,
        ont_bulk.tx.dge = ont_bulk.dge.oarfish,
        pb_bulk.tx.dge = pb_bulk.dge.oarfish,
        drna_bulk.tx.dge = ont_drna.dge.oarfish,
        ill_bulk.gene.dge = ill_bulk.gene.dge,
        ont_bulk.gene.dge = ont_bulk.gene.dge,
        pb_bulk.gene.dge = pb_bulk.gene.dge,
        drna_bulk.gene.dge = drna_bulk.gene.dge
    ), 
    file = params$output.rds
)

