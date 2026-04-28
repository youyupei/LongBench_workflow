# Wrapper for converting miniQuant outputs to DGEList objects.
# Expected layout:
#   <dataset_dir>/<sample>/abundance.tsv

library(edgeR)

ABUNDANCE_FILENAME <- "abundance.tsv"
TRANSCRIPT_COLUMN <- "Transcript_id"
COUNT_COLUMN_CANDIDATES <- c(
  "Expected_num_long_reads",
  "Expected_num_short_read_pairs",
  "TPM"
)
COUNT_TYPE_LABELS <- c(
  "Expected_num_long_reads"       = "lr",
  "Expected_num_short_read_pairs" = "sr",
  "TPM"                          = "TPM"
)

read_tx2gene <- function(tx2gene_rds) {
  tx2gene_object <- readRDS(tx2gene_rds)
  tx2gene <- rbind(
    tx2gene_object$human.G.Tx.map[, c("tx_name", "gene_id")],
    tx2gene_object$sequins.G.Tx.map[, c("tx_name", "gene_id")],
    tx2gene_object$SIRV.G.Tx.map[, c("tx_name", "gene_id")]
  )
  tx2gene[!duplicated(tx2gene$tx_name), c("tx_name", "gene_id"), drop = FALSE]
}

read_sample_counts <- function(abundance_tsv) {
  abundance_df <- read.delim(
    abundance_tsv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(abundance_df) == 0) {
    stop("Empty file: ", abundance_tsv)
  }

  available_columns <- COUNT_COLUMN_CANDIDATES[COUNT_COLUMN_CANDIDATES %in% colnames(abundance_df)]
  if (length(available_columns) == 0) {
    stop(
      "Could not identify any count column in ", abundance_tsv,
      ". Columns: ", paste(colnames(abundance_df), collapse = ", ")
    )
  }

  transcript_ids <- abundance_df[[TRANSCRIPT_COLUMN]]

  count_list <- list()
  for (col in available_columns) {
    counts <- as.numeric(abundance_df[[col]])
    names(counts) <- transcript_ids
    count_list[[col]] <- counts
  }
  count_list
}


build_count_matrix <- function(dataset_dir) {
  sample_dirs <- list.dirs(dataset_dir, full.names = TRUE, recursive = FALSE)
  if (length(sample_dirs) == 0) {
    stop("No sample directories found under ", dataset_dir)
  }

  sample_names <- basename(sample_dirs)
  sample_order <- order(sample_names)
  sample_dirs <- sample_dirs[sample_order]
  sample_names <- sample_names[sample_order]

  all_sample_counts <- vector("list", length(sample_dirs))
  names(all_sample_counts) <- sample_names

  for (sample_index in seq_along(sample_dirs)) {
    abundance_tsv <- file.path(sample_dirs[sample_index], ABUNDANCE_FILENAME)
    if (!file.exists(abundance_tsv)) {
      stop("Missing ", ABUNDANCE_FILENAME, ": ", abundance_tsv)
    }
    cat("Reading ", abundance_tsv, "\n")
    all_sample_counts[[sample_index]] <- read_sample_counts(abundance_tsv)
  }

  count_types   <- names(all_sample_counts[[1]])
  transcript_ids <- names(all_sample_counts[[1]][[count_types[1]]])

  count_matrices <- list()
  for (count_type in count_types) {
    mat <- matrix(
      0,
      nrow = length(transcript_ids),
      ncol = length(sample_names),
      dimnames = list(transcript_ids, sample_names)
    )
    for (sample_index in seq_along(sample_names)) {
      mat[, sample_index] <- all_sample_counts[[sample_index]][[count_type]]
    }
    storage.mode(mat) <- "numeric"
    count_matrices[[count_type]] <- mat
  }

  count_matrices
}


make_dge <- function(count_matrix, bulk_meta) {
  dge <- DGEList(counts = count_matrix)
  rownames(dge$samples) <- colnames(dge$counts)

  sample_info <- dge$samples[, setdiff(colnames(dge$samples), "group"), drop = FALSE]
  sample_info <- merge(sample_info, bulk_meta, by = "row.names", all.x = TRUE)
  sample_info <- sample_info[, setdiff(colnames(sample_info), "Row.names"), drop = FALSE]
  dge$samples <- sample_info[match(colnames(dge$counts), sample_info$sample), , drop = FALSE]
  dge$genes <- data.frame(feature_id = rownames(dge$counts), row.names = rownames(dge$counts))

  dge
}


build_gene_counts <- function(tx_counts, tx2gene = NULL) {
  transcript_ids <- rownames(tx_counts)
  tx_map <- merge(
    data.frame(tx_name = transcript_ids, stringsAsFactors = FALSE),
    tx2gene,
    by = "tx_name",
    all.x = TRUE,
    sort = FALSE
  )
  tx_map <- tx_map[match(transcript_ids, tx_map$tx_name), , drop = FALSE]
  gene_ids <- tx_map$gene_id
  
  keep_rows <- !is.na(gene_ids) & gene_ids != ""
  gene_counts <- rowsum(tx_counts[keep_rows, , drop = FALSE], gene_ids[keep_rows])
  storage.mode(gene_counts) <- "numeric"

  gene_counts
}


create_dataset_bundle <- function(dataset_dir, dataset_name, bulk_meta, tx2gene = NULL) {
  tx_count_matrices <- build_count_matrix(dataset_dir)

  output_bundle <- list()
  for (count_type in names(tx_count_matrices)) {
    count_label <- COUNT_TYPE_LABELS[[count_type]]
    tx_counts   <- tx_count_matrices[[count_type]]
    gene_counts <- build_gene_counts(tx_counts, tx2gene)
    output_bundle[[paste0(dataset_name, ".", count_label, ".tx.dge")]]   <- make_dge(tx_counts, bulk_meta)
    output_bundle[[paste0(dataset_name, ".", count_label, ".gene.dge")]] <- make_dge(gene_counts, bulk_meta)
  }
  output_bundle
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript wrapper_MiniQuant.R ",
    "<dataset_dir> <dataset_name> <bulk_meta_csv> <output_rds> [tx2gene_rds]"
  )
}

dataset_dir <- args[1]
dataset_name <- args[2]
bulk_meta_csv <- args[3]
output_rds <- args[4]
tx2gene_rds <- if (length(args) >= 5) args[5] else NA_character_

# Test data
# Rscript /vast/projects/LongBench/analysis/workflow/sub_workflows/hybrid_bulk/scripts/wrapper_MiniQuant.R /vast/projects/LongBench/analysis/hybrid_bulk/result/miniquant_output/dRNA_bulk dRNA_bulk /vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt /vast/projects/LongBench/analysis/hybrid_bulk/result/miniquant_output/miniquant_bulk_dge.dRNA_bulk.rds /vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds
# dataset_dir <- '/vast/projects/LongBench/analysis/hybrid_bulk/result/miniquant_output/dRNA_bulk'
# dataset_name <- 'dRNA_bulk'
# bulk_meta_csv <- '/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt'
# output_rds <- '/vast/projects/LongBench/analysis/hybrid_bulk/result/miniquant_output/miniquant_bulk_dge.dRNA_bulk.rds'
# tx2gene_rds <- '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds'

if (!dir.exists(dataset_dir)) {
  stop("dataset_dir not found: ", dataset_dir)
}
if (!file.exists(bulk_meta_csv)) {
  stop("bulk metadata not found: ", bulk_meta_csv)
}

bulk_meta <- read.csv(bulk_meta_csv, stringsAsFactors = FALSE)
rownames(bulk_meta) <- bulk_meta$sample

output_bundle <- create_dataset_bundle(
  dataset_dir = dataset_dir,
  dataset_name = dataset_name,
  bulk_meta = bulk_meta,
  tx2gene = read_tx2gene(tx2gene_rds)
)

## Reformating the transcript names
# Extract and store original feature IDs before normalization
extract_and_normalize_ids <- function(dge) {
  # Store original full ID in $genes (including annotations)
  if (!"full_id" %in% colnames(dge$genes)) {
    dge$genes$full_id <- rownames(dge)
  }
  # Simplify rownames: remove pipe annotations and version numbers
  rownames(dge) <- sub("\\|.*", "", rownames(dge))      # Remove pipe-separated info
  dge
}
for (dge_name in names(output_bundle)) {
  output_bundle[[dge_name]] <- extract_and_normalize_ids(output_bundle[[dge_name]])
}

out_dir <- dirname(output_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(output_bundle, output_rds)
cat("Saved: ", output_rds, "\n")
cat("Entries: ", paste(names(output_bundle), collapse = ", "), "\n")