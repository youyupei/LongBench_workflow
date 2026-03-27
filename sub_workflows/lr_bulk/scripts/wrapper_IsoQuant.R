# This is a wrapper script for reformating the output of IsoQuant to be compatible with the downstream analysis. 
# Includes:
# 1. Getting gene-level and transcript-level quantification results from IsoQuant output.
# 2. Reading the metadata file to get the sample information and adding it to the output files.
# 3. Saving the reformatted output as DGElist objects for downstream analysis.


# Load required libraries
library(edgeR)
library(tidyverse)

# Setup meta data if not already defined
bulk.meta <- read.csv("/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt")
rownames(bulk.meta) <- bulk.meta$sample

# =====================================================================
# Function to read IsoQuant output and convert to DGEList
# =====================================================================

#' Read IsoQuant gene/transcript counts and create DGEList
#'
#' @param isoquant_dir Directory containing IsoQuant output directories (one per sample)
#' @param level "gene" or "transcript" - which quantification level to read
#' @param sample_prefix_regex Regular expression to extract sample names from directory paths
#'
#' @return DGEList object with counts and sample metadata
#'
#' @details
#' IsoQuant outputs count tables under each sample's OUT directory.
#' This function:
#' 1. Reads the specified count file from each sample directory
#' 2. Extracts feature IDs and raw count columns
#' 3. Merges counts across samples
#' 4. Adds sample metadata from bulk.meta
#' 5. Returns a properly formatted DGEList object
get_dge_from_isoquant <- function(isoquant_dir, level = "gene", sample_prefix_regex = ".*/") {
  
	# Determine which file to read based on level
	if (level == "gene") {
		count_file <- file.path("OUT", "OUT.gene_counts.tsv")
	} else if (level == "transcript") {
		count_file <- file.path("OUT", "OUT.transcript_counts.tsv")
	} else {
		stop("level must be either 'gene' or 'transcript'")
	}
  
	# Get all sample directories
	sample_dirs <- list.dirs(isoquant_dir, full.names = TRUE, recursive = FALSE)
	if (length(sample_dirs) == 0) {
		stop("No sample directories found in: ", isoquant_dir)
	}
  
	# Extract sample names and sort
	sample_names <- basename(sample_dirs)
	sample_dirs <- sample_dirs[order(sample_names)]
	sample_names <- sample_names[order(sample_names)]
  
	counts_list <- list()
	gene_info <- NULL
  
	# Read abundance files from each sample
	for (i in seq_along(sample_dirs)) {
		sample_dir <- sample_dirs[i]

		count_path <- file.path(sample_dir, count_file)
		if (!file.exists(count_path)) {
			stop("Could not find ", count_file, " in: ", sample_dir)
		}

		cat("Reading ", count_path, " ... ")

		count_df <- read.delim(count_path, stringsAsFactors = FALSE, check.names = FALSE)
		feature_col <- colnames(count_df)[1]
		count_col <- grep("^count$", colnames(count_df), ignore.case = TRUE, value = TRUE)[1]

		if (is.na(count_col)) {
			stop("Could not find a 'count' column in ", count_path,
					 ". Available columns: ", paste(colnames(count_df), collapse = ", "))
		}

		rownames(count_df) <- count_df[[feature_col]]
		sample_counts <- as.numeric(count_df[[count_col]])
		names(sample_counts) <- count_df[[feature_col]]
		counts_list[[sample_names[i]]] <- sample_counts

		if (is.null(gene_info)) {
			gene_info <- data.frame(
				feature_id = count_df[[feature_col]],
				row.names = count_df[[feature_col]],
				stringsAsFactors = FALSE
			)
		}

		cat("done\n")
	}
  
	# Convert list to matrix
	all_ids <- unique(unlist(lapply(counts_list, names)))
	counts_matrix <- matrix(0, nrow = length(all_ids), ncol = length(counts_list))
	rownames(counts_matrix) <- all_ids
	colnames(counts_matrix) <- names(counts_list)
  
	for (i in seq_along(counts_list)) {
		counts_matrix[names(counts_list[[i]]), i] <- counts_list[[i]]
	}
  
	storage.mode(counts_matrix) <- "numeric"

	# Prepare gene annotation
	if (!is.null(gene_info)) {
		gene_data <- gene_info[rownames(counts_matrix), , drop = FALSE]
	} else {
		gene_data <- data.frame(row.names = rownames(counts_matrix))
	}
  
	# Create DGEList
	dge <- DGEList(counts = counts_matrix, genes = gene_data)
  
	# Add sample metadata
	# Extract sample names from column names (remove prefix if specified)
	sample_names_clean <- colnames(dge$counts) %>% 
		sub(sample_prefix_regex, '', .)
  
	# Update sample names
	rownames(dge$samples) <- sample_names_clean
	colnames(dge$counts) <- sample_names_clean
  
	# Merge with metadata (same style as existing Oarfish helper)
	dge$samples <- merge(
		dge$samples %>% select(-group),
		bulk.meta,
		by = "row.names",
		all.x = TRUE
	) %>% select(-Row.names)

	# Ensure sample order matches counts
	dge$samples <- dge$samples[match(colnames(dge$counts), dge$samples$sample), ]
  
	return(dge)
}


# =====================================================================
# Function to create combined DGEList for both gene and transcript levels
# =====================================================================

#' Process IsoQuant results and create gene and transcript level DGELists
#'
#' @param isoquant_dir Directory containing IsoQuant output
#' @param sample_prefix_regex Regular expression for sample name extraction
#'
#' @return List containing $gene_dge and $tx_dge DGEList objects
#'
create_isoquant_deglists <- function(isoquant_dir, 
																		 sample_prefix_regex = ".*/") {
  
	cat("Processing IsoQuant results from: ", isoquant_dir, "\n")
  
	# Read gene-level quantification
	cat("\n=== Gene-level quantification ===\n")
	gene_dge <- get_dge_from_isoquant(isoquant_dir, 
																		 level = "gene",
																		 sample_prefix_regex = sample_prefix_regex)
  
	# Read transcript-level quantification
	cat("\n=== Transcript-level quantification ===\n")
	tx_dge <- get_dge_from_isoquant(isoquant_dir, 
																	 level = "transcript",
																	 sample_prefix_regex = sample_prefix_regex)
  
	cat("\n=== Summary ===\n")
	cat("Gene-level DGEList:\n")
	cat("  Genes: ", nrow(gene_dge), "\n")
	cat("  Samples: ", ncol(gene_dge), "\n")
	cat("  Sample names: ", paste(colnames(gene_dge), collapse = ", "), "\n")
  
	cat("\nTranscript-level DGEList:\n")
	cat("  Transcripts: ", nrow(tx_dge), "\n")
	cat("  Samples: ", ncol(tx_dge), "\n")
  
	return(list(gene_dge = gene_dge, tx_dge = tx_dge))
}


#' Create a bulk-style DGE object list from IsoQuant outputs
#'
#' @param isoquant_root_dir Either a dataset directory (e.g. .../IsoQuant_output/ont_bulk)
#'        or the top-level IsoQuant directory containing multiple datasets
#'        (e.g. .../IsoQuant_output)
#' @param sample_prefix_regex Regular expression for sample name extraction
#'
#' @return Named list with keys like ont_bulk.gene.dge / ont_bulk.tx.dge
#'
create_isoquant_bulk_dge_bundle <- function(isoquant_root_dir,
																		sample_prefix_regex = ".*/") {

	dataset_dirs <- list.dirs(isoquant_root_dir, full.names = TRUE, recursive = FALSE)
	if (length(dataset_dirs) == 0) {
		stop("No subdirectories found in: ", isoquant_root_dir)
	}

	# Detect whether input is already a dataset dir (contains sample dirs with OUT files)
	is_single_dataset <- any(file.exists(file.path(dataset_dirs, "OUT", "OUT.gene_counts.tsv")))

	if (is_single_dataset) {
		# isoquant_root_dir itself is a dataset directory
		dataset_map <- setNames(list(isoquant_root_dir), basename(isoquant_root_dir))
	} else {
		# isoquant_root_dir is a top-level directory containing multiple datasets
		dataset_map <- setNames(as.list(dataset_dirs), basename(dataset_dirs))
	}

	out <- list()
	for (dataset_name in names(dataset_map)) {
		dataset_dir <- dataset_map[[dataset_name]]
		cat("\n==============================\n")
		cat("Dataset: ", dataset_name, "\n")
		cat("Path: ", dataset_dir, "\n")
		cat("==============================\n")

		dge_pair <- create_isoquant_deglists(dataset_dir, sample_prefix_regex)
		out[[paste0(dataset_name, ".gene.dge")]] <- dge_pair$gene_dge
		out[[paste0(dataset_name, ".tx.dge")]] <- dge_pair$tx_dge
	}

	return(out)
}


# =====================================================================
# Main execution (if run as script)
# =====================================================================

if (!interactive()) {
  
	# Parse command line arguments
	args <- commandArgs(trailingOnly = TRUE)
  
	if (length(args) < 2) {
		cat("Usage: Rscript wrapper_IsoQuant.R <isoquant_dir> <output_prefix> [sample_regex]\n")
		cat("  <isoquant_dir>: Directory containing IsoQuant output subdirectories\n")
		cat("  <output_prefix>: Prefix for output RDS files\n")
		cat("  [sample_regex]: Optional regex to extract sample names (default: '.*/')\n")
		quit(status = 1)
	}
  
	isoquant_dir <- args[1]
	output_prefix <- args[2]
	sample_regex <- if (length(args) >= 3) args[3] else ".*/"
  
	# Validate input directory
	if (!dir.exists(isoquant_dir)) {
		cat("Error: Directory does not exist: ", isoquant_dir, "\n")
		quit(status = 1)
	}
  
	# Create output directory if needed
	output_dir <- dirname(output_prefix)
	if (!dir.exists(output_dir)) {
		dir.create(output_dir, recursive = TRUE)
	}
  
	# Process IsoQuant results
	# If isoquant_dir points to a dataset folder, output is one pair (gene + tx).
	# If isoquant_dir points to top-level IsoQuant_output, output includes all datasets.
	dge_list <- create_isoquant_bulk_dge_bundle(isoquant_dir, sample_regex)
  
	# Save results
	cat("\nSaving results...\n")
  
	list_output <- paste0(output_prefix, ".rds")
  
	saveRDS(dge_list, file = list_output)
	cat("Saved complete DGE bundle to: ", list_output, "\n")
	cat("Bundle entries: ", paste(names(dge_list), collapse = ", "), "\n")
  
	cat("\nDone!\n")
}