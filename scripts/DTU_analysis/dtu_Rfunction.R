

my_run_drimseq <- function(counts, tx2gene, pd, id_col = NULL, cond_col, cond_levels = NULL, 
                           filtering_strategy = "bulk", add_pseudocount = FALSE, 
                           BPPARAM = BiocParallel::SerialParam(), force_dense = TRUE, 
                           subset_feature = NULL, subset_sample = NULL, 
                           carry_over_metadata = TRUE, filter_only = FALSE, ...) {
  if (methods::is(counts, "Seurat")) {
    assertthat::assert_that(requireNamespace("Seurat", quietly = TRUE), 
                            msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
    assertthat::assert_that(utils::packageVersion("Seurat") >= "3.0.0", 
                            msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
    assertthat::assert_that(counts@version >= "3.0.0", 
                            msg = "The provided 'counts' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
    if (is.vector(tx2gene)) {
      meta_df <- counts[[counts@active.assay]]@meta.features
      assertthat::assert_that(ncol(meta_df) > 1, 
                              msg = "No feature-level meta data in active assay of seurat object. Was 'tx2gene' provided in 'combine_to_matrix()'?\nAlternatively provide a real tx2gene dataframe.")
      assertthat::assert_that(all(tx2gene %in% colnames(meta_df)), 
                              msg = "Not all provided 'tx2gene' colnames are present in the feature-level meta data.")
      tx2gene <- move_columns_to_front(meta_df, tx2gene)
    }
    counts <- Seurat::GetAssayData(counts)
  }
  assertthat::assert_that(methods::is(counts, "matrix") | methods::is(counts, "sparseMatrix"), 
                          msg = "Counts must be a (sparse) matrix.")
  assertthat::assert_that(methods::is(tx2gene, "data.frame"), 
                          msg = "Tx2gene must be a data frame.")
  assertthat::assert_that(methods::is(pd, "data.frame"), 
                          msg = "pd must be a data frame.")
  assertthat::assert_that(ncol(tx2gene) > 1, 
                          msg = "'tx2gene' should at least have two columns [feature | gene --- in that order].")
  assertthat::assert_that(is.null(id_col) || (is.character(id_col) && length(id_col) == 1 && id_col %in% colnames(pd)), 
                          msg = "id_col should be a single column name of pd or NULL.")
  assertthat::assert_that((is.character(cond_col) && length(cond_col) == 1 && cond_col %in% colnames(pd)), 
                          msg = paste0("Could not find ", cond_col, " in column names of pd."))
  assertthat::assert_that(is.null(cond_levels) || length(cond_levels) == 2, 
                          msg = "'cond_levels' should be of length two or NULL.")
  # Allow 'none' as a valid filtering_strategy
  assertthat::assert_that(filtering_strategy %in% c("bulk", "sc", "own", "none"), 
                          msg = "Please select a valid filtering strategy ('bulk', 'sc', 'own', or 'none').")
  assertthat::assert_that(is.logical(add_pseudocount), 
                          msg = "`add_pseudocount` must be `TRUE` or `FALSE`.")
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), 
                          msg = "Please provide a valid BiocParallelParam object.")
  assertthat::assert_that(is.logical(force_dense), 
                          msg = "`force_dense` must be `TRUE` or `FALSE`.")
  assertthat::assert_that(is.null(subset_feature) | length(subset_feature) > 0, 
                          msg = "`subset_feature` must be `NULL` or of length>=1.")
  assertthat::assert_that(is.null(subset_sample) | length(subset_sample) > 0, 
                          msg = "`subset_sample` must be `NULL` or of length>=1.")
  assertthat::assert_that(is.logical(carry_over_metadata), 
                          msg = "`carry_over_metadata` must be `TRUE` or `FALSE`.")
  assertthat::assert_that((length(intersect(rownames(counts), tx2gene[[1]])) > 0), 
                          msg = paste0(
                            "The provided counts names and tx2gene names do not match.\n\tCounts names: ",
                            paste0(rownames(utils::head(counts, n = 5)), collapse = ", "), "\n\tTx2gene names: ", 
                            paste0(utils::head(tx2gene, n = 5)[[1]], collapse = ", ")
                          ))
  assertthat::assert_that(is.logical(filter_only), 
                          msg = "The 'filter_only' paramter must be TRUE or FALSE.")
  
  if (!is.null(subset_feature) | !is.null(subset_sample)) {
    if (is.null(subset_feature)) {
      subset_feature <- TRUE
    }
    if (is.null(subset_sample)) {
      subset_sample <- TRUE
    }
    if (is.character(subset_feature)) {
      assertthat::assert_that(all(subset_feature %in% rownames(counts)), 
                              msg = "Invalid 'subset_feature' names provided.")
    }
    if (is.character(subset_sample)) {
      assertthat::assert_that(all(subset_sample %in% colnames(counts)), 
                              msg = "Invalid 'subset_sample' names provided.")
    }
    counts <- counts[subset_feature, subset_sample, drop = FALSE]
  }
  
  assertthat::assert_that(all(rownames(counts) %in% tx2gene[[1]]), 
                          msg = "Not all rownames of the counts are present in the first tx2gene column. You may want to reorder the tx2gene columns with 'move_columns_to_front()' or use 'subset_feature' to subset the counts.")
  
  tx2gene <- rapply(tx2gene, as.character, classes = "factor", how = "replace")
  tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]), ]
  assertthat::assert_that(nrow(tx2gene) == nrow(counts))
  message("Using tx2gene columns:\n\t", colnames(tx2gene)[[1]], " ---> 'feature_id'\n\t", colnames(tx2gene)[[2]], " ---> 'gene_id'")
  colnames(tx2gene)[c(1, 2)] <- c("feature_id", "gene_id")
  colnames(tx2gene) <- make.names(colnames(tx2gene), unique = TRUE)
  
  if (is.null(cond_levels)) {
    cond_levels <- unique(pd[[cond_col]])
  }
  assertthat::assert_that(length(cond_levels) == 2, 
                          msg = "More than two levels found in 'cond_col'. Please specify the two levels you want to compare in 'cond_levels'.")
  message("\nComparing in '", cond_col, "': '", cond_levels[1], "' vs '", cond_levels[2], "'")
  
  if (is.null(id_col)) {
    assertthat::assert_that(all(rownames(pd) %in% colnames(counts)), 
                            msg = "Provided id_col does not match with sample names in counts.")
    samp <- data.frame(
      "sample_id" = rownames(pd), "condition" = as.character(pd[[cond_col]]),
      pd[, -c(which(colnames(pd) == cond_col)), drop = FALSE],
      row.names = NULL, stringsAsFactors = FALSE
    )
  } else {
    samp <- data.frame(
      "sample_id" = pd[[id_col]], "condition" = as.character(pd[[cond_col]]),
      pd[, -c(which(colnames(pd) %in% c(id_col, cond_col))), drop = FALSE],
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  samp$condition <- factor(samp$condition, levels = cond_levels)
  samp <- samp[samp$sample_id %in% colnames(counts), , drop = FALSE]
  counts <- counts[, samp$sample_id, drop = FALSE]
  
  # exclude samples not in comparison
  exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
  if (length(exclude) != 0) {
    message("Excluding ", ifelse(length(exclude) < 10, paste(exclude, collapse = " "), paste(length(exclude), "cells/samples")), " for this comparison!")
    samp <- samp[!is.na(samp$condition), ]
    counts <- counts[, !(colnames(counts) %in% exclude), drop = FALSE]
  }
  message("\nProceed with cells/samples: ", paste0(utils::capture.output(table(samp$condition)), collapse = "\n"))
  assertthat::assert_that(length(levels(samp$condition)) == 2, 
                          msg = "No two sample groups left for comparison. Aborting!")
  assertthat::assert_that(all(table(samp$condition) > 0), 
                          msg = "No sample in each group left for comparison. Aborting!")
  
  message("\nFiltering...\n")
  filter_opt_list <- list(
    "min_samps_gene_expr" = 0,
    "min_samps_feature_expr" = 0, 
    "min_samps_feature_prop" = 0,
    "min_gene_expr" = 0, 
    "min_feature_expr" = 0, 
    "min_feature_prop" = 0,
    "run_gene_twice" = FALSE
  )
  
  if (filtering_strategy == "none") {
    message("Skipping filtering step as filtering_strategy is set to 'none'.")
    # Leave counts unchanged.
  } else {
    switch(filtering_strategy,
           sc = {
             smallest_group <- min(table(samp$condition)) * 0.05
             filter_opt_list <- utils::modifyList(filter_opt_list, list(
               "min_samps_feature_prop" = smallest_group,
               "min_feature_prop" = 0.05, "run_gene_twice" = TRUE
             ))
           },
           bulk = {
             smallest_group <- min(table(samp$condition)) * 0.5
             filter_opt_list <- utils::modifyList(
               filter_opt_list,
               list(
                 "min_samps_gene_expr" = smallest_group,
                 "min_gene_expr" = 5,
                 "min_samps_feature_prop" = smallest_group,
                 "min_feature_prop" = 0.05, "run_gene_twice" = TRUE
               )
             )
           },
           own = {
             filter_opt_list <- utils::modifyList(filter_opt_list, list(...))
           }
    )
    gc(verbose = FALSE)
    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
    counts <- do.call(sparse_filter, args = c(list("counts" = counts, "tx2gene" = tx2gene, "BPPARAM" = BPPARAM), filter_opt_list), quote = TRUE)
    BiocParallel::bpprogressbar(BPPARAM) <- FALSE
  }
  
  # Add a check here for non-finite values
  if (any(!is.finite(counts))) {
    stop("Non-finite values detected in the counts matrix. Consider applying a filtering strategy or cleaning your data before running the analysis.")
  }
  
  if (filter_only) {
    return(counts)
  }
  
  tx2gene <- tx2gene[match(rownames(counts), tx2gene$feature_id), ]
  
  if (methods::is(counts, "sparseMatrix") & force_dense) {
    counts <- tryCatch(
      {
        as.matrix(counts)
      },
      error = function(cond) {
        message(cond)
        stop(
          "Your sparse count matrix is probably too big and a non-sparse representation would need too much memory.",
          "\nTry subsetting or filtering the sparse matrix beforehand.\n\nOperation would require approximately ",
          format(structure(as.double(nrow(counts)) * as.double(ncol(counts)) * 8, class = "object_size"), units = "auto"), " of memory."
        )
      }
    )
  }
  
  drim <- sparseDRIMSeq::sparse_dmDSdata(tx2gene = tx2gene, counts = counts, samples = samp)
  exp_in_tx <- getFromNamespace("ratio_expression_in", "DTUrtle")(drim, "tx", BPPARAM = BPPARAM)
  exp_in_gn <- getFromNamespace("ratio_expression_in", "DTUrtle")(drim, "gene", BPPARAM = BPPARAM)
  
  #exp_in_tx <- as.numeric(exp_in_tx)
  #exp_in_gn <- as.numeric(exp_in_gn)
  
  # Diagnostic checks for non-finite values in ratio expressions
  # exp_in_tx_mat <- as.matrix(exp_in_tx)
  # if (any(!is.finite(exp_in_tx_mat))) {
  #   warning("Non-finite values detected in transcript-level expression ratios (exp_in_tx).")
  #   print(summary(exp_in_tx_mat))
  #   problematic_tx <- rownames(exp_in_tx_mat)[apply(exp_in_tx_mat, 1, function(x) any(!is.finite(x)))]
  #   message("Problematic transcript IDs: ", paste(problematic_tx, collapse = ", "))
  # }
  # 
  # exp_in_gn_mat <- as.matrix(exp_in_gn)
  # if (any(!is.finite(exp_in_gn_mat))) {
  #   warning("Non-finite values detected in gene-level expression ratios (exp_in_gn).")
  #   print(summary(exp_in_gn_mat))
  #   problematic_genes <- rownames(exp_in_gn_mat)[apply(exp_in_gn_mat, 1, function(x) any(!is.finite(x)))]
  #   message("Problematic gene IDs: ", paste(problematic_genes, collapse = ", "))
  # }
  # 
  
  
  if (carry_over_metadata & ncol(tx2gene) > 2) {
    tx2gene <- tx2gene[match(rownames(exp_in_tx), tx2gene$feature_id), ]
    tx2gene <- tx2gene[, !apply(tx2gene, 2, function(x) any(is.na(x)))]
    if (nrow(tx2gene) == nrow(exp_in_tx) & ncol(tx2gene) > 2) {
      exp_in_tx <- cbind(exp_in_tx, tx2gene[, -c(1, 2)], stringsAsFactors = FALSE)
      add_to_gene_columns <- check_unique_by_partition(tx2gene[, -c(1, 2)], drim@counts@partitioning)
      if (!is.null(add_to_gene_columns)) {
        tx2gene <- get_by_partition(df = tx2gene, columns = add_to_gene_columns, partitioning = drim@counts@partitioning, FUN = unique, BPPARAM = BPPARAM)
        exp_in_gn <- cbind(exp_in_gn, tx2gene[match(rownames(exp_in_gn), tx2gene[[1]]), -c(1)], stringsAsFactors = FALSE)
      }
    }
  }
  design_full <- stats::model.matrix(~condition, data = samp)
  
  filter_messages <- getFromNamespace("filter_messages", "DTUrtle")
  message("\nPerforming statistical tests...\n")
  drim_test <- filter_messages(sparseDRIMSeq::dmPrecision(drim, design = design_full, prec_subset = 1, BPPARAM = BPPARAM, add_uniform = add_pseudocount, verbose = 1))
  drim_test <- filter_messages(sparseDRIMSeq::dmFit(drim_test, design = design_full, BPPARAM = BPPARAM, add_uniform = add_pseudocount, verbose = 1))
  drim_test <- filter_messages(sparseDRIMSeq::dmTest(drim_test, coef = 2, BPPARAM = BPPARAM, verbose = 1))
  group <- factor(samp$condition, levels = cond_levels, ordered = TRUE)
  
  exp_in_gn <- rapply(exp_in_gn, as.character, classes = "factor", how = "replace")
  exp_in_tx <- rapply(exp_in_tx, as.character, classes = "factor", how = "replace")
  
  return_obj <- list(
    "meta_table_gene" = exp_in_gn, 
    "meta_table_tx" = exp_in_tx, 
    "meta_table_sample" = samp,
    "drim" = drim_test, 
    "design_full" = design_full, 
    "group" = group,
    "used_filtering_options" = list("DRIM" = filter_opt_list),
    "add_pseudocount" = add_pseudocount
  )
  class(return_obj) <- append("dturtle", class(return_obj))
  return(return_obj)
}
