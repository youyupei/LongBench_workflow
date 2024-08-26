# function for preprocessing the FALMES output
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)
library(knitr)
library(DT)
library(glue)
library(tibble)


# function loads the gene quantification matrix from flames output and returns a Seurat object with gene names to gene ids map df stored in so@misc$gene_info_map
#' Load gene quantification data and create a Seurat object
#'
#' This function reads a gene quantification file, converts gene IDs to gene names,
#' filters out invalid gene names, and creates a Seurat object.
#'
#' @param tenX_out The file path to the 10X output file.
#' @param project The name of the project. Default is "singlecell".
#' @param min.cells The minimum number of cells a gene must be expressed in to be included. Default is 3.
#' @param min.features The minimum number of features (genes) a cell must express to be included. Default is 0.
#'
#' @return A Seurat object containing the gene expression data.
#'
#' @examples
#' # Load gene quantification data and create a Seurat object
#' seurat_obj <- load_gene_quant("path/to/gene_quantification.csv", project = "my_project", min.cells = 5, min.features = 10)
#'
load_longbench_result <- function(tenX_out, project = "ill_sc", min.cells = 3, min.features = 0, assay = "Gene_quant") {
    so <- Read10X(data.dir=tenX_out)
    # set min.cells and min.features
    so <- CreateSeuratObject(counts = so, project = project, min.cells = min.cells, min.features = min.features, assay=assay)
  
    # Get mito and ribo gene percentage
    so[["percent.mt"]] <- PercentageFeatureSet(so, features = grep("^MT", rownames(so), value = T)) %>% round(2)
    so[["percent.ribo"]] <- PercentageFeatureSet(so, features = grep("^RP[SL]", rownames(so), value = T)) %>% round(2)
    so <- add_qc_state(so, "raw_gene_quant")
    return(so)
}
  

# REMOVE DOUBLETS using scDblFinder
remove.doublet <- function(so, thread=8, RNGseed=2024, assay = NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(so)
  }
  bp <- BiocParallel::MulticoreParam(thread, RNGseed=RNGseed)
  dbl <- scDblFinder::scDblFinder(Seurat::GetAssayData(so, layer = "counts", assay = assay),
                   clusters = FALSE,
                   processing="default", #normal scater-based normalization and PCA
                  # processing="normFeatures", #(uses normalized features, without PCA
                  returnType = c("sce", "table", "full", "counts"),
                   verbose = TRUE,
                   BPPARAM = bp)
  so <- subset(so, cells=colnames(dbl)[dbl$scDblFinder.class=='singlet'])
  so@misc$doublet.pct <- mean(dbl$scDblFinder.class=="doublet")
  so <- add_qc_state(so, "doublet_removal")
  return(so)
}


add_qc_state <- function(so, filter, assay=NULL) {
  if (is.null(assay)) {
    assay <- DefaultAssay(so)
  }
  if (!is.null(so@misc$qc_status)) {
    so@misc$qc_status <- add_row(
      so@misc$qc_status,
      filter = filter,
      n_Cells = so %>% ncol(),
      median.read.count = so@meta.data[["nCount_{assay}" %>% glue]] %>% median
    )
  } else {
    so@misc$qc_status <- tibble(
      filter = filter,
      n_Cells = so %>% ncol(),
      median.read.count = so@meta.data[["nCount_{assay}" %>% glue]] %>% median
    )
  }
  return(so)
}

# Normalisation, Dimensionality reduction and Clustering
get_clusters <- function(so, assay=NULL, ndim=30,res=0.5,seed = 2024) {
  if (is.null(assay)) {
    assay <- DefaultAssay(so)
  }
  so <- SCTransform(so, assay = assay) %>% 
      RunPCA(.,  reduction.name=glue("pca.{assay}")) %>% 
      RunUMAP(reduction=glue('pca.{assay}') ,dims = 1:ndim, reduction.name =glue( 'umap.{assay}'), reduction.key = glue('{assay}_UMAP_'))
      so <- FindNeighbors(so, dims = 1:ndim, reduction=glue("pca.{assay}"))
      so <- FindClusters(so, resolution = 0.5, seed = seed, cluster.name=glue("seurat_clusters_{assay}_res{res}"))
  return(so)
}


remove_outliers <- function(so, assay=NULL){
  if (is.null(assay)) {
    assay <- DefaultAssay(so)
  }

  p.before <- VlnPlot(SetIdent(so, value = 'orig.ident'), features = c("nFeature_{assay}" %>% glue, "nCount_{assay}"%>% glue, "percent.mt",'percent.ribo'), ncol = 4, alpha=0.2)
  is.outlier <- scuttle::isOutlier(so@meta.data[["nFeature_{assay}" %>% glue]] %>% log, nmad=3) |
                  scuttle::isOutlier(so@meta.data[["nCount_{assay}" %>% glue]]  %>% log, nmad=3) |
                  scuttle::isOutlier(so$percent.mt, nmad=3)
  so <- subset(x=so, cells=names(is.outlier)[!is.outlier])
  so <- add_qc_state(so, "outlier_removal")
  p.after <- VlnPlot(SetIdent(so, value = 'orig.ident'), features = c("nFeature_{assay}" %>% glue, "nCount_{assay}"%>% glue, "percent.mt",'percent.ribo'), ncol = 4, alpha=0.2)
  print(p.before / p.after)
  return(so)
}
