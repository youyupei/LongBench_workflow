library(cluster, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)

get_silhouette_score <- function(seurat_obj, cluster_col, dims, reduction) {
  # remove cells with no cluster assignment (remove NA)
  seurat_obj <- seurat_obj[,!is.na(seurat_obj[[cluster_col]][,1])]
  dist.matrix <- dist(x = Embeddings(object = seurat_obj[[reduction]])[, dims])
  clusters <- seurat_obj[[cluster_col]][,1]
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  return(sil[, 3])
}




# Load the Seurat object
cell_lines <- c()
silhouette_scores <- c()
sample_names <- c()
so.pb.sc <- readRDS("/vast/projects/LongBench/analysis/lr_sc_sn/result/reports/RDS/pb_sc_annotated.rds")

so.pb.sc.sil <- get_silhouette_score(
  so.pb.sc,
  cluster_col = "cell_lines",
  dims = c(1:20),
  reduction = "pca.Gene_quant"
)
silhouette_scores <- c(silhouette_scores, so.pb.sc.sil)
cell_lines <- c(cell_lines, so.pb.sc$cell_lines)
sample_names <- c(sample_names, rep("PB SC", length(so.pb.sc.sil)))
rm(so.pb.sc)

so.pb.sn <- readRDS("/vast/projects/LongBench/analysis/lr_sc_sn/result/reports/RDS/pb_sn_annotated.rds")

so.pb.sn.sil <- get_silhouette_score(
  so.pb.sn,
  cluster_col = "cell_lines",
  dims = c(1:20),
  reduction = "pca.Gene_quant"
)
silhouette_scores <- c(silhouette_scores, so.pb.sn.sil)
cell_lines <- c(cell_lines, so.pb.sn$cell_lines)
sample_names <- c(sample_names, rep("PB SN", length(so.pb.sn.sil)))
rm(so.pb.sn)

so.ont.sc <- readRDS("/vast/projects/LongBench/analysis/lr_sc_sn/result/reports/RDS/ont_sc_clean_annotated.rds")
so.ont.sc.sil <- get_silhouette_score(
  so.ont.sc,
  cluster_col = "cell_lines",
  dims = c(1:20),
  reduction = "pca.Gene_quant"
)
silhouette_scores <- c(silhouette_scores, so.ont.sc.sil)
cell_lines <- c(cell_lines, so.ont.sc$cell_lines)
sample_names <- c(sample_names, rep("ONT SC", length(so.ont.sc.sil)))
rm(so.ont.sc)

so.ont.sn <- readRDS("/vast/projects/LongBench/analysis/lr_sc_sn/result/reports/RDS/ont_sn_clean_annotated.rds")
so.ont.sn.sil <- get_silhouette_score(
  so.ont.sn,
  cluster_col = "cell_lines",
  dims = c(1:20),
  reduction = "pca.Gene_quant"
)
silhouette_scores <- c(silhouette_scores, so.ont.sn.sil)
cell_lines <- c(cell_lines, so.ont.sn$cell_lines)
sample_names <- c(sample_names, rep("ONT SN", length(so.ont.sn.sil)))
rm(so.ont.sn)

# plot silhouette score per assay per cell line
sil_df <- data.frame(
  assay = sample_names,
  silhouette_score = silhouette_scores,
  cell_line = cell_lines
)

# get the mean silhouette score per cell line
plot_df <- sil_df %>% group_by(assay, cell_line) %>% summarise(median_silhouette_score = median(silhouette_score))

# box plot
p1 <- ggplot(sil_df, aes(x = assay, y = silhouette_score, fill = cell_line)) +
geom_boxplot(outlier.shape = NA)  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-0.1, 0.75) +
  labs(title = "Silhouette score", x = "Cell line", y = "Silhouette score")



# boxplot of median silhouette score per cell line
# boxplot of median silhouette score per cell line with dots
p2 <- ggplot(plot_df, aes(x = assay, y = median_silhouette_score)) +
    # boxplot of median silhouette score per cell line without showing outliers
    geom_boxplot(outlier.shape = NA) +
    # add geom_jitter to show the data points and color by cell line
    geom_jitter(aes(color = cell_line), width = 0.1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Median Silhouette score of each cell line", x = "Cell line", y = "Median Silhouette score")

# save the plots p1 | p2 as pdf
pdf("/vast/projects/LongBench/analysis/lr_sc_sn/result/plots/silhouette_score.pdf")
print(p1)
print(p2)
dev.off()
