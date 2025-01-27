library(limma)
library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(glue)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(knitr)
library(kableExtra)
library(rtracklayer)
library(scales)
library(UpSetR)
library(ggtree)
library(aplot)
library(tximport)
library(ComplexUpset)
source("Rfunctions.R")

# setup some conflict functions
rename <- dplyr::rename
select <- dplyr::select

# Set the color palette
color_palette <- c(
    PacBio = "#df1995",
    ONT = "#00789b",
    ONT_1 = "#04476c", # ONT cDNA
    ONT_2 = "#24cdcd", # ONT dRNA
    Illumina = "#e88b20"
)



# Set up environment
calcTxNum <- function(G.Tx.map){
  geneid = rep(NA, nrow(G.Tx.map))
  txcount = table(G.Tx.map$gene_id)
  txnum = as.numeric(txcount[match(G.Tx.map$gene_id, names(txcount))])
}
# SIR
params <- list(
    random_seed = 2024,
    cache_dir = NULL,
    gtf_C = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_20210507.gtf",
    gtf_I = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_I_20210507.gtf",
    gtf_O = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_O_20210507.gtf",
    sirv_csv = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_set4_concentration.csv",
    ont_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/oarfish/ont_bulk",
    pb_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/oarfish/pb_bulk",
    ont_drnd_oarfish_dir = "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/oarfish/dRNA_bulk",
    ill_bulk_salmon_dir = "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/salmon/salmon_quant",
    bulk_meta = "/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt"
)

bulk.meta <- read.csv(params$bulk_meta)
rownames(bulk.meta) <- bulk.meta$sample
SIRV_tx.tmp.load <- GenomicFeatures::makeTxDbFromGFF(params$gtf_C, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
SIRV.G.Tx.map.C <- tibble(
  gene_id = SIRV_tx.tmp.load$gene_id %>% unlist,
  tx_name = SIRV_tx.tmp.load$tx_name
)
rm(SIRV_tx.tmp.load)
SIRV.G.Tx.map.C$tx_count <- SIRV.G.Tx.map.C %>% calcTxNum

SIRV_tx.tmp.load <- GenomicFeatures::makeTxDbFromGFF(params$gtf_O, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
SIRV.G.Tx.map.O <- tibble(
  gene_id = SIRV_tx.tmp.load$gene_id %>% unlist,
  tx_name = SIRV_tx.tmp.load$tx_name
)
rm(SIRV_tx.tmp.load)
SIRV.G.Tx.map.O$tx_count <- SIRV.G.Tx.map.O %>% calcTxNum


SIRV_tx.tmp.load <- GenomicFeatures::makeTxDbFromGFF(params$gtf_I, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
SIRV.G.Tx.map.I <- tibble(
  gene_id = SIRV_tx.tmp.load$gene_id %>% unlist,
  tx_name = SIRV_tx.tmp.load$tx_name
)
rm(SIRV_tx.tmp.load)
SIRV.G.Tx.map.I$tx_count <- SIRV.G.Tx.map.I %>% calcTxNum

# Join three SIRV Tx maps based on tx_name
SIRV.G.Tx.map <- full_join(SIRV.G.Tx.map.C, SIRV.G.Tx.map.O, by = join_by(gene_id == gene_id, tx_name == tx_name), suffix = c(".C", ".O")) %>%
    full_join(SIRV.G.Tx.map.I, by = join_by(gene_id == gene_id, tx_name == tx_name)) %>%
    rename(tx_count.I = tx_count) %>%
    mutate(
        in.C = !is.na(tx_count.C),
        in.O = !is.na(tx_count.O),
        in.I = !is.na(tx_count.I)
    )

# Get a Venen diagram of the columns in.C, in.O, in.I
# Create the upset plot
# ComplexUpset::upset(
#   SIRV.G.Tx.map,
#   intersect = c("in.C", "in.O", "in.I"), # Columns to intersect
#   name = "Overlap",
#   width_ratio = 0.2
# )

# Create a venn diagram
# Create lists for the Venn diagram
list_C <- SIRV.G.Tx.map$tx_name[SIRV.G.Tx.map$in.C]
list_O <- SIRV.G.Tx.map$tx_name[SIRV.G.Tx.map$in.O]
list_I <- SIRV.G.Tx.map$tx_name[SIRV.G.Tx.map$in.I]

# Combine into a list
venn_list <- list("Complete annotation" = list_C, "Redundent annotation" = list_O, "Incomplete annotation" = list_I)

# Generate the Venn diagram
venn.plot <- VennDiagram::venn.diagram(
    x = venn_list,
    filename = NULL, # Set to NULL to plot directly in RStudio
    col = "#00000000",
    fill = c("#4ead05", "#ff7300", "#0000ff"),
    alpha = 0.5,
    cex = 1.5,
    cat.pos = 0,
    cat.dist = -0.015
)
# save the plot
# save the plot
svg("/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/SIRV_venn.svg")
grid::grid.draw(venn.plot)
dev.off()

# single process
combined_tx2gene <- SIRV.G.Tx.map[, c("tx_name", "gene_id")]
single_process_oarfish <- function(quant_dir, method = "oarfish") {
    if (method == "oarfish") {
        tx.dge <- get_dge_from_oarfish(quant_dir, ".*/")
        gene.dge <- get_dge_from_txi(file.path(quant_dir %>% list.dirs(full.names = TRUE, recursive = FALSE), ".quant"), ".*/", "oarfish", dropInfReps=TRUE)
    } else if (method == "salmon") {
        tx.dge <- get_dge_from_salmon(quant_dir, ".*/")
        gene.dge <- get_dge_from_txi(file.path(quant_dir %>% list.dirs(full.names = TRUE, recursive = FALSE), "quant.sf"), ".*/", "salmon", dropInfReps=TRUE)
    } else {
        stop("method should be either 'oarfish' or 'salmon'")
    }

    rst <- list(
        file = quant_dir,
        tx.dge = tx.dge,
        gene.dge = gene.dge
    )
}

# read and put together the results
protocols <- c("ont_bulk", "pb_bulk", "dRNA_bulk")
annotation <- c("C", "O", "I")
combination <- expand.grid(protocols = protocols, annotation = annotation)
lr_quant_dirs <- glue::glue(
    "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/oarfish/{combination$protocols}/{combination$annotation}"
)

sr_quant_dirs <- glue::glue(
    "/vast/projects/LongBench/analysis/main_workflow/result/annotation_redundency_analysis/salmon/salmon_quant/{annotation}"
)

lr_rst <- map(lr_quant_dirs, single_process_oarfish, method = "oarfish")
names(lr_rst) <- lr_quant_dirs %>%
      map_chr(~ paste0(tail(strsplit(.x, "/")[[1]], 2), collapse = "."))

lr_rst <- lapply(names(lr_rst), function(x) {
    lr_rst[[x]]$gene.dge$samples$protocols <- strsplit(x, "\\.")[[1]][1]
    lr_rst[[x]]$gene.dge$samples$annotation <- strsplit(x, "\\.")[[1]][2]
    lr_rst[[x]]$gene.dge$genes$annotation <- strsplit(x, "\\.")[[1]][2]
    lr_rst[[x]]$tx.dge$samples$genes <- strsplit(x, "\\.")[[1]][1]
    lr_rst[[x]]$tx.dge$samples$protocols <- strsplit(x, "\\.")[[1]][1]
    lr_rst[[x]]$tx.dge$samples$annotation <- strsplit(x, "\\.")[[1]][2]
    lr_rst[[x]]$tx.dge$genes$protocols <- strsplit(x, "\\.")[[1]][1]
    lr_rst[[x]]$tx.dge$genes$annotation <- strsplit(x, "\\.")[[1]][2]
    lr_rst[[x]]$tx.dge$ltpm <- lr_rst[[x]]$tx.dge$counts %>% cpm(prior.count=1, log=TRUE)
    lr_rst[[x]]
})
names(lr_rst) <- lr_quant_dirs %>%
      map_chr(~ paste0(tail(strsplit(.x, "/")[[1]], 2), collapse = "."))



sr_rst <- map(sr_quant_dirs, single_process_oarfish, method = "salmon")
names(sr_rst) <- paste0("ill_bulk.", annotation)
sr_rst <- lapply(names(sr_rst), function(x) {
    sr_rst[[x]]$gene.dge$samples$protocols <- strsplit(x, "\\.")[[1]][1]
    sr_rst[[x]]$gene.dge$samples$annotation <- strsplit(x, "\\.")[[1]][2]
    sr_rst[[x]]$gene.dge$genes$annotation <- strsplit(x, "\\.")[[1]][2]
    sr_rst[[x]]$tx.dge$samples$genes <- strsplit(x, "\\.")[[1]][1]
    sr_rst[[x]]$tx.dge$samples$protocols <- strsplit(x, "\\.")[[1]][1]
    sr_rst[[x]]$tx.dge$samples$annotation <- strsplit(x, "\\.")[[1]][2]
    sr_rst[[x]]$tx.dge$genes$protocols <- strsplit(x, "\\.")[[1]][1]
    sr_rst[[x]]$tx.dge$genes$annotation <- strsplit(x, "\\.")[[1]][2]
    # scale the counts by the effective length
    length_norm_count  <- sr_rst[[x]]$tx.dge$counts / sr_rst[[x]]$tx.dge$genes$EffectiveLength
    sr_rst[[x]]$tx.dge$ltpm <- length_norm_count %>% cpm(prior.count = 1, log = TRUE)
    sr_rst[[x]]
})
names(sr_rst) <- paste0("ill_bulk.", annotation)
rst <- c(lr_rst, sr_rst)

# overdisp_d
# tx_overdisp_d <- map(rst, ~ .x$tx.dge$genes)

# TODO: overdispersion vs annotation abundance


# How the quantification affected by the annotation completeness and redundancy
# ## GENE level
# rst$ill_bulk.I$gene.dge$counts
# 
# Tx_count_d <- SIRV.G.Tx.map.C %>%
#     group_by(gene_id) %>%
#     summarise(tx_count = n())
#     
# 
# merged_gene_count <- map(names(rst), ~ (rst[[.x]]$gene.dge$counts / Tx_count_d$tx_count) %>%
#     data.frame() %>%
#     rownames_to_column("Gene_id") %>%
#     pivot_longer(-Gene_id, names_to = "Sample", values_to = "Count") %>%
#     mutate(protocol = .x) 
# ) %>% do.call(rbind, .)
# 
# ## Normalised gene count per sample
# merged_gene_count %>% group_by(protocol, Sample) %>% summarise(mean = mean(Count)) %>% kable()

## TX level
merged_txi_count <- map(names(rst), ~ {
    rst[[.x]]$tx.dge$counts %>%
        data.frame() %>%
        rownames_to_column("Gene_id") %>%
        pivot_longer(-Gene_id, names_to = "Sample", values_to = "Count") %>%
        # group_by(Sample) %>%
        # summarise(mad = mad(Count)) %>%
        mutate(protocol = strsplit(.x, "\\.")[[1]][1], annotation = strsplit(.x, "\\.")[[1]][2]) %>%
        left_join(
            rst[[.x]]$tx.dge$ltpm %>%
                data.frame() %>%
                rownames_to_column("Gene_id") %>%
                pivot_longer(-Gene_id, names_to = "Sample", values_to = "ltpm")
        )
}) %>% do.call(rbind, .)

merged_txi_count <- merged_txi_count %>% mutate(tx_anno = case_when(
    merged_txi_count$Gene_id %in% list_I ~ "Incomplete",
    merged_txi_count$Gene_id %in% list_C ~ "Complete",
    merged_txi_count$Gene_id %in% list_O ~ "Redundent"
))

## order the level
merged_txi_count$annotation <- factor(merged_txi_count$annotation, levels = c("I", "C", "O"))
merged_txi_count$tx_anno <- factor(merged_txi_count$tx_anno, levels = c( "Redundent", "Complete","Incomplete"))
merged_txi_count$protocol <- factor(merged_txi_count$protocol, levels = c("ill_bulk", "ont_bulk", "pb_bulk", "dRNA_bulk"))

# TOTAL COUNT
cat_total <- merged_txi_count %>%
    group_by(protocol, annotation, tx_anno) %>%
    summarise(total = sum(Count)) %>%
    group_by(protocol) %>%
    mutate(norm_total = total / sum(total[annotation == "C"]))

cat_total %>% ggplot(aes(x = annotation, y = norm_total, fill = tx_anno)) +
    geom_bar(stat = "identity") +
    facet_grid(~protocol) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Total counts per annotation completeness", x = "Annotation completeness", y = "Normalised Total counts")


p1 <- merged_txi_count %>%
    filter( protocol == "ill_bulk") %>%
    ggplot(aes(x = annotation, y = ltpm, fill = tx_anno)) +
    geom_boxplot() +
    stat_summary(
    fun.data = function(y) {
      data.frame(
        y = quantile(y, 0.25) - 0.8,  # Slightly above the upper whisker
        label = paste0("MAD: ", round(mad(y), 2), "\n", "SD: ", round(sd(y), 2))
      )
    },
    geom = "text",
    position = position_dodge(width = 0.75),  # Properly align annotations
    size = 3
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Illumina", x = "Annotation completeness", y = "Log2 (TPM+1)") 

p2 <- merged_txi_count %>%
    filter( protocol == "ont_bulk") %>%
    ggplot(aes(x = annotation, y = ltpm, fill = tx_anno)) +
    geom_boxplot() +
    stat_summary(
    fun.data = function(y) {
      data.frame(
        y = quantile(y, 0.25) - 0.8,  # Slightly above the upper whisker
        label = paste0("MAD: ", round(mad(y), 2), "\n", "SD: ", round(sd(y), 2))
      )
    },
    geom = "text",
    position = position_dodge(width = 0.75),  # Properly align annotations
    size = 3
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ONT cDNA", x = "Annotation completeness", y = "Log2 (TPM+1)")

p3 <- merged_txi_count %>%
    filter(Sample=="H146", protocol == "dRNA_bulk") %>%
    ggplot(aes(x = annotation, y = ltpm, fill = tx_anno)) +
    geom_boxplot() +
    stat_summary(
    fun.data = function(y) {
      data.frame(
        y = quantile(y, 0.25) - 0.8,  # Slightly above the upper whisker
        label = paste0("MAD: ", round(mad(y), 2), "\n", "SD: ", round(sd(y), 2))
      )
    },
    geom = "text",
    position = position_dodge(width = 0.75),  # Properly align annotations
    size = 3
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ONT dRNA", x = "Annotation completeness", y = "Log2 (TPM+1)")

p4 <- merged_txi_count %>%
    filter(Sample=="H146", protocol == "pb_bulk") %>%
    ggplot(aes(x = annotation, y = ltpm, fill = tx_anno)) +
    geom_boxplot() +
    stat_summary(
    fun.data = function(y) {
      data.frame(
        y = quantile(y, 0.25) - 0.8,  # Slightly above the upper whisker
        label = paste0("MAD: ", round(mad(y), 2), "\n", "SD: ", round(sd(y), 2))
      )
    },
    geom = "text",
    position = position_dodge(width = 0.75),  # Properly align annotations
    size = 3
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Pacbio", x = "Annotation completeness", y = "Log2 (TPM+1)")

(p1 | p2) / (p3 | p4)


