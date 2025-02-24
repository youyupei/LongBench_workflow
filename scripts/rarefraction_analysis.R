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
library(purrr)
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


# function for splitting the transcript origin
is.sequins <- function(dge) {grepl("^R", rownames(dge))}
is.human <- function(dge) {grepl("^ENS", rownames(dge))}
is.sirv_or_ercc <- function(dge) {!is.human(dge) & !is.sequins(dge)} # this is the whole SIRV set4 (ERCC + SIRV E0 + Long SIRV)
is.sirv_E0 <- function(dge) {grepl("^SIRV\\d{3}$", rownames(dge))} 
is.sirv_all <- function(dge) {grepl("^SIRV", rownames(dge))} # including all SIRV E0 and Long SIRV
is.long_sirv <- function(dge) {grepl("^SIRV\\d{4}", rownames(dge))}
is.ercc <- function(dge) {!is.human(dge) & !is.sequins(dge) & !is.sirv_all(dge)}


# Set up environment
calcTxNum <- function(G.Tx.map){
  geneid = rep(NA, nrow(G.Tx.map))
  txcount = table(G.Tx.map$gene_id)
  txnum = as.numeric(txcount[match(G.Tx.map$gene_id, names(txcount))])
}
# SIRV

params <- list(
    random_seed = 2024,
    cache_dir = NULL,
    sirv_ercc_gtf = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf",
    sequins_gtf = "/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_annotation_2.4.gtf",
    sequins_tsv = "/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_isoforms_2.4.tsv",
    human_gtf = "/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf",
    sirv_csv = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_set4_concentration.csv",
    ont_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/ont_bulk",
    pb_bulk_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/pb_bulk",
    ont_drnd_oarfish_dir = "/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/dRNA_bulk",
    ill_bulk_salmon_dir = "/vast/projects/LongBench/analysis/sr_bulk/result/salmon/salmon_quant",
    bulk_meta = "/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt"
)
bulk.meta <- read.csv(params$bulk_meta)
rownames(bulk.meta) <- bulk.meta$sample
SIRV_tx <- GenomicFeatures::makeTxDbFromGFF(params$sirv_ercc_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
SIRV.G.Tx.map <- tibble(
  gene_id = SIRV_tx$gene_id %>% unlist,
  tx_name = SIRV_tx$tx_name
)
rm(SIRV_tx)
SIRV.G.Tx.map$tx_count <- SIRV.G.Tx.map %>% calcTxNum

# Sequins
sequins_tx  <- GenomicFeatures::makeTxDbFromGFF(params$sequins_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
sequins.G.Tx.map <- tibble(
  gene_id = sequins_tx$gene_id %>% unlist,
  tx_name = sequins_tx$tx_name
)
rm(sequins_tx)
sequins.G.Tx.map$tx_count <- sequins.G.Tx.map %>% calcTxNum

# Human
human_tx  <- GenomicFeatures::makeTxDbFromGFF(params$human_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
human.G.Tx.map <- tibble(
  gene_id = human_tx$gene_id %>% unlist,
  tx_name = human_tx$tx_name
)
rm(human_tx)
human.G.Tx.map$tx_count <- human.G.Tx.map %>% calcTxNum

# get the tx2gene mapping
combined_tx2gene <- rbind(human.G.Tx.map, sequins.G.Tx.map, SIRV.G.Tx.map)[, c("tx_name", "gene_id")]


# single process

single_process_oarfish <- function(quant_dir, method = "oarfish", length_interval=c(0,500,1000,2000, 3000, Inf)) {
    if (method == "oarfish") {
        tx.dge <- get_dge_from_oarfish(quant_dir, ".*/")
        gene.dge <- get_dge_from_txi(file.path(quant_dir %>% list.dirs(full.names = TRUE, recursive = FALSE), ".quant"), ".*/", "oarfish", dropInfReps=TRUE)
    } else if (method == "salmon") {
        tx.dge <- get_dge_from_salmon(quant_dir, ".*/")
        gene.dge <- get_dge_from_txi(file.path(quant_dir %>% list.dirs(full.names = TRUE, recursive = FALSE), "quant.sf"), ".*/", "salmon", dropInfReps=TRUE)
    } else {
        stop("method should be either 'oarfish' or 'salmon'")
    }

    # filter the dge for human genes only
    tx.dge <- tx.dge[is.human(tx.dge), ]
    gene.dge <- gene.dge[is.human(gene.dge), ]
    

    # Create labels dynamically
    labels <- paste(head(length_interval, -1), tail(length_interval, -1), sep = "-")
    tx.length <- tx.dge$genes$Length %>% cut(length_interval, labels =labels)
    gene.length <- gene.dge$genes$Length %>% cut(length_interval, labels =labels)

    rst <- data.frame(
      file = quant_dir,
      total.tx.ident = colSums(tx.dge$counts >= 5) %>% mean,
      total.gene.ident = colSums(gene.dge$counts >= 10) %>% mean
    )
    # add count by length bins
    for (x in labels) {
      rst[glue("{x}.tx.ident")] <- colSums(tx.dge[tx.length == x, ]$counts >= 5) %>% mean
      rst[glue("{x}.gene.ident")] <- colSums(gene.dge[gene.length == x, ]$counts >= 10) %>% mean
    }
    return(rst)
}

# read and put together the results
protocols <- c("ont_bulk", "pb_bulk", "dRNA_bulk")
sample_rate <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
combination <- expand.grid(protocols = protocols, sample_rate = sample_rate)
lr_quant_dirs <- glue::glue(
    "/vast/projects/LongBench/analysis/main_workflow/result/rarefraction_analysis/oarfish/{combination$protocols}/{combination$sample_rate}"
)

sr_quant_dirs <- glue::glue(
    "/vast/projects/LongBench/analysis/main_workflow/result/rarefraction_analysis/salmon/{sample_rate}"
)

lr_rst <- map(lr_quant_dirs, single_process_oarfish, method = "oarfish") %>%
    do.call(rbind, .)

sr_rst <- map(sr_quant_dirs, single_process_oarfish, method='salmon') %>%
    do.call(rbind, .)

rst <- rbind(lr_rst, sr_rst)
rst <- rst %>%
    mutate(
        protocols = map_chr(file, ~ strsplit(.x, "/")[[1]] %>%
            tail(2) %>%
            .[1]),
        sample_rate = map_chr(file, ~ strsplit(.x, "/")[[1]] %>%
            tail(1)) %>% as.numeric()
    ) %>%
    mutate(
        protocols = case_when(
            protocols == "salmon" ~ "Illumina",
            TRUE ~ protocols
        )
    ) %>% select(-file)

# save the results as a csv file
write.csv(rst, "/vast/projects/LongBench/analysis/main_workflow/result/rarefraction_analysis/rarefraction_tx_g_detection.csv", row.names = FALSE)


# plotting
# get read length table
read_number_table <- '/vast/projects/LongBench/analysis/figures/qc/read_number_table.txt'
read_number_table <- read.table(read_number_table, header = TRUE, sep = "\t")
# remove rows with "_sn" or "_sc" as suffix at the sample columns
read_number_table <- read_number_table %>%
                        filter(!grepl("_s(c|n)", sample)) %>%
                        rename(protocols = datatype) %>%
                        group_by(protocols) %>%
                        summarise(mean_read_count = mean(read_count))


# read csv
data <- read.csv(
  "/vast/projects/LongBench/analysis/main_workflow/result/rarefraction_analysis/rarefraction_tx_g_detection.csv",
  check.names = FALSE
)
data <- data %>%
            mutate(protocols = case_when(
                protocols == "Illumina" ~ "Illumina",
                protocols == "ont_bulk" ~ "ONT cDNA",
                protocols == "dRNA_bulk" ~ "ONT dRNA",
                protocols == "pb_bulk" ~ "PacBio"
            )) %>%
            left_join(read_number_table, by = "protocols") %>%
            mutate("Average read count" = sample_rate * mean_read_count)
# set level order of protocols
data$protocols <- factor(data$protocols, levels = c("Illumina","ONT cDNA", "ONT dRNA", "PacBio"))


# plot the rarefraction curve
p1 <- ggplot(data, aes(x = `Average read count`, y = total.tx.ident, color = protocols)) +
  geom_point() +
  geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Average read count",
    y = "Transcripts Detected",
    title = "Rarefraction Curve (Transcript)"
  )

# plot the rarefraction curve
p2 <- ggplot(data, aes(x = `Average read count`, y = total.gene.ident, color = protocols)) +
  geom_point() +
  geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Average read count",
    y = "Genes Detected",
    title = "Rarefaction Curve (Gene)"
  )

p1 / p2

# split based on length

lapply(labels, function(l) {
  ggplot(data, aes(x = `Average read count`, y = .data[[paste0(l, ".tx.ident")]], color = protocols)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "Average read count",
      y = "Transcripts Detected",
      title = "Rarefaction Curve (Transcript length {l})" %>% glue()
    )
}) %>% patchwork::wrap_plots(ncol = 3, guides = "collect") & theme(legend.position = "bottom")


lapply(labels, function(l) {
  ggplot(data, aes(x = `Average read count`, y = .data[[paste0(l, ".gene.ident")]], color = protocols)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "Average read count",
      y = "Genes Detected",
      title = "Rarefaction Curve (Gene length  {l})" %>% glue()
    )
}) %>% patchwork::wrap_plots(ncol = 3, guides = "collect") & theme(legend.position = "bottom")


lapply(labels, function(l) {
  ggplot(data, aes(x = `Average read count`, y = .data[[paste0(l, ".tx.ident")]] /.data[[paste0(l, ".gene.ident")]], color = protocols)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "Average read count",
      y = "Transcript per genes detected",
      title = "Rarefaction Curve (Gene length {l})" %>% glue()
    )
}) %>% patchwork::wrap_plots(ncol = 3, guides = "collect") & theme(legend.position = "bottom")