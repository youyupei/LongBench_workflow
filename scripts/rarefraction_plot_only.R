library(tidyverse)
library(ggplot2)
library(glue)
library(patchwork)

rename <- dplyr::rename
select <- dplyr::select

color_palette <- c(
    PacBio = "#df1995",
    ONT_2 = "#00789b",
    ONT_1 = "#04476c",
    Illumina = "#e88b20"
)
labels <- c("0-500", "500-1000", "1000-2000", "2000-3000", "3000-Inf")

fig_dir <- "/vast/projects/LongBench/analysis/figures/rarefraction_analysis"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_number_table <- read.table(
    "/vast/projects/LongBench/analysis/figures/qc/read_number_table.txt",
    header = TRUE, sep = "\t"
)
read_number_table <- read_number_table %>%
    filter(!grepl("_s(c|n)", sample)) %>%
    rename(protocols = datatype) %>%
    group_by(protocols) %>%
    summarise(
        mean_read_count = mean(read_count),
        mean_read_length = mean(gb * 1e9 / read_count)
    )

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
    mutate(
        "Average read count" = sample_rate * mean_read_count,
        "Average gigabases" = sample_rate * mean_read_count * mean_read_length / 1e9
    )
data$protocols <- factor(data$protocols, levels = c("Illumina", "ONT cDNA", "ONT dRNA", "PacBio"))

# helper to build a single rarefraction plot
make_plot <- function(xvar, yvar, ylabel, title) {
    ggplot(data, aes(x = .data[[xvar]], y = {{ yvar }}, color = protocols)) +
        geom_point() + geom_line() +
        scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        labs(x = xvar, y = ylabel, title = title)
}

# read-count axis: transcript / gene / tx-per-gene
p1 <- ggplot(data, aes(x = `Average read count`, y = total.tx.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average read count", y = "Transcripts Detected", title = "Rarefraction Curve (Transcript)")

p2 <- ggplot(data, aes(x = `Average read count`, y = total.gene.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average read count", y = "Genes Detected", title = "Rarefaction Curve (Gene)")

p3 <- ggplot(data, aes(x = `Average read count`, y = total.tx.ident / total.gene.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average read count", y = "Transcripts per Gene Detected", title = "Rarefaction Curve (Transcripts per Gene)")

ggsave(file.path(fig_dir, "rarefraction_reads_tx_gene_txpergene.pdf"), p1 / p2 / p3, width = 6, height = 12)

# gigabases axis: transcript / gene / tx-per-gene
p1_gb <- ggplot(data, aes(x = `Average gigabases`, y = total.tx.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average gigabases sequenced", y = "Transcripts Detected", title = "Rarefraction Curve (Transcript)")

p2_gb <- ggplot(data, aes(x = `Average gigabases`, y = total.gene.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average gigabases sequenced", y = "Genes Detected", title = "Rarefaction Curve (Gene)")

p3_gb <- ggplot(data, aes(x = `Average gigabases`, y = total.tx.ident / total.gene.ident, color = protocols)) +
    geom_point() + geom_line() +
    scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
    theme_minimal() + theme(legend.position = "bottom") +
    labs(x = "Average gigabases sequenced", y = "Transcripts per Gene Detected", title = "Rarefaction Curve (Transcripts per Gene)")

ggsave(file.path(fig_dir, "rarefraction_gb_tx_gene_txpergene.pdf"), p1_gb / p2_gb / p3_gb, width = 6, height = 12)

# length-stratified plots — theme applied per panel, no & operator
make_len_plots <- function(xvar, yvar_suffix, ylabel) {
    lapply(labels, function(l) {
        ggplot(data, aes(x = .data[[xvar]], y = .data[[paste0(l, yvar_suffix)]], color = protocols)) +
            geom_point() + geom_line() +
            scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
            theme_minimal() + theme(legend.position = "bottom") +
            labs(x = xvar, y = ylabel, title = glue("Rarefaction Curve (length {l})"))
    }) %>% wrap_plots(ncol = 3, guides = "collect")
}

ggsave(file.path(fig_dir, "rarefraction_reads_tx_by_length.pdf"),
       make_len_plots("Average read count", ".tx.ident", "Transcripts Detected"), width = 15, height = 10)
ggsave(file.path(fig_dir, "rarefraction_gb_tx_by_length.pdf"),
       make_len_plots("Average gigabases", ".tx.ident", "Transcripts Detected"), width = 15, height = 10)

ggsave(file.path(fig_dir, "rarefraction_reads_gene_by_length.pdf"),
       make_len_plots("Average read count", ".gene.ident", "Genes Detected"), width = 15, height = 10)
ggsave(file.path(fig_dir, "rarefraction_gb_gene_by_length.pdf"),
       make_len_plots("Average gigabases", ".gene.ident", "Genes Detected"), width = 15, height = 10)

make_len_ratio_plots <- function(xvar) {
    lapply(labels, function(l) {
        ggplot(data, aes(
            x = .data[[xvar]],
            y = .data[[paste0(l, ".tx.ident")]] / .data[[paste0(l, ".gene.ident")]],
            color = protocols
        )) +
            geom_point() + geom_line() +
            scale_color_manual(values = color_palette[c("Illumina", "ONT_1", "ONT_2", "PacBio")] %>% unname()) +
            theme_minimal() + theme(legend.position = "bottom") +
            labs(x = xvar, y = "Transcripts per gene detected", title = glue("Rarefaction Curve (length {l})"))
    }) %>% wrap_plots(ncol = 3, guides = "collect")
}

ggsave(file.path(fig_dir, "rarefraction_reads_txpergene_by_length.pdf"),
       make_len_ratio_plots("Average read count"), width = 15, height = 10)
ggsave(file.path(fig_dir, "rarefraction_gb_txpergene_by_length.pdf"),
       make_len_ratio_plots("Average gigabases"), width = 15, height = 10)

message("All figures saved to: ", fig_dir)
