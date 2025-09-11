library(tidyverse)
library(arrow)
library(patchwork)
library(Biostrings)
library(magrittr)
library(ComplexHeatmap)

library(conflicted)

conflicts_prefer(
    dplyr::select,
    dplyr::rename,
    dplyr::filter,
    base::intersect,
    base::setdiff
)

bed_files_df <- tibble(
    path = fs::dir_ls("data/modkit_pileup", glob = "*.bed"),
    sample_id = fs::path_file(fs::path_ext_remove(path))
)

sample_anno <- tibble::tribble(
    ~sample, ~index, ~sequins, ~group,
    "H146",1,"B","sclc_a",
    "H69",2,"B","sclc_a",
    "SHP77",5,"B","sclc_a",
    "H526",3,"A","sclc_p",
    "H211",4,"A","sclc_p",
    "H1975",6,"A","luad",
    "H2228",7,"A","luad",
    "HCC827",8,"A","luad",
)

mod_data_db <- open_dataset("data/modkit_pileup_parquet")

# filter and read the relevant columns from the parquet file
mod_data <- mod_data_db %>%
    filter(str_detect(chrom, "chr")) %>%
    filter(n_valid_cov > 20) %>%
    select(run, sample, chrom, end, strand, n_valid_cov, fraction_modified) %>%
    rename(
        pos = end,
        percent_modified = fraction_modified
    )

# calculate the number of samples with any modification
pos_mod_stats <- mod_data %>%
    summarise(
        n_with_mod = sum(percent_modified > 20),
        mean_cov = mean(n_valid_cov),
        .by = c(chrom, pos)
    )

# filter to positions that have some non-zero modification percents in 2 or more samples
pos_with_mods <- pos_mod_stats %>%
    filter(n_with_mod >= 2) %>%
    select(chrom, pos)

mod_data_pos_filtered <- left_join(pos_with_mods, mod_data, by = join_by(chrom, pos)) %>%
    collect() %>%
    mutate(
        sample_id = paste(run, sample, sep = "_")
    )

sample_mod_data <- mod_data_pos_filtered %>%
    filter(sample == "SHP77") %>%
    select(chrom, pos, strand, run, percent_modified) %>%
    pivot_wider(names_from = run, values_from = percent_modified) %>%
    relocate(run1, .before = run2)

cor(sample_mod_data$run1, sample_mod_data$run2, use = "pairwise.complete.obs")

sample_corr <- function(data, sample1, sample2) {
    filtered_data <- data %>%
        filter(sample_id == sample1 | sample_id == sample2) %>%
        select(chrom, pos, strand, sample_id, percent_modified) %>%
        pivot_wider(names_from = sample_id, values_from = percent_modified) %>%
        relocate(!!sample1, .before = !!sample2)

    cor(filtered_data[[sample1]], filtered_data[[sample2]], use = "pairwise.complete.obs")
}

unique_samples <- unique(mod_data_pos_filtered$sample_id)
unique_samples <- unique_samples[
    order(
        str_split(unique_samples, "_", simplify = TRUE)[, 2],
        str_split(unique_samples, "_", simplify = TRUE)[, 1]
    )
]

corr_matrix <- matrix(0, nrow = length(unique_samples), ncol = length(unique_samples)) %>%
    set_rownames(unique_samples) %>%
    set_colnames(unique_samples)

for (i in seq_along(unique_samples)) {
    for (j in seq_along(unique_samples)) {
        if (i != j) {
            corr_matrix[i, j] <- sample_corr(mod_data_pos_filtered, unique_samples[i], unique_samples[j])
        }
    }
}

for (i in seq_along(unique_samples)) {
    corr_matrix[i, i] <- 1
}

# > sample_df
# # A tibble: 16 × 6
#    sample_id   sample run   index sequins group
#    <chr>       <chr>  <chr> <dbl> <chr>   <chr>
#  1 run1_H146   H146   run1      1 B       sclc_a
#  2 run2_H146   H146   run2      1 B       sclc_a
#  3 run1_H1975  H1975  run1      6 A       luad
#  4 run2_H1975  H1975  run2      6 A       luad
#  5 run1_H211   H211   run1      4 A       sclc_p
#  6 run2_H211   H211   run2      4 A       sclc_p
#  7 run1_H2228  H2228  run1      7 A       luad
#  8 run2_H2228  H2228  run2      7 A       luad
#  9 run1_H526   H526   run1      3 A       sclc_p
# 10 run2_H526   H526   run2      3 A       sclc_p
# 11 run1_H69    H69    run1      2 B       sclc_a
# 12 run2_H69    H69    run2      2 B       sclc_a
# 13 run1_HCC827 HCC827 run1      8 A       luad
# 14 run2_HCC827 HCC827 run2      8 A       luad
# 15 run1_SHP77  SHP77  run1      5 B       sclc_a
# 16 run2_SHP77  SHP77  run2      5 B       sclc_a

sample_df <- tibble(
    sample_id = unique_samples,
    sample = str_split(unique_samples, "_", simplify = TRUE)[, 2],
    run = str_split(unique_samples, "_", simplify = TRUE)[, 1]
) %>%
left_join(sample_anno, by = join_by(sample))

set1 <- c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF", "#984EA3FF", "#FF7F00FF", "#FFFF33FF", "#A65628FF", "#F781BFFF", "#999999FF")
dark2 <- c("#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF", "#66A61EFF", "#E6AB02FF", "#A6761DFF", "#666666FF")
pair1 <- c("#A6CEE3FF", "#1F78B4FF", "#B2DF8AFF", "#33A02CFF", "#FB9A99FF", "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF", "#CAB2D6FF", "#6A3D9AFF", "#FFFF99FF", "#B15928FF")

column_anno <- HeatmapAnnotation(
    sample = sample_df$sample,
    group = sample_df$group,
    run = sample_df$run,
    col = list(
        # samples use color brewer Set1
        sample = c(
            "H146" = set1[1],
            "H1975" = set1[2],
            "H211" = set1[3],
            "H2228" = set1[4],
            "H526" = set1[5],
            "H69" = set1[6],
            "HCC827" = set1[7],
            "SHP77" = set1[8]
        ),

    group = c(
            "sclc_a" = dark2[1],
            "sclc_p" = dark2[2],
            "luad" = dark2[3]
        ),
        run = c(
            "run1" = pair1[3],
            "run2" = pair1[4]
        )
    )
)

col_fun = circlize::colorRamp2(c(0.8, 0.9, 1), c("blue", "white", "red"))

Heatmap(
    corr_matrix,
    col = col_fun,
    top_annotation = column_anno
)
