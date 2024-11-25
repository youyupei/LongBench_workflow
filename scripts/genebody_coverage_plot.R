# scripts <- c("/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H146.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H69.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H526.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H211.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_SHP77.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H1975.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H2228.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_HCC827.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H146.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H69.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H526.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H211.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_SHP77.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H1975.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H2228.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_HCC827.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H146.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H69.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H526.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H211.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_SHP77.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H1975.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H2228.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_HCC827.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_sc_sn/result/qc/RSeQC/pb_sc.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_sc_sn/result/qc/RSeQC/pb_sn.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_sc_sn/result/qc/RSeQC/ont_sn.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/lr_sc_sn/result/qc/RSeQC/ont_sc.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H146.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H69.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H526.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H211.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/SHP77.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H1975.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H2228.geneBodyCoverage.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/HCC827.geneBodyCoverage.r")

scripts <- c(snakemake@input$lr_bulk |> unlist(), snakemake@input$sr_bulk |> unlist())
v_name <- c(snakemake@params$lr_sample_name |> unlist(), snakemake@params$sr_sample_name |> unlist())

# v_name <- c("ont_bulk_H146", "ont_bulk_H69", "ont_bulk_H526", "ont_bulk_H211", "ont_bulk_SHP77", "ont_bulk_H1975", "ont_bulk_H2228", "ont_bulk_HCC827", "pb_bulk_H146", "pb_bulk_H69", "pb_bulk_H526", "pb_bulk_H211", "pb_bulk_SHP77", "pb_bulk_H1975", "pb_bulk_H2228", "pb_bulk_HCC827", "dRNA_bulk_H146", "dRNA_bulk_H69", "dRNA_bulk_H526", "dRNA_bulk_H211", "dRNA_bulk_SHP77", "dRNA_bulk_H1975", "dRNA_bulk_H2228", "dRNA_bulk_HCC827", "pb_sc", "pb_sn","ont_sc", "ont_sn", "Illumina_H146", "Illumina_H69", "Illumina_H526", "Illumina_H211", "Illumina_SHP77", "Illumina_H1975", "Illumina_H2228", "Illumina_HCC827")
variables <- list()
for (i in seq_along(scripts)) {
    # Read the first line from the file
    first_line <- readLines(scripts[i], n = 1)
    
    # Use regular expression to extract the variable name and values
    var_name <- v_name[i]
    values <- sub(".*c\\((.*)\\)", "\\1", first_line)     # Extract values between c(...)
    
    # Convert values to a numeric vector
    values <- as.numeric(strsplit(values, ",")[[1]])
    
    # Assign to the variables list
    variables[[var_name]] <- values
}

# Plot a combined graph with each variable as a curve using ggplot2
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# color
color_palette <- c(
  PacBio = "#df1995",
  PacBio_1 = "#d5408d",
  ONT = "#00789b",
  ONT_1 = "#04476c",
  ONT_2 = "#24cdcd",
  Illumina = "#e88b20"
)


df <- data.frame(variables) %>%
    mutate(Percentile = seq(1, nrow(.))) %>%
    melt(id.vars = "Percentile") %>%
    rename(sample_names=variable) %>%
    mutate(datatype = ifelse(grepl("Illumina", sample_names), "Illumina",
        ifelse(grepl("pb", sample_names), "PacBio",
          ifelse(grepl("dRNA", sample_names), "ONT dRNA", "ONT cDNA")
        )
      )
    )


# Reorder factor levels
df$sample <- factor(df$sample)
df$datatype <- factor(df$datatype, levels = c("Illumina", "PacBio", "ONT dRNA", "ONT cDNA"))
df$sample_names <- factor(df$sample_names)

# order the samples
prefix_order <- c("Illumina", "pb", "dRNA", "ont")
# Sort levels by custom prefix order
sorted_levels <- unlist(lapply(prefix_order, function(p) {
  grep(paste0("^", p, "_*"), levels(df$sample), value = TRUE)
}))
df$sample <- factor(df$sample, levels = sorted_levels)


# Plot
p1 <- ggplot(df, aes(x = Percentile, y = value, color = datatype, group=sample_names)) + 
    geom_line(alpha=0.6) + 
    labs(title = "Gene Body Coverage", x = "Gene body percentile (5'->3')", y = "Coverage") + 
    theme_minimal()+ 
    scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)

# Take average of samples with same datatype
df_avg <- df %>%
    group_by(datatype, Percentile) %>%
    summarise(value = mean(value))

# Plot
p2 <- ggplot(df_avg, aes(x = Percentile, y = value, color = datatype)) + 
    geom_line(alpha=0.8, linewidth=2) +
    labs(title = "Gene Body Coverage", x = "Gene body percentile (5'->3')", y = "Coverage") + 
    theme_minimal()+ 
    scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)

# Save the plot
p <- p1 / p2
ggsave(snakemake@output[[1]], p1/p2, width = 10, height = 8)