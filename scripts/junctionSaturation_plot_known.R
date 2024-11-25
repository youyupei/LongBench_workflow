rl_files <- c(snakemake@input$lr_length |> unlist(), snakemake@input$sr_length |> unlist())
scripts <- c(snakemake@input$lr_bulk |> unlist(), snakemake@input$sr_bulk |> unlist())
# scripts <- c("/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H146.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H69.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H526.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H211.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_SHP77.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H1975.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_H2228.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/ont_bulk_HCC827.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H146.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H69.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H526.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H211.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_SHP77.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H1975.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_H2228.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/pb_bulk_HCC827.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H146.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H69.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H526.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H211.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_SHP77.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H1975.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_H2228.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/lr_bulk/result/qc/RSeQC/dRNA_bulk_HCC827.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H146.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H69.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H526.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H211.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/SHP77.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H1975.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/H2228.junctionSaturation_plot.r", "/vast/projects/LongBench/analysis/sr_bulk/result/qc/RSeQC/HCC827.junctionSaturation_plot.r")
v_name <- c(snakemake@params$lr_sample_name |> unlist(), snakemake@params$sr_sample_name |> unlist())
read_number_table <- snakemake@input$read_number_table
read_length_table <- snakemake@input$read_length_table
variables <- list()
# read the read number table, get the average read count from rows whose "sample" matches the sample name
read_number_table <- read.table(read_number_table, header = TRUE, sep = "\t")
read_length_table <- read.table(read_length_table, header = TRUE, sep = "\t")
read_number <- read_number_table[match(v_name, read_number_table$sample), ]$read_count

# for each v_name, get the read length
read_length <- c()
for (i in v_name) {
  # if i start with "illumna", then read length is 100
  if (grepl("Illumina", i)) {
    read_length <- c(read_length, 100)
  } else {
    read_length <- c(read_length, read_length_table[read_length_table$sample == i, ]$mean_length)
  }
}





for (i in seq_along(scripts)) {
    # Read the third line from the file
    third_line <- readLines(scripts[i], n = 3)[3]
    
    # Use regular expression to extract the variable name and values
    var_name <- v_name[i]
    values <- sub(".*c\\((.*)\\)", "\\1", third_line)     # Extract values between c(...)
    
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


interval <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
df <- data.frame(variables) %>%
    mutate(x=interval) %>%
    melt(id.vars = "x") %>%
    rename(sample_names=variable) %>%
    mutate(datatype = ifelse(grepl("Illumina", sample_names), "Illumina",
        ifelse(grepl("pb", sample_names), "PacBio",
          ifelse(grepl("dRNA", sample_names), "ONT dRNA", "ONT cDNA")
        )
      )
    )


df$read_number <- read_number[match(df$sample_names, v_name)] * df$x / 100
df$total_bases <- df$read_number * read_length[match(df$sample_names, v_name)]

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
p1 <- ggplot(df, aes(x =x, y = value, color = datatype, group=sample_names)) + 
    geom_line(alpha=0.6) + 
    labs(title = "Saturation of known splice junctions", x = "Percentage of total reads", y = "Number of splicing junctions") + 
    theme_minimal()+ 
    scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)


p1_2 <- ggplot(df, aes(x = read_number, y = value, color = datatype, group = sample_names)) +
  geom_line(alpha = 0.6) +
  labs(title = "Saturation of known splice junctions", x = "Total reads", y = "Number of splicing junctions") +
  theme_minimal() +
  scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname())

p1_3 <- ggplot(df, aes(x = total_bases, y = value, color = datatype, group = sample_names)) +
  geom_line(alpha = 0.6) +
  labs(title = "Saturation of known splice junctions", x = "Total bases", y = "Number of splicing junctions") +
  theme_minimal() +
  scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname())

# Take average of samples with same datatype
df_avg <- df %>%
    group_by(datatype, x) %>%
    summarise(value = mean(value), read_number = mean(read_number), total_bases = mean(total_bases))

# Plot
p2 <- ggplot(df_avg, aes(x =  x, y = value, color = datatype)) + 
    geom_line(alpha=0.8, linewidth=2) +
    labs(title = "Saturation of known splice junctions", x = "Percentage of total reads", y = "Number of splicing junctions") + 
    theme_minimal()+ 
    scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)

p2_2 <- ggplot(df_avg, aes(x = read_number, y = value, color = datatype)) +
  geom_line(alpha = 0.8, linewidth = 2) +
  labs(title = "Saturation of known splice junctions", x = "Total reads", y = "Number of splicing junctions") +
  theme_minimal() +
  scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname())

p2_3 <- ggplot(df_avg, aes(x = total_bases, y = value, color = datatype)) +
  geom_line(alpha = 0.8, linewidth = 2) +
  labs(title = "Saturation of known splice junctions", x = "Total bases", y = "Number of splicing junctions") +
  theme_minimal() +
  scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname())

# 
# p1 <- p1_2 / p2_2
# 
# ggsave(snakemake@output[[1]], p, width = 10, height = 8)

# Open PDF device to save the plots as multipage
pdf(snakemake@output[[1]], width = 7, height = 5)
  print(p1 / p2)
  print(p1_2 / p2_2)
  print(p1_3 / p2_3)
# Close the PDF device
dev.off()