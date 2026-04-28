# config
bulk_count_file    <- snakemake@input$bulk_read_count |> unlist()
sc_blaze_summary   <- snakemake@input$sc_blaze_summary |> unlist()
sr_bulk_json       <- paste0(snakemake@input$sr_bulk_salmon |> unlist(), "/aux_info/meta_info.json")
sr_sc_out          <- snakemake@input$sr_sc_out |> unlist()
lr_bulk_base_count    <- snakemake@input$lr_bulk_base_count |> unlist()
sr_bulk_base_count    <- snakemake@input$sr_bulk_base_count |> unlist()
lr_sc_sn_base_count   <- snakemake@input$lr_sc_sn_base_count |> unlist()
sr_sc_sn_base_count   <- snakemake@input$sr_sc_sn_base_count |> unlist()

bulk_sample_name <- snakemake@params$bulk_sample_name |> unlist()
lr_sc_sample_name <- snakemake@params$lr_sc_sample_name |> unlist()
sr_sc_sample_name <- snakemake@params$sr_sc_sample_name |> unlist()
sr_bulk_sample_name <- snakemake@params$sr_sample_name |> unlist()

output_fig <- snakemake@output[[1]]

# color
color_palette <- c(
  PacBio = "#df1995",
  PacBio_1 = "#d5408d",
  ONT = "#00789b",
  ONT_1 = "#04476c",
  ONT_2 = "#24cdcd",
  Illumina = "#e88b20"
)

# Load necessary library
library(stringr)
library(tidyr)
library(tidyverse)

extract_sr_bulk_read_count <- function(json) {
  # Read the json file
  json_data <- jsonlite::fromJSON(json)
  # Extract the read count
  return(json_data$num_processed)
}

extract_sr_bulk_mapped_read_count <- function(json) {
  # Read the json file
  json_data <- jsonlite::fromJSON(json)
  # Extract the read count
  return(json_data$num_mapped)
}

extract_blaze_stats <- function(file_path) {
  # Read the file
  lines <- readLines(file_path, warn = FALSE)
  
  # Function to extract numbers and remove commas
  extract_number <- function(pattern, offset = 1) {
    match <- grep(pattern, lines)
    if (length(match) > 0) {
      num_str <- str_extract(lines[match + offset], "[\\d,]+")
      as.numeric(gsub(",", "", num_str))
    } else {
      return(NA)
    }
  }
  
  # Extract numbers using the custom function
  total_reads <- extract_number("Total number of reads:")
  return(total_reads)
}

read_count_from_file <- function(file_path) {
  return(as.numeric(readLines(file_path))) 
}

# SR cellranger output
extract_cellranger_stats <- function(file_path) {
  # Read the file
  sc_out <- file.path(file_path, "outs/metrics_summary.csv")
  n <- read.csv(sc_out, header = TRUE, sep = ",")$Number.of.Reads %>%
    gsub(",", "", .) %>%
    as.numeric()
  return(n)
  }

# Extract read counts from all files
bulk_read_counts <- sapply(bulk_count_file, read_count_from_file)
sc_read_counts <- sapply(sc_blaze_summary, extract_blaze_stats)
sr_bulk_read_counts <- sapply(sr_bulk_json, extract_sr_bulk_read_count)
sr_sc_read_counts <- sapply(sr_sc_out, extract_cellranger_stats)

# Create a data frame
sample_names <- c(bulk_sample_name, lr_sc_sample_name, sr_bulk_sample_name, sr_sc_sample_name)
df <- data.frame(
  sample = factor(sample_names),
  read_count = c(bulk_read_counts, sc_read_counts, sr_bulk_read_counts, sr_sc_read_counts),
  datatype = ifelse(grepl("Illumina", sample_names), "Illumina",
    ifelse(grepl("ill", sample_names), "Illumina",
      ifelse(grepl("pb", sample_names), "PacBio",
        ifelse(grepl("dRNA", sample_names), "ONT dRNA", "ONT cDNA")
      )
    )
  )
) %>% mutate(sample = case_when(
  sample=="ill_sc" ~ "Illumina_sc",
  sample=="ill_sn" ~ "Illumina_sn",
  TRUE ~ sample
))


df$datatype <- factor(df$datatype, levels = c("Illumina", "PacBio", "ONT dRNA", "ONT cDNA"))
# sort the data frame by datatype
df <- df[order(df$datatype), ]
# reorder the level of sample to match the datatype
# Define the custom prefix order
prefix_order <- c("Illumina", "pb", "dRNA", "ont")
# Sort levels by custom prefix order
sorted_levels <- unlist(lapply(prefix_order, function(p) {
  grep(paste0("^", p, "_*"), unique(df$sample), value = TRUE)
}))
# Reorder factor levels
df$sample <- factor(df$sample, levels = sorted_levels)
df$celllines <- purrr::map_chr(df$sample, ~str_split(.x, "_")[[1]] %>% tail(1))
df <- df %>% mutate(
  library_type = case_when(
    grepl("sc", celllines) ~ "Single Cell",
    grepl("sn", celllines) ~ "Single Nuclei",
    TRUE ~ "Bulk"
  )
)

# Build base-count lookup: raw filename stem -> total bases
read_base_count_files <- function(paths) {
  stems  <- sub("\\.total_bases$", "", basename(paths))
  counts <- sapply(paths, function(f) as.numeric(readLines(f)))
  setNames(counts, stems)
}

all_base_counts <- c(
  read_base_count_files(lr_bulk_base_count),
  read_base_count_files(sr_bulk_base_count),
  read_base_count_files(lr_sc_sn_base_count),
  read_base_count_files(sr_sc_sn_base_count)
)

# Map raw stems to df$sample names
# sr_bulk files are named ill_bulk_{cell_line} -> Illumina_{cell_line}
# sr_sc_sn files are named ill_sc/ill_sn       -> Illumina_sc/Illumina_sn
names(all_base_counts) <- dplyr::case_when(
  grepl("^ill_bulk_", names(all_base_counts)) ~ sub("^ill_bulk_", "Illumina_", names(all_base_counts)),
  names(all_base_counts) == "ill_sc"          ~ "Illumina_sc",
  names(all_base_counts) == "ill_sn"          ~ "Illumina_sn",
  TRUE                                         ~ names(all_base_counts)
)

df$gb <- all_base_counts[as.character(df$sample)] / 1e9

# Plot barplot — read counts
library(ggplot2)
library(patchwork)

bar_theme <- theme(axis.text.x = element_text(angle = 60, hjust = 1))
fill_scale <- scale_fill_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)

p1 <- ggplot(df, aes(x = sample, y = read_count / 1e6, fill = datatype)) +
  geom_bar(stat = "identity") +
  bar_theme +
  labs(title = "Read count comparison", x = NULL, y = "Read/Read pair counts (Millions)") +
  fill_scale

p2 <- ggplot(df, aes(x = sample, y = gb, fill = datatype)) +
  geom_bar(stat = "identity") +
  bar_theme +
  labs(title = "Gigabases sequenced", x = "Datasets", y = "Gigabases (Gb)") +
  fill_scale

combined <- p1 / p2 + patchwork::plot_layout(guides = "collect")

ggsave(output_fig, combined, width = 10, height = 10, units = "in", dpi = 300)

# Dual y-axis: normalize both metrics independently to [0, max] of their own range
# so each axis has its own natural scale and tick labels are meaningful.
max_reads_M <- max(df$read_count / 1e6, na.rm = TRUE)
max_gb      <- max(df$gb,           na.rm = TRUE)

# Bars on primary axis (reads in millions, natural scale).
# Points plotted as fraction of max_gb, rescaled to primary axis range.
to_primary   <- function(gb_val)  gb_val / max_gb * max_reads_M
from_primary <- function(y_val)   y_val  / max_reads_M * max_gb

p_dual <- ggplot(df, aes(x = sample)) +
  geom_bar(aes(y = read_count / 1e6, fill = datatype), stat = "identity") +
  geom_point(aes(y = to_primary(gb), color = datatype), size = 2.5) +
  geom_line(aes(y = to_primary(gb), group = 1), color = "grey40", linewidth = 0.4) +
  scale_y_continuous(
    name = "Read/Read pair counts (Millions)",
    sec.axis = sec_axis(transform = from_primary, name = "Gigabases (Gb)")
  ) +
  fill_scale +
  scale_color_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname,
                     guide = "none") +
  bar_theme +
  labs(title = "Read counts and gigabases sequenced", x = "Datasets", fill = "Platform") +
  theme(axis.title.y.right = element_text(color = "grey40"),
        axis.text.y.right  = element_text(color = "grey40"))

dual_fig <- sub("\\.svg$", "_dual_axis.svg", output_fig)
ggsave(dual_fig, p_dual, width = 12, height = 5, units = "in", dpi = 300)

# save a summary table
write.table(df, snakemake@output[[2]], sep = "\t", quote = FALSE, row.names = FALSE)