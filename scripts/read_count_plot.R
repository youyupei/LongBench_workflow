# config
bulk_count_file <- snakemake@input$bulk_read_count |> unlist()
sc_blaze_summary <- snakemake@input$sc_blaze_summary |> unlist()
sr_bulk_json <- paste0(snakemake@input$sr_bulk_salmon |> unlist(), "/aux_info/meta_info.json")
sr_sc_out <- snakemake@input$sr_sc_out |> unlist()

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

# Plot barplot
library(ggplot2)
p <- ggplot(df, aes(x = sample, y = read_count/1000000, fill = datatype)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Read count comparison", x = "Datasets", y = "Read/Read pair counts (Millions)") +
  scale_fill_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)
# Save the plot
ggsave(output_fig, p, width = 10, height = 6, units = "in", dpi = 300)

# save a summary table
write.table(df, snakemake@output[[2]], sep = "\t", quote = FALSE, row.names = FALSE)