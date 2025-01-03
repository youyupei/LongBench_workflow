# config
bulk_count_file <- snakemake@input$bulk_read_count |> unlist()
sc_blaze_summary <- snakemake@input$sc_blaze_summary |> unlist()
sr_bulk_json <- paste0(snakemake@input$sr_bulk_salmon |> unlist(), "/aux_info/meta_info.json")
print(sr_bulk_json)
bulk_sample_name <- snakemake@params$bulk_sample_name |> unlist()
sc_sample_names <- snakemake@params$sc_sample_name |> unlist()
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

# Extract read counts from all files
bulk_read_counts <- sapply(bulk_count_file, read_count_from_file)
sc_read_counts <- sapply(sc_blaze_summary, extract_blaze_stats)
sr_bulk_read_counts <- sapply(sr_bulk_json, extract_sr_bulk_read_count)

# Create a data frame
sample_names <- c(bulk_sample_name, sc_sample_names, sr_bulk_sample_name)
df <- data.frame(
  sample = factor(sample_names),
  read_count = c(bulk_read_counts, sc_read_counts, sr_bulk_read_counts),
  datatype = ifelse(grepl("Illumina", sample_names), "Illumina",
    ifelse(grepl("pb", sample_names), "PacBio",
      ifelse(grepl("dRNA", sample_names), "ONT dRNA", "ONT cDNA")
    )
  )
)
df$datatype <- factor(df$datatype, levels = c("Illumina", "PacBio", "ONT dRNA", "ONT cDNA"))
# sort the data frame by datatype
df <- df[order(df$datatype), ]
# reorder the level of sample to match the datatype
# Define the custom prefix order
prefix_order <- c("Illumina", "pb", "dRNA", "ont")
# Sort levels by custom prefix order
sorted_levels <- unlist(lapply(prefix_order, function(p) {
  grep(paste0("^", p, "_*"), levels(df$sample), value = TRUE)
}))
# Reorder factor levels
df$sample <- factor(df$sample, levels = sorted_levels)

# Plot barplot
library(ggplot2)
p <- ggplot(df, aes(x = sample, y = read_count, fill = datatype)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Read count comparison", x = "Sample", y = "Read/Read pair count") +
  scale_fill_manual(values = color_palette[c("Illumina", "PacBio", "ONT", "ONT_1")] %>% unname)
# Save the plot
ggsave(output_fig, p, width = 10, height = 6, units = "in", dpi = 300)


# save a summary table
write.table(df, snakemake@output[[2]], sep = "\t", quote = FALSE, row.names = FALSE)