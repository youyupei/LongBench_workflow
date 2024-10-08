# config
bulk_count_file <- snakemake@input$bulk_read_count |> unlist()
sc_blaze_summary <- snakemake@input$sc_blaze_summary |> unlist()
bulk_sample_name <- snakemake@params$bulk_sample_name |> unlist()
sc_sample_names <- snakemake@params$sc_sample_name |> unlist()
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

# Create a data frame
sample_names <- c(bulk_sample_name, sc_sample_names)
datatype <- ifelse(grepl("pb", sample_names), "PacBio",
              ifelse(grepl("dRNA", sample_names), "ONT dRNA", "ONT cDNA")
)
# fix levels of datatype
datatype <- factor(datatype, levels = c("PacBio", "ONT dRNA", "ONT cDNA"))

df <- data.frame(
  sample = sample_names,
  read_count = c(bulk_read_counts, sc_read_counts),
  datatype = datatype
)

# Plot barplot
library(ggplot2)
p <- ggplot(df, aes(x = sample, y = read_count, fill = datatype)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Read count comparison", x = "Sample", y = "Read count") +
  scale_fill_manual(values = color_palette[c("PacBio", "ONT", "ONT_1")] %>% unname)
# Save the plot
ggsave(output_fig, p, width = 10, height = 6, units = "in", dpi = 300)