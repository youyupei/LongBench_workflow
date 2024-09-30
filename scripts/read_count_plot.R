# config
bulk_count_file <- snakemake@input$bulk_read_count |> unlist()
sc_blaze_summary <- snakemake@input$sc_blaze_summary |> unlist()
bulk_sample_name <- snakemake@params$bulk_sample_name |> unlist()
sc_sample_names <- snakemake@params$sc_sample_name |> unlist()
output_fig <- snakemake@output[[1]]

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
df <- data.frame(
  sample = c(bulk_sample_name, sc_sample_names),
  read_count = c(bulk_read_counts, sc_read_counts)
)

# Plot barplot
library(ggplot2)
p <- ggplot(df, aes(x = sample, y = read_count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Read count comparison", x = "Sample", y = "Read count")

# Save the plot
ggsave(output_fig, p, width = 10, height = 6, units = "in", dpi = 300)