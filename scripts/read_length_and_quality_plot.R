# Load necessary libraries
library(ggplot2)
library(scales)
library(viridis)
library(gghalves)
# setup input and output directories
# List of file paths and corresponding sample names
file_paths <- unlist(snakemake@input)
sample_names <- snakemake@params$sample_id
output_fig <- snakemake@output[[1]]

# create directory if not exist
dir.create(dirname(output_fig), recursive = TRUE, showWarnings = FALSE)
# Function to read data from gzipped files
read_sample_data <- function(file_path, sample_name) {
  data <- read.table(gzfile(file_path), header = TRUE)
  data$sample <- sample_name
  return(data)
}



# Read and combine data
all_data <- do.call(rbind, lapply(seq_along(file_paths), function(i) {
  read_sample_data(file_paths[i], sample_names[i])
}))

# Boxplot for quals
p1 <- ggplot(all_data, aes(x = sample, y = quals, fill = sample)) +
  geom_half_violin(trim = FALSE, side = "l") +  # Half violin plot
  geom_half_boxplot(outliers = FALSE, color = 'darkgrey', side = "r") +  # Half boxplot
  theme_minimal() +
  labs(title = "Quality Scores Across Samples", x = "Sample", y = "Quality Score") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()

# Boxplot for lengths
p2 <- ggplot(all_data, aes(x = sample, y = lengths, fill = sample)) +
  geom_half_violin(trim = FALSE, side = "l") +  # Half violin plot
  geom_half_boxplot(outliers = FALSE, color = 'darkgrey', side = "r") +  # Half boxplot
  theme_minimal() +
  labs(title = "Read Lengths Across Samples", x = "Sample", y = "Read Length") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,5000)
  scale_fill_viridis_d()

# arrange the plots in a grid and save
library(gridExtra)


rst <- grid.arrange( p1,p2, ncol=2)

# save the plots
# increase the left margin

ggsave(output_fig, rst, width = 32, height = 8, units = "in", dpi = 150)



