# Load necessary libraries
library(ggplot2)
library(scales)
# library(viridis)
library(gghalves)
library(patchwork)
library(dplyr)

# setup input and output directories
# List of file paths and corresponding sample names
file_paths <- unlist(snakemake@input)
sample_names <- snakemake@params$sample_id
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

# Specify datatype
all_data$datatype <- ifelse(grepl("pb", all_data$sample), "PacBio",
              ifelse(grepl("dRNA", all_data$sample), "ONT dRNA", "ONT cDNA")
)
# fix levels of datatype
all_data$datatype <- factor(all_data$datatype, levels = c("PacBio", "ONT dRNA", "ONT cDNA"))

# Boxplot for quals
p1 <- ggplot(all_data, aes(x = sample, y = quals, fill = datatype)) +
  geom_half_violin(trim = FALSE, side = "l" ,width = 1.8 ) + # Half violin plot
  geom_half_boxplot(outliers = FALSE, color = "darkgrey", side = "r", width = 0.5) + # Half boxplot
  theme_minimal() +
  labs(title = "Quality Scores Across Samples", x = "Sample", y = "Quality Score") +
  theme(
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = color_palette[c("PacBio", "ONT", "ONT_1")] |> unname())

# Boxplot for lengths
p2 <- ggplot(all_data, aes(x = sample, y = lengths, fill = datatype)) +
  geom_half_violin(trim = FALSE, side = "l",width = 1.8) + # Half violin plot
  geom_half_boxplot(outliers = FALSE, color = "darkgrey", side = "r", width = 0.5) + # Half boxplot
  theme_minimal() +
  labs(title = "Read Lengths Across Samples", x = "Sample", y = "Read Length") +
  theme(
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(0, 5000) +
  scale_fill_manual(values = color_palette[c("PacBio", "ONT", "ONT_1")] |> unname())

# arrange the plots in a grid and save
library(gridExtra)

rst <- p1 / (p2 + theme(legend.position = "none")) +
        patchwork::plot_layout(
          guides = "collect"
        ) & theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# save the plots
# increase the left margin

ggsave(output_fig, rst, width = 25, height = 13, units = "in", dpi = 150)


# save a summary table
summary_table <- all_data %>%
  group_by(sample, datatype) %>%
  summarise(
    mean_qual = mean(quals),
    median_qual = median(quals),
    mean_length = mean(lengths),
    median_length = median(lengths)
  ) %>%
  arrange(datatype, sample)

write.table(summary_table, snakemake@output[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
