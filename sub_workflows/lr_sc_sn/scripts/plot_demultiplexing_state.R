# config
file_paths <- snakemake@input
sample_names <- snakemake@params[["sample_names"]]
output_fig <- snakemake@output[[1]]

# Load necessary library


library(stringr)

extract_stats <- function(file_path) {
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
  unambiguous_reads <- extract_number("Reads with unambiguous polyT and adapter positions found:")
  high_quality_bc <- extract_number("Reads with unambiguous polyT and adapter positions found:", offset = 2)
  no_polyT_adapter <- extract_number("no polyT and adapter positions found:")
  both_end_fail <- extract_number("polyT and adapter positions found in both end")
  multiple_polyT_adapter <- extract_number("multiple polyT and adapter found in one end")
  identified_cells <- extract_number("Identified # of cells:", 0)
  reads_in_cells <- extract_number("Total reads in cells:", 0)
  
  # Create a dataframe row
  data.frame(
    Total_Reads = total_reads,
    Unambiguous_Reads = unambiguous_reads,
    High_Quality_BC = high_quality_bc,
    No_PolyT_Adapter = no_polyT_adapter,
    Both_End_Fail = both_end_fail,
    Multiple_PolyT_Adapter = multiple_polyT_adapter,
    Identified_Cells = identified_cells,
    Reads_In_Cells = reads_in_cells,
    stringsAsFactors = FALSE
  )
}

# List all text files in the directory


# Apply the function to all files and combine the results into a dataframe
result_df <- do.call(rbind, lapply(file_paths, extract_stats))
rownames(result_df) <- sample_names




# Plot
library(ggplot2)
library(reshape2)

# Prepare data for the plot
plot_data <- data.frame(
  Sample = rownames(result_df),
  Reads_In_Cells = result_df$Reads_In_Cells,
  Other_Reads = result_df$Total_Reads - result_df$Reads_In_Cells
)

# Convert to long format for ggplot
plot_data_melt <- melt(plot_data, id.vars = "Sample")

# Reorder factor levels to flip the order
plot_data_melt$variable <- factor(plot_data_melt$variable, levels = c("Other_Reads","Reads_In_Cells"))

# Create stacked bar plot with flipped order
plot <- ggplot(plot_data_melt, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Reads and Reads in Cells", x = "Sample", y = "Number of Reads") +
  scale_fill_manual(values = c("Other_Reads" = "lightgray", "Reads_In_Cells" = "lightblue"), 
                    name = "Reads", labels = c("Other Reads", "Reads in Cells")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 32, face = "bold"),      
    axis.title.x = element_text(size = 28),                  
    axis.title.y = element_text(size = 28),               
    axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 24),                   
    legend.title = element_text(size = 24),                
    legend.text = element_text(size = 24)                       
  )

ggsave(output_fig, plot = plot, width = 20, height = 6)