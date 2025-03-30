# read from file:
library(tidyverse)

# Get all file names
file_paths <- snakemake@input %>% unlist
# Read all files into a tibble
merged_data <- map_dfr(file_paths, ~ read_tsv(.x, col_names = FALSE) %>%
    mutate(file = basename(.x)))


colnames(merged_data) <- c("Count", "file_name")
merged_data <- merged_data %>% 
  separate(file_name, into = c('platform', "protocol", "cell_line"), sep = "_", extra = "drop") %>% 
  mutate(sample = paste(platform, protocol, cell_line, sep = "_"))


# add subsample information

merged_data %>%
    group_by(platform, protocol) %>%
    mutate(
        count.proportion = Count / sum(Count)
    ) %>% write_csv(snakemake@output[[1]])