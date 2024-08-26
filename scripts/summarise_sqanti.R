# This script contains code modified from the SQANTI3 github under GPLv3 license
        ### Author: Lorena de la Fuente, Elizabeth Tseng & Francisco J Pardo-Palacios
        ### Last Modified: 09/23/2020 by etseng@pacb.com

#********************** Taking arguments from python script

# take multiple sqanti folders and input
# dirs <- snakemake@input |> unlist()
# sample_names <- snakemake@params$sample_names
# utilities.path <- snakemake@params$utilities_path
# output <- snakemake@output[[1]]
# test
dirs <- c("/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/qc/sqanti3/ont_sn",
            "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/qc/sqanti3/ont_sc",
            "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/qc/sqanti3/ont_sn_clean",
            "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/qc/sqanti3/ont_sc_clean",
            "/home/users/allstaff/you.yu/LongBench/analysis/lr_bulk/result/qc/sqanti3/ont_bulk_H146"
            )
sample_names <- c("ont_sn", "ont_sc", "ont_sn_clean", "ont_sc_clean", "ont_bulk_H146")
utilities.path <- '/home/users/allstaff/you.yu/LongBench/software/SQANTI3/utilities'
output <- '~/test.pdf'

parse_sqanti_ourdir <- function(dir){
  # get the classification file and junction file
  # search for file end with classification.txt and junctions.txt
  files <- list.files(dir, full.names=TRUE)
  class.file <- files[grep("classification.txt", files)]
  junction.file <- files[grep("junctions.txt", files)]
  return(list(class.file=class.file, junction.file=junction.file))
}


get_structural_category_table <- function(class.file, sample.name){
    data.class <- read.table(class.file, header=T, as.is=T, sep="\t")
    rownames(data.class) <- data.class$isoform

    xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
    xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
    subc.levels=c("alternative_3end",'alternative_3end5end', "alternative_5end","reference_match", "3prime_fragment","internal_fragment", "5prime_fragment","combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention","no_combination_of_known_junctions", "mono-exon_by_intron_retention", "at_least_one_novel_splicesite", "mono-exon", "multi-exon")
    subc.labels=c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end", "Reference match", "3' fragment", "Internal fragment", "5' fragment", "Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Not comb. of annot. junctions", "Mono-exon by intron ret.", "At least 1 annot. don./accept.", "Mono-exon", "Multi-exon")
    coding.levels=c("coding", "non_coding")
    coding.labels=c("Coding", "Non coding")

    data.class$structural_category = factor(data.class$structural_category,
                                            labels = xaxislabelsF1,
                                            levels = xaxislevelsF1,
                                            ordered=TRUE)
    data.class$subcategory = factor(data.class$subcategory,
                                            labels = subc.labels,
                                            levels = subc.levels,
                                            ordered=TRUE)
    data.class$coding = factor(data.class$coding,
                                    labels = coding.labels,
                                    levels = coding.levels,
                                    ordered=TRUE)

    # Table 1: Number of isoforms in each structural category
    table1MD <- as.data.frame(table(data.class$structural_category))
    rownames(table1MD) <- table1MD[,1]
    table1MD <- table1MD['Freq']
    colnames(table1MD) <- c(sample.name)
    return(table1MD)
}

structural_category_table <- do.call(cbind, lapply(1:length(dirs), function(i){
    dir <- dirs[i]
    sample.name <- sample_names[i]
    files <- parse_sqanti_ourdir(dir)
    class.file <- files$class.file
    table1MD <- get_structural_category_table(class.file, sample.name)
    return(table1MD)
}))


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Convert row names to a column before pivoting
structural_category_table <- structural_category_table %>%
  rownames_to_column(var = "Category")

# Convert to long format
long_df <- pivot_longer(structural_category_table, cols = -Category, names_to = "Sample", values_to = "Count")

# Add a 'Group' column to differentiate FSM and Others
long_df <- long_df %>% group_by(Sample) %>%
  mutate("Full splice match" = ifelse(Category == "FSM", "Yes", "No"), Proportion = Count / sum(Count))

# Make FSM last to be on top and Yes to be on top of No
long_df$Category <- factor(long_df$Category, levels = c("ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron", "FSM"))
long_df$`Full splice match` <- factor(long_df$`Full splice match`, levels = c("No", "Yes"))
# Define the fill colors to simulate hierarchy
fill_colors <- c("Yes" = "#FF9100", "No" = "#e2dbdb00")
# Plot
p <- ggplot(long_df, aes(x = Sample, y = Proportion)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = Category)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 3, alpha = 0, aes(color = `Full splice match`)) +
  guides(fill = guide_legend(title = "Structural Category", 
                             title.position = "top", 
                             title.hjust = 0.5,
                             keyheight=2)) + 
    guides(color = guide_legend(
    title = "Full splice match", 
    title.position = "top",
    title.hjust = 0.5,
    override.aes = list(fill = "#baaaaa", alpha = 1),
    keyheight=2
    ))+
  scale_fill_viridis_d() +
  scale_color_manual(values = fill_colors) +
  theme_minimal(base_size = 15) +
  labs(title = "Proportion of Structural Categories Across Samples",
       x = "Sample",
       y = "Proportion")

# Save the plot
ggsave(output, plot = p, width = 10, height = 10, dpi = 300)