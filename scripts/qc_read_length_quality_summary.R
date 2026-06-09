library(dplyr)

file_paths   <- unlist(snakemake@input)
sample_names <- unlist(snakemake@params$sample_id)
output_table <- snakemake@output[[1]]

read_sample_data <- function(file_path, sample_name) {
  data <- read.table(gzfile(file_path), header = TRUE)
  data$sample <- sample_name
  data
}

all_data <- do.call(rbind, lapply(seq_along(file_paths), function(i) {
  read_sample_data(file_paths[i], sample_names[i])
}))

all_data$datatype <- ifelse(grepl("pb", all_data$sample), "PacBio",
  ifelse(grepl("dRNA", all_data$sample), "ONT dRNA",
  ifelse(grepl("^ill_", all_data$sample), "Illumina", "ONT cDNA")))
all_data$datatype <- factor(
  all_data$datatype,
  levels = c("PacBio", "ONT dRNA", "ONT cDNA", "Illumina")
)

whisker_low  <- function(x) { q1 <- quantile(x, 0.25); min(x[x >= q1 - 1.5 * IQR(x)]) }
whisker_high <- function(x) { q3 <- quantile(x, 0.75); max(x[x <= q3 + 1.5 * IQR(x)]) }

summary_table <- all_data |>
  group_by(sample, datatype) |>
  summarise(
    qual_min          = min(quals),
    qual_whisker_low  = whisker_low(quals),
    qual_Q1           = quantile(quals, 0.25),
    qual_median       = median(quals),
    qual_mean         = mean(quals),
    qual_Q3           = quantile(quals, 0.75),
    qual_whisker_high = whisker_high(quals),
    qual_max          = max(quals),
    length_min          = min(lengths),
    length_whisker_low  = whisker_low(lengths),
    length_Q1           = quantile(lengths, 0.25),
    length_median       = median(lengths),
    length_mean         = mean(lengths),
    length_Q3           = quantile(lengths, 0.75),
    length_whisker_high = whisker_high(lengths),
    length_max          = max(lengths),
    .groups = "drop"
  ) |>
  arrange(datatype, sample)

dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)
write.table(summary_table, output_table, sep = "\t", quote = FALSE, row.names = FALSE)
