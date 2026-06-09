library(dplyr)
library(stringr)
library(jsonlite)

# inputs
bulk_count_file     <- unlist(snakemake@input$bulk_read_count)
sc_blaze_summary    <- unlist(snakemake@input$sc_blaze_summary)
sr_bulk_json        <- paste0(unlist(snakemake@input$sr_bulk_salmon), "/aux_info/meta_info.json")
sr_sc_out           <- unlist(snakemake@input$sr_sc_out)
lr_bulk_base_count  <- unlist(snakemake@input$lr_bulk_base_count)
sr_bulk_base_count  <- unlist(snakemake@input$sr_bulk_base_count)
lr_sc_sn_base_count <- unlist(snakemake@input$lr_sc_sn_base_count)
sr_sc_sn_base_count <- unlist(snakemake@input$sr_sc_sn_base_count)

bulk_sample_name    <- unlist(snakemake@params$bulk_sample_name)
lr_sc_sample_name   <- unlist(snakemake@params$lr_sc_sample_name)
sr_sc_sample_name   <- unlist(snakemake@params$sr_sc_sample_name)
sr_bulk_sample_name <- unlist(snakemake@params$sr_sample_name)

output_table <- snakemake@output[[1]]

# extraction helpers
extract_blaze_total <- function(f) {
  lines <- readLines(f, warn = FALSE)
  idx <- grep("Total number of reads:", lines)
  if (length(idx) == 0) return(NA_real_)
  as.numeric(gsub(",", "", str_extract(lines[idx + 1], "[\\d,]+")))
}

extract_cellranger_reads <- function(dir) {
  f <- file.path(dir, "outs/metrics_summary.csv")
  as.numeric(gsub(",", "", read.csv(f)$Number.of.Reads))
}

read_base_count_files <- function(paths) {
  stems  <- sub("\\.total_bases$", "", basename(paths))
  counts <- sapply(paths, function(f) as.numeric(readLines(f, warn = FALSE)))
  setNames(counts, stems)
}

# read counts per source
bulk_counts   <- sapply(bulk_count_file,  function(f) as.numeric(readLines(f, warn = FALSE)))
sc_counts     <- sapply(sc_blaze_summary, extract_blaze_total)
sr_bulk_counts <- sapply(sr_bulk_json,   function(f) fromJSON(f)$num_processed)
sr_sc_counts  <- sapply(sr_sc_out,       extract_cellranger_reads)

# combined read count data frame
sample_names <- c(bulk_sample_name, lr_sc_sample_name, sr_bulk_sample_name, sr_sc_sample_name)
df <- data.frame(
  sample     = sample_names,
  read_count = c(bulk_counts, sc_counts, sr_bulk_counts, sr_sc_counts),
  stringsAsFactors = FALSE
) |>
  mutate(
    datatype = case_when(
      grepl("^Illumina", sample) | grepl("^ill_", sample) ~ "Illumina",
      grepl("pb",  sample) ~ "PacBio",
      grepl("dRNA", sample) ~ "ONT dRNA",
      TRUE ~ "ONT cDNA"
    ),
    datatype = factor(datatype, levels = c("Illumina", "PacBio", "ONT dRNA", "ONT cDNA")),
    library_type = case_when(
      grepl("_sc$", sample) ~ "Single Cell",
      grepl("_sn$", sample) ~ "Single Nuclei",
      TRUE ~ "Bulk"
    )
  )

# gigabases
all_base_counts <- c(
  read_base_count_files(lr_bulk_base_count),
  read_base_count_files(sr_bulk_base_count),
  read_base_count_files(lr_sc_sn_base_count),
  read_base_count_files(sr_sc_sn_base_count)
)
names(all_base_counts) <- case_when(
  grepl("^ill_bulk_", names(all_base_counts)) ~ sub("^ill_bulk_", "Illumina_", names(all_base_counts)),
  names(all_base_counts) == "ill_sc"          ~ "Illumina_sc",
  names(all_base_counts) == "ill_sn"          ~ "Illumina_sn",
  TRUE                                         ~ names(all_base_counts)
)
df$gb <- all_base_counts[df$sample]

df <- df |>
  arrange(datatype, sample) |>
  select(sample, datatype, library_type, read_count, gb)

dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)
write.table(df, output_table, sep = "\t", quote = FALSE, row.names = FALSE)
