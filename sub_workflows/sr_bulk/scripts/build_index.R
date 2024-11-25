# Load Rsubread library
library(Rsubread)

# Path to the reference sequences in FASTA format
reference_sequences <- snakemake@input[["genome"]]
output_dir <- snakemake@output[["index"]]

# Check if the reference file exists
if (file.exists(reference_sequences)) {
  # make sure the output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # Build index using file.path to construct the full path to the index basename
  buildindex(basename = output_dir, reference = reference_sequences)
} else {
  stop("Reference sequence file does not exist at the specified path.")
}
