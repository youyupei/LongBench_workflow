# Load the Rsubread package
library(Rsubread)


# Path to the annotation file
assembly <- snakemake@input[["gtf"]]
index <- snakemake@input[["index"]]
forward_reads <- snakemake@input[["R1"]]
reverse_reads <- snakemake@input[["R2"]]
threads <- snakemake@resources[['cpus_per_task']]
output_bams <- snakemake@output[["bam"]]
# make sure the output directory exists
dir.create(dirname(output_bams), showWarnings = FALSE, recursive = TRUE)


# Perform alignment using Rsubread's subjunc function
align_files <- subjunc(
  index = index,           # Ensure the path to the index is correct
  readfile1 = forward_reads,
  readfile2 = reverse_reads,
  output_file = output_bams,
  nthreads = threads,                   # Adjust the number of threads as needed
  annot.ext = assembly,
  reportAllJunctions = TRUE,
  useAnnotation = TRUE,
  input_format = "gzFASTQ",
  isGTF = TRUE,
  unique = TRUE
)

# Save the alignment result to an RDS file
# saveRDS(align_files, 'aligned_files.RDS')

# Print a message indicating that the process is complete
# cat("Alignment complete. Aligned files saved to 'aligned_files.RDS'.\n")
