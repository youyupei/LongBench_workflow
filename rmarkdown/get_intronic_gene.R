library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(purrr)
library(pbapply)
library(tibble)
gtf_file <- snakemake@input[["gtf"]]
annotations <- rtracklayer::import(gtf_file)

# Filter for gene and exon annotations
genes <- annotations[annotations$type == "gene"]
exons <- annotations[annotations$type == "exon"]

# for each gene, get the introns
exon_split <- split(exons, exons$gene_id)
library(furrr)
# Set up parallel processing with future_map (e.g., using all available cores)
plan(multisession)
options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB

introns  <- exon_split %>% names %>% 
    future_map(~ {
        intr <- gaps(exon_split[[.x]])
        if (length(intr) == 1) {
            return(NULL)
        } else {
            intr$gene_name <- .x
            return(intr[-1])
        }
    }, .progress = TRUE)

merged_introns <- do.call(c, introns)
overlaps <- findOverlaps(exons, merged_introns, ignore.strand=TRUE, type='within',  select="all", minoverlap =50)
overlaps_tibble <- tibble(
    exon_gene = exons$gene_id[queryHits(overlaps)],
    intron_gene = merged_introns$gene_name[subjectHits(overlaps)],
) %>%
    filter(exon_gene != intron_gene) %>%
    distinct()


# gene exon count
exons_count <-  exon_split %>% future_map(length,  .progress = TRUE) %>% unlist %>% enframe(name = "Gene", value = "Exon_Count")

saveRDS(list(
    exon_intron_overlaps = overlaps_tibble,
    exons_count = exons_count
), file = snakemake@output[["rds"]])