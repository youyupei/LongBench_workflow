library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(purrr)
library(pbapply)

# Load GTF
gtf_file <- "/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf"
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

overlaps_tibble$exon_gene <- sub("\\..*","", overlaps_tibble$exon_gene)
overlaps_tibble$intron_gene <- sub("\\..*","", overlaps_tibble$intron_gene)






# Build intron ranges for each gene (exons -> introns)
introns_by_gene <- lapply(exons_by_gene, function(gr) {
  gr <- reduce(gr)  # merge overlapping exons
  gaps <- gaps(gr)  # non-exonic regions within the span
  gaps[strand(gr)[1]]  # keep strand
})

# Convert to GRanges with gene_id
exons_all <- unlist(exons_by_gene, use.names = FALSE)
introns_all <- unlist(introns_by_gene, use.names = FALSE)
exons_all$gene_id <- rep(names(exons_by_gene), lengths(exons_by_gene))
introns_all$gene_id <- rep(names(introns_by_gene), lengths(introns_by_gene))

## 1. Genes whose exons overlap other genes' introns
hits_exon_intron <- findOverlaps(exons_all, introns_all, ignore.strand = TRUE)
overlap_exon_intron <- tibble(
  gene1 = exons_all$gene_id[queryHits(hits_exon_intron)],
  gene2 = introns_all$gene_id[subjectHits(hits_exon_intron)]
) %>% filter(gene1 != gene2)

## 2. Genes whose exons overlap other genes' exons
hits_exon_exon <- findOverlaps(exons_all, exons_all, ignore.strand = TRUE)
overlap_exon_exon <- tibble(
  gene1 = exons_all$gene_id[queryHits(hits_exon_exon)],
  gene2 = exons_all$gene_id[subjectHits(hits_exon_exon)]
) %>% filter(gene1 != gene2)

## 3. Genes whose exons overlap any part of other genes (exon or intron)
hits_exon_any <- findOverlaps(exons_all, c(exons_all, introns_all), ignore.strand = TRUE)
overlap_exon_any <- tibble(
  gene1 = exons_all$gene_id[queryHits(hits_exon_any)],
  gene2 = c(exons_all$gene_id, introns_all$gene_id)[subjectHits(hits_exon_any)]
) %>% filter(gene1 != gene2)