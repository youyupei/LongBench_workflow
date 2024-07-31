library(FLAMES)
library(dplyr)

# take command line input for input bam file, input allele stat file, output file and threads
args <- commandArgs(trailingOnly = TRUE)
bam_path <- args[1]
allele_stat_path <- args[2]
barcode_list <- args[3]
threads <- as.numeric(args[4])
output_rds <- args[5]



all_mut <- read.csv(
  gzfile(allele_stat_path), 
  comment.char='#')

snps_df <- FLAMES::sc_mutations(
  bam_path = bam_path,
  seqnames = all_mut %>% .$ch,
  positions = all_mut %>% .$pos+1,
  indel = F,
  barcodes = read.csv(barcode_list, header=F)$V1,
  threads = threads)


saveRDS(snps_df, file = paste0(output_rds, '_snps_df.Rds'))
# create a mutation assay
mut_df <- all_mut
result_list <- list()
for (i in 1:nrow(mut_df)) {
  ref <- mut_df[i,]$REF
  alt <- mut_df[i,]$ALT
  g_pos <- mut_df[i,]$pos+1
  ref_d <- snps_df %>% 
    subset(pos==g_pos & allele == ref) %>% 
    mutate(!!paste0("REF_", mut_df[i,]$gene_name,'_', g_pos):=allele_count) %>% 
    dplyr::select(barcode, !!paste0("REF_", mut_df[i,]$gene_name,'_', g_pos))
  
  alt_d <- snps_df %>% 
    subset(pos==g_pos & allele == alt) %>% 
    mutate(!!paste0("ALT_", mut_df[i,]$gene_name,'_', g_pos):=allele_count) %>% 
    dplyr::select(barcode, !!paste0("ALT_", mut_df[i,]$gene_name,'_', g_pos))
  result_list <- append(result_list, list(ref_d))
  result_list <- append(result_list, list(alt_d))
  
} 
mut_count <- purrr::reduce(result_list, left_join, by="barcode") %>% as.data.frame()
rownames(mut_count) <- mut_count$barcode
mut_count$barcode <- NULL

saveRDS(mut_count, file = output_rds)
