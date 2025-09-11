library(tidyverse)
library(dplyr)
library(tidyr)

params <- list(
    random_seed = 2024,
    sirv_ercc_gtf = '/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf',
    sequins_gtf = '/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_annotation_2.4.gtf',
    sequins_tsv = '/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_isoforms_2.4.tsv',
    human_gtf = '/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf',
    output.rds = '/home/users/allstaff/you.yu/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds'
)

calcTxNum <- function(G.Tx.map){
  geneid = rep(NA, nrow(G.Tx.map))
  txcount = table(G.Tx.map$gene_id)
  txnum = as.numeric(txcount[match(G.Tx.map$gene_id, names(txcount))])
}
# SIRV
SIRV_tx <- GenomicFeatures::makeTxDbFromGFF(params$sirv_ercc_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
SIRV.G.Tx.map <- tibble(
  gene_id = SIRV_tx$gene_id %>% unlist,
  tx_name = SIRV_tx$tx_name
)
rm(SIRV_tx)
SIRV.G.Tx.map$tx_count <- SIRV.G.Tx.map %>% calcTxNum

# Sequins
sequins_tx  <- GenomicFeatures::makeTxDbFromGFF(params$sequins_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
sequins.G.Tx.map <- tibble(
  gene_id = sequins_tx$gene_id %>% unlist,
  tx_name = sequins_tx$tx_name
)
rm(sequins_tx)
sequins.G.Tx.map$tx_count <- sequins.G.Tx.map %>% calcTxNum

# Human
human_tx  <- GenomicFeatures::makeTxDbFromGFF(params$human_gtf, format = "gtf") %>% GenomicFeatures::transcripts(columns = c("tx_id", "tx_name", "gene_id"))
human.G.Tx.map <- tibble(
  gene_id = human_tx$gene_id %>% unlist,
  tx_name = human_tx$tx_name
)
rm(human_tx)
human.G.Tx.map$tx_count <- human.G.Tx.map %>% calcTxNum

# --------------------------------------------------------------
# OUTPUT
# --------------------------------------------------------------
# Save the DGE objects
saveRDS(
    list(
        SIRV.G.Tx.map = SIRV.G.Tx.map,
        sequins.G.Tx.map = sequins.G.Tx.map,
        human.G.Tx.map = human.G.Tx.map
    ),
    file = output.rds
)