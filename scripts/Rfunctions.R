# Setup meta data if not already defined
bulk.meta <- read.csv("/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt")
rownames(bulk.meta) <- bulk.meta$sample


catchOarfish <- function(paths, verbose = TRUE, prefix=".") {
    NSamples <- length(paths)
    OK <- requireNamespace("jsonlite", quietly = TRUE)
    if (!OK) 
        stop("jsonlite package required but is not installed (or can't be loaded)")
    OK <- requireNamespace("readr", quietly = TRUE)
    if (!OK) 
        stop("readr package required but is not installed (or can't be loaded)")
    Arrow <- requireNamespace("arrow", quietly = TRUE)
    if (!Arrow) 
        stop("arrow package required but is not installed (or can't be loaded)")
    ResampleType <- rep_len("bootstrap", NSamples)
    for (j in 1L:NSamples) {
        if (verbose) 
            cat("Reading ", paths[j], ", ", sep = "")
        MetaFile <- file.path(paths[j],  paste0(prefix, "meta_info.json"))
        QuantFile <- file.path(paths[j], paste0(prefix,"quant"))
        BootFile <- file.path(paths[j], paste0(prefix,"infreps.pq"))
        if (!file.exists(QuantFile)) 
            stop("quant file not found at specified path")
        Meta <- jsonlite::fromJSON(MetaFile)
        NBoot <- Meta$num_bootstraps
        if (is.null(NBoot)) 
            stop("Can't find number of bootstraps")
        if (verbose) 
            cat(NBoot, "bootstrap", "samples\n")
        if (j == 1L) {
            Quant1 <- suppressWarnings(readr::read_tsv(QuantFile, 
                col_types = "cdd", progress = FALSE))
            NTx <- nrow(Quant1)
            Counts <- matrix(0, NTx, NSamples)
            DF <- rep_len(0L, NTx)
            OverDisp <- rep_len(0, NTx)
            Counts[, 1L] <- Quant1$num_reads
        }
        else {
            Quant <- suppressWarnings(readr::read_tsv(QuantFile, 
                col_types = "__d", progress = FALSE))
            Counts[, j] <- Quant$num_reads
        }
        if (NBoot > 0L) {
            Boot <- arrow::read_parquet(
                            BootFile,
                            col_select = NULL,
                            as_data_frame = TRUE
                            )
            M <- rowMeans(Boot)
            i <- (M > 0)
            OverDisp[i] <- OverDisp[i] + rowSums((Boot[i, ] - 
                M[i])^2)/M[i]
            DF[i] <- DF[i] + NBoot - 1L
        }
    }
    i <- (DF > 0L)
    if (sum(i) > 0L) {
        OverDisp[i] <- OverDisp[i]/DF[i]
        DFMedian <- median(DF[i])
        DFPrior <- 3
        OverDispPrior <- median(OverDisp[i])/qf(0.5, df1 = DFMedian, 
            df2 = DFPrior)
        if (OverDispPrior < 1) 
            OverDispPrior <- 1
        OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i])/(DFPrior + 
            DF[i])
        OverDisp <- pmax(OverDisp, 1)
        OverDisp[!i] <- OverDispPrior
    }
    else {
        OverDisp[] <- NA_real_
        OverDispPrior <- NA_real_
    }
    Quant1 <- as.data.frame(Quant1, stringsAsFactors = FALSE)
    dimnames(Counts) <- list(Quant1$tname, paths)
    row.names(Quant1) <- Quant1$tname
    Quant1$tname <- NULL
    Quant1$num_reads <- NULL
    Quant1$Length <- Quant1$len
    Quant1$len <- NULL
    Quant1$Overdispersion <- OverDisp
    list(counts = Counts, annotation = Quant1, overdispersion.prior = OverDispPrior, 
        resample.type = ResampleType)
}



catchMyKallisto <- function(paths) {
    NSamples <- length(paths)
    # read bus_quant_tcc.tsv from each sample
    for (j in 1L:NSamples) {
        QuantFile <- file.path(paths[j], "bus_quant_tcc.tsv")
        if (!file.exists(QuantFile)) {
            print(QuantFile)
            stop("quant file not found at specified path")
        }
        if (j == 1L) {
            Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,
                col_types = "cd", progress = FALSE
            ))
            NTx <- nrow(Quant1)
            Counts <- matrix(0, NTx, NSamples)
            DF <- rep_len(0L, NTx)
            OverDisp <- rep_len(0, NTx)
            Counts[, 1L] <- Quant1$bus_counts
            Quant1 <- as.data.frame(Quant1, stringsAsFactors = FALSE)
            # read length from file.path(paths[j], "transcript_lengths.txt")
            Quant1$Length <- read.table(file.path(paths[j], "transcript_lengths.txt"))$V2
        } else {
            Quant <- suppressWarnings(readr::read_tsv(QuantFile,
                col_types = "_d", progress = FALSE
            ))
            Counts[, j] <- Quant$bus_counts
        }
    }
   
    dimnames(Counts) <- list(Quant1$transcript_id, paths)
    row.names(Quant1) <- Quant1$transcript_id
    Quant1$transcript_id <- NULL
    
    list(counts = Counts, annotation = Quant1)
    
}

# Load dge from output directory of salmon, oarfish or kallisto

get_dge_from_salmon <- function(dir, sample_prefix_regex) {
  #To Fix: at the moment, this function requires the bulk.meta to be defined in the global environment
  oarfish.dirs <- file.path(dir, list.dirs(dir, full.names = FALSE, recursive = FALSE))
  catched.salmon <- edgeR::catchSalmon(oarfish.dirs, verbose = TRUE)
  rst.dge <- DGEList(counts=catched.salmon$counts, gene=catched.salmon$annotation)
  rownames(rst.dge$samples) <- rownames(rst.dge$samples) %>% sub(sample_prefix_regex, '', .)
  colnames(rst.dge$counts) <- colnames(rst.dge$counts) %>% sub(sample_prefix_regex, '', .)
  rst.dge$samples <- merge(rst.dge$samples %>% select(-group), bulk.meta, by = "row.names", all.x = TRUE) %>% select(-Row.names)
  rst.dge$samples <- rst.dge$samples[match(colnames(rst.dge$counts), rst.dge$samples$sample),]
  return(rst.dge)
}

get_dge_from_oarfish <- function(dir, sample_prefix_regex) {
  #To Fix: at the moment, this function requires the bulk.meta to be defined in the global environment
  oarfish.dirs <- file.path(dir, list.dirs(dir, full.names = FALSE, recursive = FALSE))
  catched.oarfish <- catchOarfish(oarfish.dirs, verbose = TRUE)
  rst.dge <- DGEList(counts=catched.oarfish$counts, gene=catched.oarfish$annotation)
  rownames(rst.dge$samples) <- rownames(rst.dge$samples) %>% sub(sample_prefix_regex, '', .)
  colnames(rst.dge$counts) <- colnames(rst.dge$counts) %>% sub(sample_prefix_regex, '', .)
  rst.dge$samples <- merge(rst.dge$samples %>% select(-group), bulk.meta, by = "row.names", all.x = TRUE) %>% select(-Row.names)
  rst.dge$samples <- rst.dge$samples[match(colnames(rst.dge$counts), rst.dge$samples$sample),]
  return(rst.dge)
}

get_dge_from_kallisto <- function(dir, sample_prefix_regex) {
   #To Fix: at the moment, this function requires the bulk.meta to be defined in the global environment
  dirs <- file.path(dir, list.dirs(dir, full.names = FALSE, recursive = FALSE))
  catched.kallisto <- catchMyKallisto(dirs)
  catched.kallisto
  rst.dge <- DGEList(counts=catched.kallisto$counts, gene=catched.kallisto$annotation)
  rownames(rst.dge$samples) <- rownames(rst.dge$samples) %>% sub(sample_prefix_regex, '', .)
  colnames(rst.dge$counts) <- colnames(rst.dge$counts) %>% sub(sample_prefix_regex, '', .)
  rst.dge$samples <- merge(rst.dge$samples %>% select(-group), bulk.meta, by = "row.names", all.x = TRUE) %>% select(-Row.names)
  rst.dge$samples <- rst.dge$samples[match(colnames(rst.dge$counts), rst.dge$samples$sample),]
  return(rst.dge)
}
# Function to import DGE from tximport
get_dge_from_txi <- function(quant_dir, sample_prefix_regex, type, dropInfReps = FALSE) {
  #To Fix: at the moment, this function requires the bulk.meta to be defined in the global environment
  txi <- tximport(quant_dir,
                  type = type,
                  tx2gene = combined_tx2gene,
                  ignoreAfterBar = TRUE,
                  countsFromAbundance = "no",
                  dropInfReps=dropInfReps) # raw counts
  counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)
  sample_names <- gsub(sample_prefix_regex, '', basename(dirname(quant_dir)))
  colnames(counts) <- sample_names
  rst.dge <- DGEList(counts = counts)
  rst.dge$genes <- data.frame(Length=txi$length %>% rowMeans())
  rownames(rst.dge$samples) <- sub(sample_prefix_regex, '', rownames(rst.dge$samples))
  colnames(rst.dge$counts) <- sub(sample_prefix_regex, "", colnames(rst.dge$counts))
  
  rst.dge$samples <- merge(rst.dge$samples %>% select(-group), bulk.meta, by = "row.names", all.x = TRUE) %>% select(-Row.names)
  rst.dge$samples <- rst.dge$samples[match(colnames(rst.dge$counts), rst.dge$samples$sample),]
  return(rst.dge)
}

