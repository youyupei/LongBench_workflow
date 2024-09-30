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