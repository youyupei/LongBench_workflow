# depends: FLAMES, dplyr, ggplot2

# R
library(dplyr)
library(ggplot2)
library(GenomicFeatures)


# take commandline arguments for bamfil and output filename
args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
gtf_file <- args[2]
output_file <- args[3]


plot_coverage <- function(bam, isoform = NULL, length_bins = c(0, 1,
  2, 5, 10, Inf), weight_fn = "read_counts", gtf=NULL) {
  
  transcript_info <- transcript_coverage(bam, isoform, gtf, length_bins, weight_fn)

  if (!is.null(isoform)) {
    p <- transcript_info |>
      tidyr::as_tibble(rownames = "transcript") |>
      tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
      dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
      ggplot2::ggplot(aes(x = x, y = coverage, color = transcript)) + geom_line()
    return(p)
  }

  mean_coverage <- transcript_info |>
    dplyr::group_by(length_bin) |>
    dplyr::summarise(dplyr::across(paste0("coverage_", 1:100), ~stats::weighted.mean(.,
      w = weight)))

  p <- mean_coverage |>
    tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
    dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
    ggplot2::ggplot(aes(x = x, y = coverage, color = length_bin)) + geom_line()

  return(p)
}


transcript_coverage <- function(bam, isoform = NULL, gtf=NULL, length_bins = c(0, 1,
  2, 5, 10, Inf), weight_fn = "read_counts") {
  
  if (!is(bam, "GAlignments")) {
    bam <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
  }

  if (!is.null(isoform)) {
    bam <- bam[GenomicAlignments::seqnames(bam) %in% isoform]
  }

  read_counts <- table(GenomicAlignments::seqnames(bam))
  transcript_names <- names(read_counts)
  transcript_info <- data.frame(tr_length =  GenomeInfoDb::seqlengths(bam)[transcript_names], 
    read_counts = as.data.frame(read_counts[transcript_names])$Freq)
  transcript_info$length_bin <- cut(transcript_info$tr_length/1000, length_bins)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage()


  if (!is.null(gtf)){
    gtf <- GenomicFeatures::makeTxDbFromGFF(gtf)

    print("Getting UTRs:")
    threeUTRs.data <- GenomicFeatures::threeUTRsByTranscript(gtf, use.names=T)
    fiveUTRs.data <- GenomicFeatures::fiveUTRsByTranscript(gtf, use.names=T)


    print("Calculating UTR length:")
    threeUTR.len <- threeUTRs.data  |> width() |> sum()
    fiveUTR.len <- fiveUTRs.data  |> width() |> sum()


    print("Removing transcripts with no UTR information:")
    # remove transcripts with no 3' UTR or 5' UTR information
    shared_names <- intersect(names(threeUTR.len), names(fiveUTR.len))
    shared_names <- intersect(shared_names, names(cover))
    transcript_names <- shared_names
    transcript_info <- transcript_info[transcript_names,]
    threeUTR.len <- threeUTR.len[shared_names]
    fiveUTR.len <- fiveUTR.len[shared_names]
    cover <- cover[shared_names]


    print("Removing UTRs:")
    tx_names  <- names(cover)
    cover <-
      lapply(tx_names, function(x) {
        cover[[x]][fiveUTR.len[x]: (length(cover[[x]])-threeUTR.len[x])]
      })
    names(cover) <- tx_names
  }

  print("Binning coverage:")
  cover <- cover |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    subset(select = transcript_names) |>
    t() |>
    as.data.frame()
  colnames(cover) <- paste0("coverage_", 1:100)
  cover <- cover[transcript_names, ]

  if (weight_fn == "sigmoid") {
    weight_fn <- function(mat, read_counts) {
      sigmoid <- function(x) {
        exp(x)/(exp(x) + 1)
      }
      sigmoid((read_counts - 2000)/500)
    }
  } else if (weight_fn == "read_counts") {
    weight_fn <- function(mat, read_counts) {read_counts}
  }

  cover <- cover/transcript_info$read_counts  # scale by read counts
  transcript_info <- cbind(transcript_info, cover)
  transcript_info$weight <- weight_fn(mat, transcript_info$read_counts)
  transcript_info[is.na(transcript_info)] <- 0
  return(transcript_info)
}

# check wether the output filenames contains "CDS"
if (grepl('CDS', output_file)) {
  p <- plot_coverage(bam_file, length_bins = c(seq(0,4), Inf),weight_fn='read_counts', gtf = gtf_file)
} else {
  p <- plot_coverage(bam_file, length_bins = c(seq(0,4), Inf),weight_fn='read_counts', gtf = NULL)
}
  ggsave(output_file, p, dpi = 300)

