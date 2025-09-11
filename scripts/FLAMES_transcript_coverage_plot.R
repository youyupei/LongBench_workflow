# Note: FLAMES v2.3.4 or higher is required for this script

list(
  tar_target(lr_bulk_bam,
    list.files(
      "/vast/projects/LongBench/analysis/lr_bulk/result/TranscriptAlignment",
      pattern = "\\.bam$", full.names = TRUE
    ) %>%
      setNames(., basename(.))
  ),
  tar_target(lr_bulk_coverage,
    lr_bulk_bam |>
      BiocParallel::bplapply(
        FLAMES::get_coverage,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = 12, stop.on.error = TRUE
        )
      ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_12c128g")
    )
  ),
  tar_target(lr_bulk_coverage_filtered,
    lr_bulk_coverage |>
      BiocParallel::bplapply(
        FLAMES::filter_coverage,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = 12, stop.on.error = TRUE
        )
      ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "slurm_12c128g")
    )
  ),
  tar_target(lr_bulk_coverage_spikins,
    lr_bulk_coverage |>
      by_protocol() |>
      lapply(\(x) {
        bind_rows(x, .id = "sample") |>
          group_by(transcript) |>
          summarise(dplyr::across(
            paste0("coverage_", 1:100),
            # do a weighted mean of the coverage values, weighted by the read_counts column
            # ~ stats::weighted.mean(., w = read_counts, na.rm = TRUE)
            ~ mean(., na.rm = TRUE) # simple mean of the coverage values
          ))
      }) |>
      bind_rows(.id = "protocol") |>
      mutate(spikin_type = case_when(
        str_detect(transcript, "^R") ~ "Sequins",
        str_detect(transcript, "^ENS") ~ "Human",
        str_detect(transcript, "^SIRV\\d{3}$") ~ "SIRV_E0",
        str_detect(transcript, "^SIRV\\d{4}") ~ "SIRV_Long",
        str_detect(transcript, "^SIRV") ~ "SIRV",
        TRUE ~ "ERCC"
      ))
  ),
  tar_target(lr_bulk_coverage_spikins_plot,
    format = "file",
    {
      out_file <- file.path("output", "lr_bulk_coverage_spikins_plot.pdf")
      p <- lr_bulk_coverage_spikins |>
        dplyr::filter(!spikin_type == "Human") |>
        # simple mean for each spikin
        group_by(spikin_type, protocol) |>
        summarise(dplyr::across(
          paste0("coverage_", 1:100),
          ~ mean(., na.rm = TRUE)
        )) |>
        ungroup() |>
        pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
        mutate(x = as.numeric(str_remove(x, "coverage_"))) |>
        ggplot(aes(x = x, y = coverage, color = spikin_type)) +
        geom_line() +
        facet_wrap(~protocol) +
        theme_minimal() +
        labs(
          title = "Coverage of spike-ins",
          x = "% Position in transcript",
          y = "Coverage"
        ) +
        guides(color = guide_legend(title = "Spikin type"))
      ggsave(out_file, p, width = 16, height = 9)
      out_file
    }
  ),

  tar_target(
    lr_bulk_coverage_filtered_averaged,
    lr_bulk_coverage_filtered |>
      by_protocol() |>
      lapply(\(x) {
        bind_rows(x, .id = "sample") |>
          group_by(transcript) |>
          summarise(
            # do a weighted mean of the coverage values, weighted by the read_counts column
            # ~ stats::weighted.mean(., w = read_counts, na.rm = TRUE)
            dplyr::across(
              paste0("coverage_", 1:100),
              ~ mean(., na.rm = TRUE) # simple mean of the coverage values
            ),
            read_counts = sum(read_counts, na.rm = TRUE),
            tr_length = dplyr::first(tr_length)
          )
      }) |>
      bind_rows(.id = "protocol"),
    deployment = "main"
  ),

  tar_target(
    lr_bulk_coverage_filtered_weighted,
    lr_bulk_coverage_filtered |>
      by_protocol() |>
      lapply(\(x) {
        bind_rows(x, .id = "sample") |>
          group_by(transcript) |>
          summarise(
            dplyr::across(
              paste0("coverage_", 1:100),
              ~ stats::weighted.mean(., w = read_counts, na.rm = TRUE)
            ),
            read_counts = sum(read_counts, na.rm = TRUE),
            tr_length = dplyr::first(tr_length)
          )
      }) |>
      bind_rows(.id = "protocol"),
    deployment = "main"
  ),

  tar_target(lr_bulk_coverage_plot,
    format = "file",
    {
      out_file <- file.path("output", "lr_bulk_coverage_plot.pdf")

      lr_bulk_coverage_filtered_weighted |>
        mutate(
          length_bin = cut(
            tr_length,
            breaks = c(0, 1000, 2500, Inf),
            labels = c("0-1kb", "1-2.5kb", "2.5kb+")
          )
        ) |>
        group_by(length_bin, protocol) |>
        dplyr::summarise(
          dplyr::across(
            paste0("coverage_", 1:100),
            # ~ mean(., na.rm = TRUE)
            ~ stats::weighted.mean(., w = read_counts, na.rm = TRUE)
          ),
          .groups = "drop"
        ) |>
        tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
        dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
        ggplot2::ggplot(aes(x = x, y = coverage, color = protocol)) +
        ggplot2::geom_line() +
        ggplot2::xlab("% Position in transcript") +
        ggplot2::ylab("Weighted average coverage") +
        ylim(0.25, 1) +
        facet_wrap(~length_bin, nrow = 1) +
        ggplot2::theme_minimal()
      ggsave(out_file, width = 16, height = 6)
      out_file
    }
  ),

  tar_target(lr_bulk_average_coverage_plot,
    format = "file",
    {
      out_file <- file.path("output", "lr_bulk_average_coverage_plot.pdf")

      lr_bulk_coverage_filtered |>
        by_protocol() |>
        lapply(\(x) {
          bind_rows(x, .id = "sample") |>
            filter(read_counts >= 10) |>
            mutate(transcript = paste0(transcript, "_", sample))
        }) |>
        bind_rows(.id = "protocol") |>
        mutate(
          length_bin = cut(
            tr_length,
            breaks = c(0, 1000, 2500, Inf),
            labels = c("0-1kb", "1-2.5kb", "2.5kb+")
          )
        ) |>
        group_by(length_bin, protocol) |>
        dplyr::summarise(
          dplyr::across(
            paste0("coverage_", 1:100),
            ~ mean(., na.rm = TRUE)
            # ~ stats::weighted.mean(., w = read_counts, na.rm = TRUE)
          ),
          .groups = "drop"
        ) |>
        tidyr::pivot_longer(paste0("coverage_", 1:100), names_to = "x", values_to = "coverage") |>
        dplyr::mutate(x = as.numeric(gsub("coverage_", "", x))) |>
        ggplot2::ggplot(aes(x = x, y = coverage, color = protocol)) +
        ggplot2::geom_line() +
        ggplot2::xlab("% Position in transcript") +
        ggplot2::ylab("Weighted average coverage") +
        ylim(0.25, 1) +
        facet_wrap(~length_bin, nrow = 1) +
        ggplot2::theme_minimal()
      ggsave(out_file, width = 16, height = 6)
      out_file
    }
  )
)