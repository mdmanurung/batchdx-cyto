## cyCombine detection, non-interactive.
##
## Native plots come from cyCombine::detect_batch_effect_express / detect_batch_effect.
## Neither returns a structured value, so the tables handed back to Python come from a
## direct call to the exported cyCombine::compute_emd rather than from parsing the
## messages those functions print.

suppressPackageStartupMessages({
  library(jsonlite)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(sub("^--file=", "", file_arg[[1]])))
  "."
}
source(file.path(script_dir(), "lib", "binmat.R"))

emd_to_rows <- function(emd) {
  rows <- list()
  for (cell_type in names(emd)) {
    for (marker in names(emd[[cell_type]])) {
      m <- emd[[cell_type]][[marker]]
      batches <- rownames(m)
      for (i in seq_along(batches)) {
        for (j in seq_along(batches)) {
          value <- m[i, j]
          if (is.na(value)) next
          rows[[length(rows) + 1L]] <- list(
            population = cell_type, marker = marker,
            batch_a = batches[i], batch_b = batches[j], value = as.numeric(value)
          )
        }
      }
    }
  }
  rows
}

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  if (!requireNamespace("cyCombine", quietly = TRUE)) {
    stop("the cyCombine R package is not installed in this environment")
  }

  input <- read_binmat(p$input_matrix)
  markers <- as.character(p$markers)
  df <- as.data.frame(input$matrix, stringsAsFactors = FALSE)
  df$batch <- as.character(p$batch)
  ## detect_batch_effect* groups per-sample medians by a `sample` column; without it
  ## the MDS step fails inside upstream.  The per-cell sample label is required.
  df$sample <- as.character(p$sample)
  df$id <- seq_len(nrow(df))
  labels <- as_opt(p$label)
  if (!is.null(labels)) df$label <- as.character(labels)

  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)
  warnings_out <- character(0)
  plots <- list()

  if (identical(as.character(p$mode), "express")) {
    cyCombine::detect_batch_effect_express(
      df = df, out_dir = p$out_dir, markers = markers,
      batch_col = "batch", seed = as.integer(p$seed), binSize = as.numeric(p$bin_size)
    )
  } else {
    cyCombine::detect_batch_effect(
      df = df, out_dir = p$out_dir, markers = markers,
      batch_col = "batch", label_col = if (is.null(labels)) "label" else "label",
      seed = as.integer(p$seed), binSize = as.numeric(p$bin_size)
    )
  }
  for (f in list.files(p$out_dir, pattern = "[.](png|pdf)$", full.names = TRUE)) {
    plots[[basename(f)]] <- f
  }

  ## Structured tables: computed, not scraped.
  if (is.null(labels)) {
    df$label <- "all"
    warnings_out <- c(warnings_out,
                      "no frozen evaluation labels were supplied; EMD is computed over all cells")
  }
  emd <- cyCombine::compute_emd(df, binSize = as.numeric(p$bin_size),
                                cell_col = "label", batch_col = "batch", markers = markers)
  emd_rows <- emd_to_rows(emd)

  counts <- table(df$batch, df$label)
  perc <- (counts / rowSums(counts)) * 100
  cluster_rows <- list()
  for (b in rownames(perc)) {
    for (cl in colnames(perc)) {
      cluster_rows[[length(cluster_rows) + 1L]] <- list(
        batch = b, population = cl, value = as.numeric(perc[b, cl])
      )
    }
  }

  emd_csv <- file.path(p$out_dir, "cycombine_emd.csv")
  utils::write.csv(do.call(rbind, lapply(emd_rows, as.data.frame)), emd_csv, row.names = FALSE)

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    emd = emd_rows,
    cluster_percentage = cluster_rows,
    plots = plots,
    tables = list(emd = emd_csv),
    warnings = warnings_out,
    cycombine_version = as.character(utils::packageVersion("cyCombine")),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
