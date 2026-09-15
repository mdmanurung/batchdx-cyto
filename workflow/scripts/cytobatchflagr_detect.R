## CytoBatchFlagR, non-interactive.
##
## The upstream sources are GPL-3 and are distributed as scripts rather than as a
## package.  They are source()d here UNMODIFIED from a checkout the caller supplies;
## nothing is copied into batchdx-cyto and nothing is refactored.  The upstream
## workflow script is deliberately not sourced, because it calls
## rstudioapi::selectDirectory and Shiny helpers; the analysis functions in
## R/CytoBatchFlagR_functions.R are called directly with configured values instead.

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

upstream_commit <- function(dir) {
  out <- tryCatch(
    system2("git", c("-C", shQuote(dir), "rev-parse", "HEAD"), stdout = TRUE, stderr = NULL),
    error = function(e) NA_character_)
  if (length(out) == 0 || is.na(out[[1]])) NA_character_ else out[[1]]
}

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  functions_file <- file.path(p$upstream_dir, "R", "CytoBatchFlagR_functions.R")
  if (!file.exists(functions_file)) {
    stop(sprintf("CytoBatchFlagR sources not found at %s", functions_file))
  }
  source(functions_file, local = FALSE)

  for (pkg in c("dplyr", "tidyr", "matrixStats")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("CytoBatchFlagR requires the '%s' package, which is not installed", pkg))
    }
  }

  input <- read_binmat(p$input_matrix)
  markers <- as.character(p$markers)
  df <- as.data.frame(input$matrix, stringsAsFactors = FALSE)
  df$control <- as.character(p$control)
  df$batch <- as.character(p$batch)

  ## The upstream helpers take the thresholds as a two-column frame aligned with
  ## `markers`; they are fixed in the declared coordinate system before this call and
  ## are never retuned per candidate.
  cutf <- data.frame(marker = markers, threshold = as.numeric(p$thresholds),
                     stringsAsFactors = FALSE)

  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)

  neg <- neg_mfi_df_all(df = df, markers = markers, cutf = cutf,
                        sample_col_name = "control", batch_col_name = "batch")
  pos <- pos_mfi_df_all(df = df, markers = markers, cutf = cutf,
                        sample_col_name = "control", batch_col_name = "batch")
  pct <- percent_pos_df_all(df = df, markers = markers, cutf = cutf,
                            sample_col_name = "control", batch_col_name = "batch")

  neg_df <- if (is.list(neg) && !is.data.frame(neg)) neg[[1]] else neg
  pos_df <- if (is.list(pos) && !is.data.frame(pos)) pos[[1]] else pos
  pct_df <- if (is.list(pct) && !is.data.frame(pct)) pct[[1]] else pct

  iqr_df <- iqr_dataframe(neg_df = neg_df, pos_df = pos_df, percent_pos_df = pct_df,
                          markers = markers, sample_col_name = "control",
                          batch_col_name = "batch", output_dir = p$out_dir)

  keep <- intersect(c("type", "marker", "control", "batch", "value", "batch_effect",
                      "flagged_batch"), names(iqr_df))
  rows <- iqr_df[, keep, drop = FALSE]
  csv <- file.path(p$out_dir, "cytobatchflagr_iqr.csv")
  utils::write.csv(rows, csv, row.names = FALSE)

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    iqr = rows,
    tables = list(iqr = csv),
    plots = list(),
    commit = upstream_commit(p$upstream_dir),
    warnings = character(0),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
