## cyCombine native metric evaluation, non-interactive.
##
## Route: cyCombine::evaluate_emd and cyCombine::evaluate_mad.  These are the two
## functions registered in the `native` namespace of
## batchdx.cyto.metrics.registry, and they are the only place this package uses an
## upstream *metric* rather than an upstream *correction*.  Both take an uncorrected
## and a corrected dataframe and reduce the marker-wise distances to a single relative
## reduction, so a result is only interpretable against the pair it was computed from.
##
## Structure of both inputs is fixed by the pinned API at fcd42dbb:
##   evaluate_emd(uncorrected, corrected, binSize, cell_col, batch_col, markers, ...)
##   evaluate_mad(uncorrected, corrected,             cell_col, batch_col, markers, ...)
## and both call check_colname() on cell_col, batch_col AND 'id'.  So the `id` column
## is a hard requirement, not a convenience: it is what the two frames are aligned on.

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

## Build one of the two dataframes evaluate_* expects.
##
## `id` is the event identity carried across the boundary in the binary matrix's
## row_ids, so the uncorrected and corrected frames are aligned by event and not by
## position.  That matters: a candidate correction preserves row order today, but
## nothing in the contract requires it, and silently comparing different cells under
## the same label would produce a plausible-looking reduction with no meaning.
as_eval_frame <- function(path, batch, sample, label, cell_col) {
  input <- read_binmat(path)
  df <- as.data.frame(input$matrix, stringsAsFactors = FALSE)
  if (nrow(df) != length(batch)) {
    stop(sprintf("frame from %s has %d rows but %d batch labels were supplied",
                 path, nrow(df), length(batch)))
  }
  df$id <- as.character(input$meta$row_ids)
  df$batch <- as.character(batch)
  # evaluate_* only needs `sample` if a caller asks for it; it is carried because
  # compute_emd (reached by evaluate_emd) groups by a column literally named "sample".
  df$sample <- as.character(sample)
  if (is.null(label)) {
    df[[cell_col]] <- "all"
  } else {
    if (length(label) != nrow(df)) {
      stop(sprintf("label vector has length %d but %s has %d rows",
                   length(label), path, nrow(df)))
    }
    df[[cell_col]] <- as.character(label)
  }
  df
}

## One row of the per-marker before/after table evaluate_emd returns, computed from
## the returned tibble rather than scraped from a plot or a log.
emd_to_rows <- function(emd) {
  if (is.null(emd) || is.null(emd$emd)) return(list())
  tab <- as.data.frame(emd$emd)
  rows <- list()
  for (i in seq_len(nrow(tab))) {
    rows[[length(rows) + 1L]] <- list(
      name = as.character(tab$Name[[i]]),
      corrected = as.numeric(tab$Corrected[[i]]),
      uncorrected = as.numeric(tab$Uncorrected[[i]]),
      reduction = as.numeric(tab$Reduction[[i]])
    )
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

  markers <- as.character(p$markers)
  cell_col <- as_opt(p$cell_col)
  if (is.null(cell_col)) cell_col <- "label"

  uncorrected <- as_eval_frame(
    p$uncorrected_matrix, p$batch, p$sample, as_opt(p$label), cell_col
  )
  corrected <- as_eval_frame(
    p$corrected_matrix, p$batch, p$sample, as_opt(p$label), cell_col
  )

  ## Event identity must line up before anything is compared.
  if (!identical(uncorrected$id, corrected$id)) {
    stop("the uncorrected and corrected frames do not carry the same event ids in the same order")
  }

  emd <- cyCombine::evaluate_emd(
    uncorrected = uncorrected, corrected = corrected,
    binSize = as.numeric(p$bin_size), cell_col = cell_col, batch_col = "batch",
    markers = markers, plots = FALSE
  )
  mad <- cyCombine::evaluate_mad(
    uncorrected = uncorrected, corrected = corrected,
    cell_col = cell_col, batch_col = "batch", markers = markers
  )

  payload <- list(
    status = "ok",
    emd_reduction = as.numeric(emd$reduction),
    # evaluate_mad returns list(score =, mad =), not a "reduction" field: its summary
    # number is the median absolute difference between the uncorrected and corrected
    # per-(population, batch, marker) MADs.  Reading a "reduction" key here would yield
    # NULL and report the metric as undefined for every input, including real ones.
    mad_score = if (is.null(mad$score)) NA_real_ else as.numeric(mad$score),
    emd_rows = emd_to_rows(emd),
    software = list(
      name = "cyCombine",
      version = as.character(utils::packageVersion("cyCombine"))
    ),
    elapsed_s = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  batchdx_write_result(parsed$result_path, payload)
}

main()
