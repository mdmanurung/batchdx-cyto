## CytoBatchFlagR, non-interactive.
##
## The upstream sources are GPL-3 and are distributed as scripts rather than as a
## package.  They are source()d here UNMODIFIED from a checkout the caller supplies;
## nothing is copied into batchdx-cyto and nothing is refactored.  The upstream
## workflow script is deliberately not sourced, because it calls
## rstudioapi::selectDirectory and Shiny helpers; the analysis functions in
## R/CytoBatchFlagR_functions.R are called directly with configured values instead.
##
## Everything in this file is original code.  It may read the upstream file, it may
## check that this machine can satisfy it, and it may pass upstream's own documented
## arguments different values; it does not alter a line of upstream source.  That is
## the distinction docs/decisions/0003-cytobatchflagr-integration.md draws, and it is
## why a trimmed copy of the upstream file is not an available shortcut: trimming it
## would be the derivative work the decision forbids.

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

## A JSON object rather than a JSON array.  jsonlite serialises a plain empty list() as
## "[]", and the adapter reads these fields as mappings, so an absent map has to be sent
## as an empty object or the reader trips over it.
empty_map <- function() stats::setNames(list(), character(0))

## The packages the upstream analysis file needs, read out of the file rather than
## transcribed here.  Upstream has no DESCRIPTION and no NAMESPACE -- the `library()`
## calls are the only machine-readable statement of its dependencies -- so a list
## written down here would be a second source of truth that can drift from the
## checkout actually in use.  Reading the file makes the preflight self-updating.
upstream_required_packages <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hits <- regmatches(lines, regexec("library\\(\\s*([A-Za-z][A-Za-z0-9._]*)", lines))
  pkgs <- unlist(lapply(hits, function(m) if (length(m) > 1L) m[[2L]] else NULL))
  sort(unique(pkgs))
}

## Cells on each side of the threshold, per control sample and batch, per marker.
##
## Two of the three upstream summarisers, neg_mfi_df_all and pos_mfi_df_all, default to
## min_cells = TRUE and then silently drop any control with fewer than 100 cells on that
## side of a marker threshold, printing the fact to stdout and nothing else
## (functions.R:1041, 1157, 1084, 1090-1102).  The third, percent_pos_df_all, has no
## min_cells argument at all (functions.R:1272-1276) and drops a control the same way by
## emitting rows only for sample_ids holding at least one cell above the threshold.
##
## A dropped control is a missing row, and a missing row is precisely the failure this
## adapter must not commit: with fewer than four values the 1.5 IQR fence cannot fire, so
## a silent drop converts a detectable batch effect into an apparent absence of one.  The
## counts are therefore computed here and returned, and the caller's coverage check turns
## a short panel into a reported ineligibility.
sign_counts <- function(df, markers, cutv) {
  ## The unit here is the control sample AS STAINED IN A GIVEN BATCH, not the control
  ## identity.  One reference aliquot measured in four batches is deliberately the same
  ## value in the `control` column for all four -- that equality is what iqr_metric_function
  ## groups on, so the four per-batch values are the ones the fence compares.  Keying on
  ## `control` alone would merge the four batches into a single value and leave nothing
  ## for the fence to compare, which is the failure this whole path exists to avoid.
  sid_all <- paste(as.character(df$batch), as.character(df$control), sep = "/")
  rows <- list()
  for (m_i in seq_along(markers)) {
    mn <- markers[m_i]
    cutoff <- cutv[m_i]
    expr <- suppressWarnings(as.numeric(df[[mn]]))
    for (sid in unique(sid_all)) {
      sel <- sid_all == sid
      e <- expr[sel]
      idx <- which(sel)[1L]
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sid,
        control   = as.character(df$control[idx]),
        batch     = as.character(df$batch[idx]),
        marker    = mn,
        n_negative = sum(e < cutoff, na.rm = TRUE),
        n_positive = sum(e >= cutoff, na.rm = TRUE),
        n_cells    = sum(is.finite(e)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  functions_file <- file.path(p$upstream_dir, "R", "CytoBatchFlagR_functions.R")
  if (!file.exists(functions_file)) {
    stop(sprintf("CytoBatchFlagR sources not found at %s", functions_file))
  }

  ## Preflight before source().  Every one of these library() calls is unconditional,
  ## so without this check a single absent package aborts the whole source() and the
  ## run reaches Python as a non-zero exit, which the adapter reports as
  ## status="failed" with a raw R error in the log.  A missing package is not a failed
  ## detection: it is a stated unavailability, and it is reported as one.
  required <- upstream_required_packages(functions_file)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    batchdx_write_result(parsed$result_path, list(
      status = "not_available",
      reason = sprintf(
        paste0("CytoBatchFlagR cannot be sourced: %d of its %d unconditionally ",
               "library() packages are not installed in this R library: %s. Upstream has ",
               "no DESCRIPTION and no NAMESPACE, so the library() block in ",
               "R/CytoBatchFlagR_functions.R is the only statement of its dependencies, ",
               "and one absent package aborts the entire source() rather than degrading ",
               "one function. Build the environment that provides them with ",
               "scripts/setup_r_env.sh. No flag was produced and none is implied: an ",
               "absent dependency is never reported as an absence of batch effects."),
        length(missing), length(required), paste(missing, collapse = ", ")),
      missing_packages = as.list(missing),
      required_packages = as.list(required),
      tables = empty_map(), plots = empty_map(),
      commit = upstream_commit(p$upstream_dir),
      warnings = character(0),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
    quit(save = "no", status = 0L)
  }

  source(functions_file, local = FALSE)

  input <- read_binmat(p$input_matrix)
  markers <- as.character(p$markers)
  df <- as.data.frame(input$matrix, stringsAsFactors = FALSE)
  df$control <- as.character(p$control)
  df$batch <- as.character(p$batch)
  ## No `sample_id` column is supplied.  Upstream builds one itself, per helper, as
  ## paste0(batch, "_", control) in neg_mfi_df_all/pos_mfi_df_all and
  ## paste0(control, "_", batch) in percent_pos_df_all (functions.R:1049, 1165, 1283).
  ## Both are per-batch identities, and letting upstream build them keeps the per-batch
  ## split intact; supplying a sample_id that is constant across batches merges all four
  ## batches into one value and leaves the fence nothing to compare.

  ## The upstream helpers take the thresholds as a two-column frame aligned with
  ## `markers`; they are fixed in the declared coordinate system before this call and
  ## are never retuned per candidate.
  cutf <- data.frame(marker = markers, threshold = as.numeric(p$thresholds),
                     stringsAsFactors = FALSE)

  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)

  counts <- sign_counts(df, markers, as.numeric(p$thresholds))

  ## min_cells = FALSE where upstream offers it, and the reason is not performance: with
  ## it left at its default any control with fewer than 100 cells on that side of a
  ## threshold vanishes from the output, and a vanished control reduces the per-marker
  ## value count below the four the Tukey fence needs, so the batch is silently reported
  ## as unflagged.  percent_pos_df_all offers no such argument (functions.R:1272-1276) and
  ## drops a control the same way regardless.  The cell counts computed above are returned
  ## in the payload as `coverage` and judged by the caller, which is what makes a shortfall
  ## reported rather than silent.
  neg <- neg_mfi_df_all(df = df, markers = markers, cutf = cutf,
                        sample_col_name = "control", batch_col_name = "batch",
                        min_cells = FALSE)
  pos <- pos_mfi_df_all(df = df, markers = markers, cutf = cutf,
                        sample_col_name = "control", batch_col_name = "batch",
                        min_cells = FALSE)
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

  counts_csv <- file.path(p$out_dir, "cytobatchflagr_cell_counts.csv")
  utils::write.csv(counts, counts_csv, row.names = FALSE)

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    iqr = rows,
    coverage = counts,
    tables = list(iqr = csv, cell_counts = counts_csv),
    plots = empty_map(),
    commit = upstream_commit(p$upstream_dir),
    required_packages = as.list(required),
    warnings = character(0),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
