## cyCombine correction, non-interactive.
##
## Route: cyCombine::normalize -> cyCombine::create_som -> cyCombine::correct_data.
## The ComBat route is pinned by the Python adapter; this script refuses to run if the
## caller asks for anything else, because cyCombine's alternative (ComBat_seq) is a
## count model and transformed cytometry expression is not counts.

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

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  if (!identical(as.character(p$method), "ComBat")) {
    stop(sprintf(
      "refusing method '%s': this adapter is pinned to the continuous-expression ComBat route",
      as.character(p$method)))
  }
  if (!requireNamespace("cyCombine", quietly = TRUE)) {
    stop("the cyCombine R package is not installed in this environment")
  }

  input <- read_binmat(p$input_matrix)
  markers <- as.character(p$markers)
  mat <- input$matrix
  if (!identical(colnames(mat), markers)) {
    stop("column names of the input matrix do not match the declared marker list")
  }

  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df$batch <- as.character(p$batch)
  df$id <- seq_len(nrow(df))
  covar_name <- as_opt(p$covar_name)
  if (!is.null(covar_name)) {
    df[[covar_name]] <- as.character(p$covar)
  }

  warnings_out <- character(0)

  set.seed(as.integer(p$seed))
  normalised <- cyCombine::normalize(df,
                                     markers = markers,
                                     norm_method = as.character(p$norm_method))
  labels <- cyCombine::create_som(normalised,
                                  markers = markers,
                                  seed = as.integer(p$seed),
                                  cluster_method = as.character(p$cluster_method),
                                  distf = as.character(p$distf),
                                  rlen = as.integer(p$rlen),
                                  xdim = as.integer(p$xdim),
                                  ydim = as.integer(p$ydim),
                                  nClus = as_opt(p$nClus))
  labels <- as.character(labels)

  ## Identifiability inside the correction partitions, not only globally: ComBat is
  ## fitted per SOM node, so a node holding one batch contributes no adjustment.
  node_batches <- tapply(df$batch, labels, function(x) length(unique(x)))
  confounded_nodes <- names(node_batches)[node_batches < 2]
  if (!is.null(covar_name)) {
    node_confound <- tapply(seq_along(labels), labels, function(idx) {
      tab <- table(df$batch[idx], df[[covar_name]][idx])
      all(rowSums(tab > 0) <= 1)
    })
    confounded_nodes <- union(confounded_nodes,
                              names(node_confound)[vapply(node_confound, isTRUE, logical(1))])
  }

  corrected <- cyCombine::correct_data(
    df,
    label = labels,
    markers = markers,
    method = "ComBat",
    covar = covar_name,
    anchor = as_opt(p$anchor),
    ref.batch = as_opt(p$ref_batch),
    parametric = isTRUE(p$parametric)
  )

  corrected <- corrected[order(corrected$id), , drop = FALSE]
  if (!identical(as.integer(corrected$id), as.integer(df$id))) {
    stop("cyCombine did not return one row per input event in a recoverable order")
  }
  out <- as.matrix(corrected[, markers, drop = FALSE])
  if (any(!is.finite(out))) {
    stop("cyCombine returned non-finite values; they are not silently repaired")
  }

  write_binmat(p$output_matrix, out, rownames(mat), markers)

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    n_cells = nrow(out),
    n_markers = ncol(out),
    n_som_nodes = length(unique(labels)),
    confounded_nodes = as.character(confounded_nodes),
    cycombine_version = as.character(utils::packageVersion("cyCombine")),
    warnings = warnings_out,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
