## Typed binary matrix exchange, R side.
##
## Mirrors batchdx.cyto.io.binary.  Uses base R only, so adapters add no R package
## dependency beyond what the upstream tool itself needs (plus jsonlite for the
## sidecar and the result payload).

BATCHDX_BINMAT_MAGIC <- "batchdx-binmat-1"

read_binmat <- function(path) {
  sidecar <- paste0(path, ".json")
  if (!file.exists(sidecar)) {
    stop(sprintf("sidecar not found for %s: %s", path, sidecar))
  }
  meta <- jsonlite::fromJSON(sidecar, simplifyVector = TRUE)
  if (!identical(meta$magic, BATCHDX_BINMAT_MAGIC)) {
    stop(sprintf("unexpected magic %s in %s", meta$magic, sidecar))
  }
  n <- as.integer(meta$n_rows) * as.integer(meta$n_cols)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  values <- readBin(con, what = "double", n = n, size = 8, endian = "little")
  if (length(values) != n) {
    stop(sprintf("expected %d values in %s, read %d", n, path, length(values)))
  }
  m <- matrix(values,
              nrow = as.integer(meta$n_rows),
              ncol = as.integer(meta$n_cols),
              byrow = TRUE)
  colnames(m) <- as.character(meta$col_names)
  rownames(m) <- as.character(meta$row_ids)
  list(matrix = m, meta = meta)
}

write_binmat <- function(path, m, row_ids, col_names) {
  stopifnot(nrow(m) == length(row_ids), ncol(m) == length(col_names))
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.vector(t(as.matrix(m))), con, size = 8, endian = "little")
  meta <- list(
    magic = BATCHDX_BINMAT_MAGIC,
    dtype = "float64",
    byte_order = "little",
    layout = "row_major",
    n_rows = nrow(m),
    n_cols = ncol(m),
    row_ids = I(as.character(row_ids)),
    col_names = I(as.character(col_names))
  )
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
             paste0(path, ".json"))
  invisible(path)
}

## Read the params file passed by batchdx.cyto.adapters.r_runner.run_r_script.
batchdx_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("usage: Rscript <script>.R <params.json> <result.json>")
  }
  list(
    params = jsonlite::fromJSON(args[[1]], simplifyVector = TRUE),
    result_path = args[[2]]
  )
}

batchdx_write_result <- function(path, payload) {
  ## auto_unbox must not turn one warning or node ID into a scalar string, and a
  ## NAME must not turn the array into a JSON object: c() and unlist(lapply(...))
  ## propagate the names of the lists they combine (a sample_id, a node id), and
  ## jsonlite serialises a named list as an object, which no reader of these fields
  ## accepts. The validator in r_runner.py requires an array of strings, so the
  ## names are stripped here, once, for every script that writes a result.
  list_fields <- c("warnings", "confounded_nodes", "fitted_labels", "unfitted_labels",
                   "missing_packages", "metaclusters", "outputs")
  for (field in intersect(names(payload), list_fields)) {
    payload[[field]] <- unname(as.list(payload[[field]]))
  }
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", pretty = TRUE), path)
  invisible(path)
}

## Normalise jsonlite's representation of an absent value to R's NULL.
as_opt <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 1 && is.na(x)) return(NULL)
  x
}
