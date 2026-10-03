## Shared setup for the R-side test suite.
##
## The R scripts under workflow/scripts/ are hand-mirrored against their Python
## counterparts — `lib/binmat.R` against `src/batchdx/cyto/io/binary.py` is the
## explicit case, and it says so in its own header. A mirror is only as good as the
## tests that hold the two in agreement, and a divergence here is the quiet kind: the
## matrix values stay plausible while the row identity or the memory order is wrong.
##
## Nothing in this file needs a package beyond `jsonlite`, which `binmat.R` itself
## already requires. `testthat` is the only addition, and it is a test-time dependency
## only. The cross-language tests additionally need a Python that can import
## `batchdx.cyto`; those skip with a named reason when there is not one, and the runner
## prints the skip count so a skip is visible rather than silent.

suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(jsonlite))

batchdx_test_root <- function() {
  ## Walk up from the working directory until both markers of the repository root are
  ## present. Hard-coding a relative path here would break the moment the suite is run
  ## from anywhere other than the root, which is exactly what a CI runner does.
  dir <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "pyproject.toml")) &&
        file.exists(file.path(dir, "workflow", "scripts", "lib", "binmat.R"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("could not locate the batchdx-cyto repository root from ", getwd())
    }
    dir <- parent
  }
}

BATCHDX_ROOT <- batchdx_test_root()
source(file.path(BATCHDX_ROOT, "workflow", "scripts", "lib", "binmat.R"))

## A fresh directory for one test. Under `withr`-less testthat 3 the idiomatic form is a
## local variable plus `withr::defer`, but withr is not a dependency here, so tests take
## `tempdir()`-rooted paths and let the session clean up.
batchdx_workdir <- function(name) {
  path <- file.path(tempdir(), paste0("binmat-", name, "-", as.integer(runif(1, 1e6, 1e7))))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

## The A and B of a binary matrix, computed without touching `write_binmat`.
##
## This is the point of the whole exercise. If the expected bytes are produced by the
## function under test, a test can only prove the function agrees with itself; a
## row-major/column-major or byte-order divergence between R and Python is a
## disagreement between two implementations, and the reference has to be a third,
## independent thing. Base R `writeBin` plus an explicitly transposed vector is that
## third thing: it is the definition of the layout, written out longhand.
batchdx_expected_row_major_bytes <- function(m) {
  ## `writeBin` against a `raw()` connection returns the encoded bytes rather than
  ## writing anywhere, which is exactly what is wanted: no temporary file, no connection
  ## to leak, and the encoding is done by base R's own writer rather than by the
  ## package's.
  writeBin(as.vector(t(as.matrix(m))), raw(), size = 8, endian = "little")
}

## Write raw bytes to a path without leaving a connection open on the error path.
batchdx_write_raw <- function(path, bytes) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(bytes, con)
  invisible(path)
}

## A Python list literal, built by hand. `deparse` would emit R's `c(...)` with a
## trailing comma, which is not valid Python.
batchdx_python_list <- function(values) {
  sprintf("[\"%s\"]", paste(as.character(values), collapse = "\", \""))
}

## A sidecar written by hand, so `read_binmat` is never fed a file the code under test
## produced.
batchdx_write_sidecar <- function(path, n_rows, n_cols, row_ids, col_names,
                                  byte_order = "little", layout = "row_major",
                                  dtype = "float64", magic = BATCHDX_BINMAT_MAGIC) {
  writeLines(
    jsonlite::toJSON(
      list(magic = magic, dtype = dtype, byte_order = byte_order, layout = layout,
           n_rows = n_rows, n_cols = n_cols,
           row_ids = as.character(row_ids), col_names = as.character(col_names)),
      auto_unbox = TRUE, pretty = TRUE
    ),
    paste0(path, ".json")
  )
}

## A matrix whose values cannot be confused under a transpose: every entry distinct.
batchdx_distinct_matrix <- function(n_rows, n_cols) {
  m <- matrix(seq_len(n_rows * n_cols) + 0.5, nrow = n_rows, ncol = n_cols, byrow = TRUE)
  m
}

## Locate a Python that can import batchdx.cyto, or NULL. `BATCHDX_PYTHON` wins so a
## runner can point at a specific interpreter; the repository's own `.venv` is next,
## because that is the environment the pytest suite uses.
batchdx_python <- function() {
  candidates <- c(Sys.getenv("BATCHDX_PYTHON"), file.path(BATCHDX_ROOT, ".venv", "bin", "python"),
                  Sys.getenv("BATCHDX_PYTHON_FALLBACK", "python3"))
  for (candidate in candidates) {
    if (!nzchar(candidate)) next
    ok <- suppressWarnings(system2(
      candidate, c("-c", shQuote("import batchdx.cyto")),
      stdout = FALSE, stderr = FALSE
    ))
    if (identical(as.integer(ok), 0L)) return(candidate)
  }
  NULL
}
