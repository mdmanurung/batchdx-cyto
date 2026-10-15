test_that("write_binmat emits row-major little-endian float64, verified without itself", {
  ## The core of the mirror. Python's `io/binary.py` writes
  ## `arr.tobytes(order="C")` — row-major, little-endian, float64, no header. If R
  ## writes column-major, or big-endian, or float32, then every value in every matrix
  ## crossing the boundary is wrong, and nothing about the failure looks like a failure:
  ## a transposed cytometry matrix is still a plausible cytometry matrix.
  ##
  ## So the expected bytes are produced here by hand, from the definition of the layout,
  ## and `write_binmat`'s output is compared against them byte for byte. The function
  ## under test is not allowed to generate its own reference.
  for (shape in list(c(3L, 7L), c(7L, 3L), c(4L, 4L), c(1L, 5L), c(5L, 1L), c(1L, 1L))) {
    m <- batchdx_distinct_matrix(shape[[1]], shape[[2]])
    path <- file.path(batchdx_workdir("layout"), "m.f64")
    write_binmat(path, m, row_ids = paste0("r", seq_len(shape[[1]])),
                 col_names = paste0("c", seq_len(shape[[2]])))

    actual <- readBin(path, what = "raw", n = file.info(path)$size)
    expected <- batchdx_expected_row_major_bytes(m)
    expect_identical(actual, expected,
                     info = sprintf("shape %dx%d is not row-major little-endian float64",
                                    shape[[1]], shape[[2]]))
  }
})

test_that("write_binmat's sidecar declares the layout it actually writes", {
  ## The sidecar is a contract, and a contract that disagrees with the bytes is worse
  ## than no contract: a reader that trusts `byte_order` and is handed the other one
  ## gets silent nonsense. Python's `read_matrix` does not check either field — it
  ## hard-codes `dtype="<f8"` at io/binary.py:67 — so this assertion is the only place
  ## the declaration is tied to reality.
  m <- batchdx_distinct_matrix(3L, 7L)
  path <- file.path(batchdx_workdir("sidecar"), "m.f64")
  write_binmat(path, m, row_ids = c("a", "b", "c"), col_names = paste0("c", 1:7))
  meta <- jsonlite::fromJSON(paste0(path, ".json"), simplifyVector = TRUE)

  expect_identical(meta$magic, BATCHDX_BINMAT_MAGIC)
  expect_identical(meta$dtype, "float64")
  expect_identical(meta$byte_order, "little")
  expect_identical(meta$layout, "row_major")
  expect_identical(as.integer(meta$n_rows), 3L)
  expect_identical(as.integer(meta$n_cols), 7L)
  expect_equal(as.integer(meta$n_rows) * as.integer(meta$n_cols) * 8, file.info(path)$size)
})

test_that("read_binmat decodes a row-major stream written by hand", {
  ## The other direction, and the one that actually runs in production: Python writes,
  ## R reads. The payload and sidecar are produced here without `write_binmat`, so a
  ## divergence in the reader cannot hide behind a matching divergence in the writer.
  m <- batchdx_distinct_matrix(3L, 7L)
  row_ids <- c("a", "b", "c")
  col_ids <- paste0("c", 1:7)
  dir <- batchdx_workdir("read")
  path <- file.path(dir, "m.f64")

  batchdx_write_raw(path, batchdx_expected_row_major_bytes(m))
  batchdx_write_sidecar(path, nrow(m), ncol(m), row_ids, col_ids)

  bm <- read_binmat(path)
  expect_identical(dim(bm$matrix), c(3L, 7L))
  expect_identical(as.vector(t(bm$matrix)), as.vector(t(m)))
  expect_identical(rownames(bm$matrix), row_ids)
  expect_identical(colnames(bm$matrix), col_ids)
})

test_that("a row-major/column-major confusion would be detected, not absorbed", {
  ## A negative control for the tests above. An all-ones or symmetric matrix survives a
  ## transpose unchanged, which is how a layout bug survives a naive round-trip test;
  ## so this asserts that the fixture actually discriminates, by checking that the
  ## transposed interpretation is a *different* matrix.
  m <- batchdx_distinct_matrix(3L, 7L)
  expect_false(identical(as.vector(m), as.vector(t(m))))
  dir <- batchdx_workdir("absorb")
  path <- file.path(dir, "m.f64")
  batchdx_write_raw(path, batchdx_expected_row_major_bytes(m))
  ## Declare the shape transposed: a reader that reshaped from the declared shape alone,
  ## rather than honouring the stream order, now returns a 7x3 matrix. The values in
  ## stream order must be untouched — that is what proves the reader read the bytes in
  ## the order they were written and then reshaped, not that it reinterpreted them.
  batchdx_write_sidecar(path, 7L, 3L, paste0("r", 1:7), paste0("c", 1:3))
  bm <- read_binmat(path)
  expect_identical(dim(bm$matrix), c(7L, 3L))
  expect_identical(as.vector(t(bm$matrix)), as.vector(t(m)))
  expect_identical(rownames(bm$matrix), paste0("r", 1:7))
})

test_that("row_ids and col_names are JSON arrays, never unboxed scalars", {
  ## A real defect, found by writing these tests.
  ##
  ## `write_binmat` serialises its sidecar with `jsonlite::toJSON(..., auto_unbox = TRUE)`.
  ## `auto_unbox` collapses any length-1 vector to a bare scalar, so a one-row matrix
  ## produces `"row_ids": "c1"` — a JSON *string* — where Python's `io/binary.py:48`
  ## does `list(meta["row_ids"])`, which iterates the string's characters and yields
  ## `['c', '1']`. Two bogus cell ids instead of one real one, silently.
  ##
  ## Severity, stated precisely rather than inflated: the only current consumer of an
  ## R-written matrix is `adapters/cycombine.py:193`, which compares `row_ids` against
  ## the dataset's cell ids at :194-197 and raises `BackendError` on a mismatch. So
  ## today this is a loud failure at the adapter boundary rather than a silently
  ## corrupted correction. It is still a defect: `read_matrix` itself does not validate,
  ## a one-marker panel unboxes `col_names` the same way, and any second consumer of an
  ## R-written matrix inherits it unguarded.
  ##
  ## The fix is one line in `binmat.R` — drop `auto_unbox = TRUE`, or wrap the two id
  ## vectors in `I()`. This file asserts the contract so the fix is verifiable.
  dir <- batchdx_workdir("unbox")
  for (case in list(
    list(name = "one row",    m = batchdx_distinct_matrix(1L, 3L), n = 1L, k = 3L),
    list(name = "one column", m = batchdx_distinct_matrix(3L, 1L), n = 3L, k = 1L),
    list(name = "one by one", m = batchdx_distinct_matrix(1L, 1L), n = 1L, k = 1L)
  )) {
    path <- file.path(dir, paste0(case$n, "x", case$k, ".f64"))
    row_ids <- paste0("cell", seq_len(case$n))
    col_names <- paste0("marker", seq_len(case$k))
    write_binmat(path, case$m, row_ids = row_ids, col_names = col_names)

    ## Read without simplifyVector, because that is what exposes the JSON type. With
    ## simplifyVector = TRUE a one-element array and a scalar are indistinguishable,
    ## which is precisely how this reached a commit.
    raw_meta <- jsonlite::fromJSON(paste0(path, ".json"), simplifyVector = FALSE)
    expect_true(is.list(raw_meta$row_ids),
                info = sprintf(paste0("%s: row_ids serialised as %s, not a JSON array; ",
                                      "Python's list(meta['row_ids']) would iterate its ",
                                      "characters"),
                               case$name, class(raw_meta$row_ids)[[1]]))
    expect_true(is.list(raw_meta$col_names),
                info = sprintf("%s: col_names serialised as %s, not a JSON array",
                               case$name, class(raw_meta$col_names)[[1]]))
    expect_identical(unlist(raw_meta$row_ids), row_ids)
    expect_identical(unlist(raw_meta$col_names), col_names)

    ## And the round trip must still hold through the public reader.
    bm <- read_binmat(path)
    expect_identical(rownames(bm$matrix), row_ids)
    expect_identical(colnames(bm$matrix), col_names)
  }
})

test_that("multi-row and multi-column sides stay arrays", {
  ## The control for the test above: the same assertion on shapes that do not trigger
  ## `auto_unbox`, so a blanket "always a list" change cannot pass by accident.
  dir <- batchdx_workdir("unbox-control")
  path <- file.path(dir, "m.f64")
  write_binmat(path, batchdx_distinct_matrix(4L, 6L),
               row_ids = paste0("cell", 1:4), col_names = paste0("marker", 1:6))
  raw_meta <- jsonlite::fromJSON(paste0(path, ".json"), simplifyVector = FALSE)
  expect_true(is.list(raw_meta$row_ids))
  expect_true(is.list(raw_meta$col_names))
  expect_length(raw_meta$row_ids, 4L)
  expect_length(raw_meta$col_names, 6L)
})

test_that("a round trip preserves every value exactly, across shapes", {
  ## Exact equality throughout: both sides are IEEE-754 float64 written as raw bytes, so
  ## any difference is a defect rather than floating-point noise. Non-square shapes are
  ## included because a square one hides a transpose.
  for (shape in list(c(2L, 3L), c(3L, 2L), c(5L, 1L), c(1L, 5L), c(6L, 4L), c(1L, 1L))) {
    m <- batchdx_distinct_matrix(shape[[1]], shape[[2]])
    path <- file.path(batchdx_workdir("roundtrip"), "m.f64")
    write_binmat(path, m, row_ids = paste0("r", seq_len(shape[[1]])),
                 col_names = paste0("c", seq_len(shape[[2]])))
    bm <- read_binmat(path)
    expect_identical(as.vector(bm$matrix), as.vector(m),
                     info = sprintf("round trip lost values at shape %dx%d",
                                    shape[[1]], shape[[2]]))
  }
})

test_that("special float64 values survive the exchange", {
  ## A cytometry export can legitimately contain NaN (a compensated channel with no
  ## events) and +/-Inf. These are the values most likely to be mangled by a dtype or
  ## endianness mistake, and least likely to be noticed, because NaN still prints as
  ## NaN.
  m <- matrix(c(NaN, Inf, -Inf, 1e308, -1e308, 5e-324, 0, -0, 1), nrow = 3, byrow = TRUE)
  path <- file.path(batchdx_workdir("special"), "m.f64")
  write_binmat(path, m, row_ids = c("a", "b", "c"), col_names = c("x", "y", "z"))
  bm <- read_binmat(path)
  expect_identical(is.nan(bm$matrix[1, 1]), TRUE)
  expect_identical(bm$matrix[1, 2], Inf)
  expect_identical(bm$matrix[1, 3], -Inf)
  expect_identical(bm$matrix[2, 1], 1e308)
  expect_identical(bm$matrix[2, 3], 5e-324)
  expect_identical(bm$matrix[3, 1], 0)
  expect_identical(bm$matrix[3, 2], 0)
  ## Byte identity is the real assertion: -0 and 0 differ in the stream and must stay
  ## distinguishable.
  expect_identical(readBin(path, what = "raw", n = file.info(path)$size),
                   batchdx_expected_row_major_bytes(m))
})

test_that("read_binmat refuses a missing sidecar, a wrong magic and a short payload", {
  dir <- batchdx_workdir("refuse")
  m <- batchdx_distinct_matrix(3L, 4L)

  ## No sidecar.
  orphan <- file.path(dir, "orphan.f64")
  batchdx_write_raw(orphan, batchdx_expected_row_major_bytes(m))
  expect_error(read_binmat(orphan), "sidecar not found")

  ## Wrong magic: a file from some other producer must not be silently reinterpreted.
  wrong <- file.path(dir, "wrong.f64")
  batchdx_write_raw(wrong, batchdx_expected_row_major_bytes(m))
  batchdx_write_sidecar(wrong, 3L, 4L, c("a", "b", "c"), c("w", "x", "y", "z"),
                        magic = "some-other-format-0")
  expect_error(read_binmat(wrong), "unexpected magic")

  ## Truncated payload: the shape in the sidecar promises more values than the file
  ## holds. Without the length check at binmat.R:22-24 a short read would be padded or
  ## silently reshaped.
  short <- file.path(dir, "short.f64")
  batchdx_write_raw(short, batchdx_expected_row_major_bytes(m)[1:20])
  batchdx_write_sidecar(short, 3L, 4L, c("a", "b", "c"), c("w", "x", "y", "z"))
  expect_error(read_binmat(short), "expected 12 values")
})

test_that("write_binmat refuses a shape that disagrees with its ids", {
  path <- file.path(batchdx_workdir("shape"), "m.f64")
  expect_error(write_binmat(path, matrix(0, nrow = 2, ncol = 2), c("only-one"), c("x", "y")))
  expect_error(write_binmat(path, matrix(0, nrow = 2, ncol = 2), c("a", "b"), "only-one"))
})

test_that("as_opt normalises an absent JSON value to NULL", {
  ## `batchdx_write_result` serialises with `null = "null"`, and the adapter scripts read
  ## optional upstream parameters through `as_opt`. A literal JSON `null` must arrive as
  ## R's NULL rather than as a length-0 vector or the string "null", or an optional
  ## parameter is passed to the backend as a value the backend never expected.
  expect_null(as_opt(NULL))
  expect_null(as_opt(NA))
  expect_null(as_opt(NA_character_))
  expect_identical(as_opt(5), 5)
  expect_identical(as_opt("x"), "x")
  expect_identical(as_opt(c(1, 2)), c(1, 2))
  ## Recorded as the implementation actually behaves, not as one might wish it did: a
  ## length-0 vector passes through rather than becoming NULL. It is unreachable from
  ## the adapter scripts, because `batchdx_write_result` serialises an absent value as
  ## JSON `null` and jsonlite reads that back as a real NULL, which is the case asserted
  ## above. Stated here so the gap is on the record if that ever changes.
  expect_identical(as_opt(character(0)), character(0))
})


# --------------------------------------------------------------------------------------
# Cross-language agreement
# --------------------------------------------------------------------------------------
#
# Everything above holds R to a contract stated in R. These hold it to the contract the
# *Python* side actually implements, which is the only thing that makes the mirror a
# mirror. They need a Python that can import batchdx.cyto; the skip below names the
# condition rather than hiding it, and `tests/run_r_tests.R` reports the count so a
# skipped cross-language test is visible in the output.

batchdx_python_or_skip <- function() {
  py <- batchdx_python()
  if (is.null(py)) {
    testthat::skip("no Python that can import batchdx.cyto was found; set BATCHDX_PYTHON")
  }
  py
}

## Write a matrix with the Python implementation and return the path, so the R side
## reads a file produced entirely by the other implementation.
##
## `paste(..., sep = "\n")` is load-bearing: the default separator is a space, which
## collapses the script into one line and produces a Python syntax error rather than a
## diagnosis of anything.
batchdx_write_with_python <- function(py, dir, m, row_ids, col_names) {
  path <- file.path(dir, "from-python.f64")
  script <- paste(
    "import sys",
    sprintf("sys.path.insert(0, %s)", shQuote(file.path(BATCHDX_ROOT, "src"))),
    "from batchdx.cyto.io.binary import write_matrix",
    sprintf("m = __import__('numpy').array(%s, dtype='float64').reshape(%d, %d)",
            batchdx_python_list(as.vector(t(as.matrix(m)))), nrow(m), ncol(m)),
    sprintf("write_matrix(%s, m, %s, %s)", shQuote(path),
            batchdx_python_list(row_ids), batchdx_python_list(col_names)),
    sep = "\n"
  )
  out <- system2(py, c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("the Python writer failed (exit ", status, "):\n", paste(out, collapse = "\n"))
  }
  path
}

test_that("a matrix written by Python's write_matrix reads back in R unchanged", {
  py <- batchdx_python_or_skip()
  for (shape in list(c(3L, 7L), c(7L, 3L), c(1L, 4L), c(4L, 1L), c(1L, 1L))) {
    m <- batchdx_distinct_matrix(shape[[1]], shape[[2]])
    row_ids <- paste0("cell", seq_len(shape[[1]]))
    col_names <- paste0("marker", seq_len(shape[[2]]))
    written <- batchdx_write_with_python(py, batchdx_workdir("xlang"), m, row_ids, col_names)

    bm <- read_binmat(written)
    expect_identical(dim(bm$matrix), c(shape[[1]], shape[[2]]))
    expect_identical(as.vector(bm$matrix), as.vector(m),
                     info = sprintf("Python wrote %dx%d and R read something else",
                                    shape[[1]], shape[[2]]))
    expect_identical(rownames(bm$matrix), row_ids)
    expect_identical(colnames(bm$matrix), col_names)
  }
})

test_that("a matrix written by R's write_binmat reads back in Python unchanged", {
  ## The return trip, and the one that carries the production payload:
  ## `cycombine_correct.R:99` calls `write_binmat` and `adapters/cycombine.py:193` reads
  ## the result. This asserts the two agree on the bytes and on the identity of the rows
  ## and columns, for the degenerate shapes where they currently do not.
  py <- batchdx_python_or_skip()
  for (shape in list(c(3L, 7L), c(4L, 1L), c(1L, 1L))) {
    m <- batchdx_distinct_matrix(shape[[1]], shape[[2]])
    row_ids <- paste0("cell", seq_len(shape[[1]]))
    col_names <- paste0("marker", seq_len(shape[[2]]))
    path <- file.path(batchdx_workdir("xlang-back"), "from-r.f64")
    write_binmat(path, m, row_ids = row_ids, col_names = col_names)

    script <- paste(
      "import sys, json",
      sprintf("sys.path.insert(0, %s)", shQuote(file.path(BATCHDX_ROOT, "src"))),
      "from batchdx.cyto.io.binary import read_matrix",
      sprintf("m, row_ids, col_names, meta = read_matrix(%s)", shQuote(path)),
      "print(json.dumps({'shape': list(m.shape), 'rows': row_ids,",
      "                  'cols': col_names, 'values': m.ravel().tolist()}))",
      sep = "\n"
    )
    out <- system2(py, c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    if (!is.null(status) && status != 0L) {
      stop("the Python reader failed (exit ", status, "):\n", paste(out, collapse = "\n"))
    }
    got <- jsonlite::fromJSON(out[[length(out)]], simplifyVector = TRUE)

    expect_identical(as.integer(got$shape), c(shape[[1]], shape[[2]]))
    expect_identical(got$rows, row_ids,
                     info = sprintf("Python read row ids %s, R wrote %s",
                                    paste(got$rows, collapse = ","),
                                    paste(row_ids, collapse = ",")))
    expect_identical(got$cols, col_names)
    ## `as.vector(t(m))`, not `as.vector(m)`: R flattens column-major by default, and
    ## using the unqualified form here is the very confusion these tests exist to catch —
    ## in the test rather than in the code, which is the appropriate place to catch it.
    expect_identical(unlist(got$values), as.vector(t(m)),
                     info = sprintf("the two implementations disagree on the bytes at %dx%d",
                                    shape[[1]], shape[[2]]))
  }
})
