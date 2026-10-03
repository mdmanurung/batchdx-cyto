#!/usr/bin/env Rscript
## Run the R-side test suite.
##
##     Rscript --vanilla tests/run_r_tests.R
##
## There are five tracked R files that `batchdx.cyto.pipeline` drives and that only one
## of — `scripts/audit_upstream.R` — was executed by CI before this. This runner is the
## single command that makes the other four testable, so wiring it into CI is a
## one-line change rather than a judgement call made inside a YAML file.
##
## Run it with `--vanilla`. A `~/.Renviron` on the build host can redirect `R_LIBS_USER`
## into a personal library, which is how an R environment becomes unattributable; the
## same hazard is recorded against `scripts/setup_r_env.sh` and `run_r_script`.
##
## Exit status: 0 only if every test passed. A skipped test is reported on stdout and
## counted in the summary, so a suite that skips everything cannot read as a pass.
##
## Environment:
##   BATCHDX_PYTHON   interpreter used by the cross-language tests. Defaults to the
##                    repository's `.venv/bin/python`, then `python3`. The cross-language
##                    tests skip, with that reason printed, if none of them can import
##                    `batchdx.cyto`.

suppressPackageStartupMessages(library(testthat))

args <- commandArgs(trailingOnly = TRUE)
filter <- if (length(args) >= 1L) args[[1]] else NULL

## Locate the repository root from this script's own path, not from the working
## directory, so the runner works from anywhere.
script_path <- NULL
for (candidate in commandArgs(trailingOnly = FALSE)) {
  if (startsWith(candidate, "--file=")) {
    script_path <- substring(candidate, nchar("--file=") + 1L)
  }
}
root <- if (is.null(script_path)) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

if (!dir.exists(file.path(root, "tests", "testthat"))) {
  stop("no tests/testthat directory under ", root)
}
setwd(root)

message("batchdx-cyto R test suite")
message("  root:     ", root)
message("  testthat: ", as.character(utils::packageVersion("testthat")))
message("  R:        ", R.version.string)
message("  python:   ", {
  py <- Sys.getenv("BATCHDX_PYTHON")
  if (nzchar(py)) py else file.path(root, ".venv", "bin", "python")
})
message("")

## Report every failure, not the first ten. A truncated report invites the reader to
## assume the tenth failure is the tenth *kind* of failure.
options(testthat.progress.max_fails = 1000L)

results <- test_dir(
  file.path(root, "tests", "testthat"),
  filter = filter,
  reporter = "summary",
  stop_on_failure = FALSE,
  stop_on_warning = FALSE
)

df <- as.data.frame(results)
skipped <- sum(df$skipped)
failed <- sum(df$failed)
errored <- sum(df$error)
passed <- sum(df$passed)

cat("\n")
cat(sprintf("passed %d, failed %d, errored %d, skipped %d, warnings %d\n",
            passed, failed, errored, skipped, sum(df$warning)))

if (skipped > 0L) {
  cat("\nSkipped tests (each names its reason in the output above):\n")
  for (i in seq_len(nrow(df))) {
    if (df$skipped[i] > 0L) {
      cat(sprintf("  %s  (%s)\n", df$file[i], df$test[i]))
    }
  }
  cat(
    "\nA skip is a missing environment, not a pass. In particular, the cross-language\n",
    "tests in test-binmat.R skip unless a Python that can import batchdx.cyto is found;\n",
    "set BATCHDX_PYTHON to run them.\n"
  )
}

if (failed > 0L || errored > 0L) quit(status = 1L)
quit(status = 0L)
