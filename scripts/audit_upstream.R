#!/usr/bin/env Rscript
## Stage 0 audit, R side.
##
## Usage: Rscript scripts/audit_upstream.R <out.json> [upstream_dir] [lock.yaml]
##
## For every R resource the audit reports the resolved version and commit, enumerates
## the exported functions, and prints the formals of every function an adapter calls.
##
## Two modes, and the mode is recorded in the output:
##   "installed"  the package is on the library path; formals come from the loaded
##                function objects.
##   "source"     the package is not installed; the pinned checkout is PARSED (never
##                evaluated) and the formals come from the assignment expressions in
##                the source at that commit.
## Neither mode infers a signature: both read one.  A function an adapter calls that
## appears in neither is reported as missing, and the script exits non-zero.

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else "docs/audit/audit_r.json"
upstream_dir <- if (length(args) >= 2) args[[2]] else ".upstream"
lock_path <- if (length(args) >= 3) args[[3]] else "config/upstream_lock.yaml"

## Minimal reader for the few fields we need, so the audit does not depend on an R
## YAML package that may not be installed.
read_adapter_calls <- function(path) {
  if (!file.exists(path)) return(list())
  lines <- readLines(path, warn = FALSE)
  out <- list()
  current <- NULL
  in_calls <- FALSE
  for (line in lines) {
    if (grepl("^  [a-z_]+:\\s*$", line)) {
      current <- trimws(sub(":", "", line))
      in_calls <- FALSE
    } else if (grepl("^    adapter_calls:\\s*$", line)) {
      in_calls <- TRUE
      out[[current]] <- character(0)
    } else if (in_calls && grepl("^      - ", line)) {
      out[[current]] <- c(out[[current]], trimws(sub("^      - ", "", line)))
    } else if (grepl("^    [a-z_]+:", line)) {
      in_calls <- FALSE
    }
  }
  out
}

git_commit <- function(dir) {
  if (!dir.exists(file.path(dir, ".git"))) return(NA_character_)
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", dir, "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_)
  if (length(out) == 0) NA_character_ else out[[1]]
}

## Parse (do not evaluate) every R file and return function name -> formals text.
parse_function_signatures <- function(r_dir) {
  files <- list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  out <- list()
  for (f in files) {
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) {
      if (!is.call(e)) next
      if (!(as.character(e[[1]])[1] %in% c("<-", "=", "assign"))) next
      lhs <- e[[2]]
      rhs <- e[[3]]
      if (!is.name(lhs) || !is.call(rhs)) next
      if (as.character(rhs[[1]])[1] != "function") next
      name <- as.character(lhs)
      formals_txt <- paste(deparse(rhs[[2]]), collapse = " ")
      formals_txt <- sub("^as\\.pairlist\\(alist\\(", "", formals_txt)
      formals_txt <- sub("^pairlist\\(", "", formals_txt)
      formals_txt <- sub("\\)\\)$", "", formals_txt)
      formals_txt <- gsub("\\s+", " ", formals_txt)
      formals_txt <- gsub(", \\.\\.\\. = $", ", ...", formals_txt)
      out[[name]] <- list(
        signature = sprintf("%s(%s)", name, formals_txt),
        arguments = names(rhs[[2]]),
        file = basename(f)
      )
    }
  }
  out
}

read_namespace_exports <- function(dir) {
  ns <- file.path(dir, "NAMESPACE")
  if (!file.exists(ns)) return(character(0))
  lines <- readLines(ns, warn = FALSE)
  exports <- regmatches(lines, regexpr("(?<=^export\\()[^)]+", lines, perl = TRUE))
  gsub('"', "", exports)
}

read_description_field <- function(dir, field) {
  d <- file.path(dir, "DESCRIPTION")
  if (!file.exists(d)) return(NA_character_)
  value <- tryCatch(unname(read.dcf(d, fields = field)[1, 1]), error = function(e) NA_character_)
  value
}

audit_resource <- function(key, dir_name, package_name, adapter_calls) {
  dir <- file.path(upstream_dir, dir_name)
  entry <- list(
    resource = key,
    checkout = dir,
    checkout_present = dir.exists(dir),
    commit = git_commit(dir),
    adapter_calls = adapter_calls
  )

  installed <- !is.null(package_name) &&
    requireNamespace(package_name, quietly = TRUE)

  if (installed) {
    entry$mode <- "installed"
    entry$version <- as.character(utils::packageVersion(package_name))
    entry$exports <- sort(getNamespaceExports(package_name))
    signatures <- list()
    for (fn in adapter_calls) {
      if (!exists(fn, envir = asNamespace(package_name), inherits = FALSE)) next
      obj <- get(fn, envir = asNamespace(package_name))
      if (!is.function(obj)) next
      signatures[[fn]] <- list(
        signature = paste(trimws(deparse(args(obj))[1:(length(deparse(args(obj))) - 1)]),
                          collapse = " "),
        arguments = names(formals(obj)),
        file = NA_character_
      )
    }
    entry$signatures <- signatures
  } else if (entry$checkout_present) {
    entry$mode <- "source"
    entry$version <- read_description_field(dir, "Version")
    entry$exports <- sort(read_namespace_exports(dir))
    parsed <- parse_function_signatures(file.path(dir, "R"))
    entry$signatures <- parsed[intersect(names(parsed), adapter_calls)]
    entry$all_top_level_functions <- sort(names(parsed))
  } else {
    entry$mode <- "unavailable"
    entry$version <- NA_character_
    entry$exports <- character(0)
    entry$signatures <- list()
  }

  entry$missing <- setdiff(adapter_calls, names(entry$signatures))
  entry
}

calls <- read_adapter_calls(lock_path)

resources <- list(
  list(key = "cycombine", dir = "cyCombine", package = "cyCombine"),
  list(key = "cytonorm", dir = "CytoNorm", package = "CytoNorm"),
  list(key = "cytobatchflagr", dir = "CytoBatchFlagR", package = NULL)
)

results <- list()
for (r in resources) {
  results[[r$key]] <- audit_resource(
    r$key, r$dir, r$package,
    if (is.null(calls[[r$key]])) character(0) else calls[[r$key]]
  )
}

payload <- list(
  generated_by = "scripts/audit_upstream.R",
  r_version = R.version.string,
  jsonlite_version = as.character(utils::packageVersion("jsonlite")),
  upstream_dir = upstream_dir,
  lock = lock_path,
  resources = results
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null"), out_path)

cat("== batchdx-cyto upstream audit (R) ==\n")
cat(R.version.string, "\n\n")
missing_any <- FALSE
for (key in names(results)) {
  e <- results[[key]]
  cat(sprintf("-- %s [%s] version=%s commit=%s\n", key, e$mode,
              ifelse(is.na(e$version), "unknown", e$version),
              ifelse(is.na(e$commit), "unknown", substr(e$commit, 1, 12))))
  if (length(e$exports)) {
    cat("   exports:", paste(e$exports, collapse = ", "), "\n")
  }
  for (fn in names(e$signatures)) {
    cat("   ", e$signatures[[fn]]$signature, "\n", sep = "")
  }
  if (length(e$missing)) {
    missing_any <- TRUE
    cat("   MISSING (referenced by an adapter, not found):",
        paste(e$missing, collapse = ", "), "\n")
  }
  cat("\n")
}
cat("audit written to", out_path, "\n")

if (missing_any) {
  cat("\nFAIL: at least one adapter-referenced function was not found in the audited sources.\n")
  quit(status = 1)
}
