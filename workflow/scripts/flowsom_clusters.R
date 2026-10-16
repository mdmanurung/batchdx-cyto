## FlowSOM mapping-similarity clustering, non-interactive.
##
## Route: FlowSOM::FlowSOM on the REFERENCE batch's cells, then FlowSOM::NewData to map
## every other batch's cells through that tree.  This is the clustering step of the
## pinned OpenProblems component src/metrics/flowsom_mapping_similarity (script.R at
## commit 37151479): the tree is built once and the second side is *mapped* into it, so
## both sides are labelled by the same clusters and their compositions are comparable.
## The pinned script reads its three parameters (num_clusters, som_xdim, som_ydim) from
## the benchmark's unintegrated file; this package has no such file, so they are declared
## adaptation parameters passed in params.json and named on the metric spec.
##
## Only the cluster ASSIGNMENTS cross the boundary back: the similarity arithmetic
## (weighted absolute population-percentage differences) is computed on the Python side,
## where it is verified against the pinned helper.R by
## tests/scientific/test_stage7_openproblems.py.  The assignments are serialised as
## strings because the r_runner result validator requires arrays of strings; Python
## converts them back to integers.

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

  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    batchdx_write_result(parsed$result_path, list(
      status = "not_available",
      reason = paste("FlowSOM is not installed in this R environment;",
                     "the mapping-similarity row cannot be computed"),
      missing_packages = "FlowSOM"
    ))
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(FlowSOM))

  reference <- read_binmat(p$reference_matrix)
  mapped <- read_binmat(p$mapped_matrix)
  if (nrow(reference$matrix) != length(p$reference_labels)) {
    stop(sprintf("reference frame has %d rows but %d labels were supplied",
                 nrow(reference$matrix), length(p$reference_labels)))
  }
  if (nrow(mapped$matrix) != length(p$mapped_labels)) {
    stop(sprintf("mapped frame has %d rows but %d labels were supplied",
                 nrow(mapped$matrix), length(p$mapped_labels)))
  }
  if (!identical(colnames(reference$matrix), colnames(mapped$matrix))) {
    stop("the reference and mapped frames do not carry the same markers in the same order")
  }

  markers <- colnames(reference$matrix)
  num_clusters <- as.integer(p$num_clusters)
  seed <- as.integer(p$seed)

  ## The tree is built on the reference side only, with the pinned seed: the mapped side
  ## must be labelled by a tree it did not contribute to, or the comparison is circular.
  fsom <- FlowSOM::FlowSOM(
    flowCore::flowFrame(reference$matrix),
    colsToUse = markers,
    nClus = num_clusters,
    xdim = as.integer(p$som_xdim),
    ydim = as.integer(p$som_ydim),
    seed = seed
  )
  mapped_fsom <- FlowSOM::NewData(fsom, flowCore::flowFrame(mapped$matrix))

  ## as.character keeps the payload an array of strings for the r_runner validator;
  ## the cluster ids themselves are the integers GetClusters returns.
  clusters_reference <- as.character(FlowSOM::GetClusters(fsom))
  clusters_mapped <- as.character(FlowSOM::GetClusters(mapped_fsom))
  if (length(clusters_reference) != nrow(reference$matrix) ||
      length(clusters_mapped) != nrow(mapped$matrix)) {
    stop("FlowSOM returned an assignment vector that does not match its input frame")
  }

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    clusters_reference = unname(as.list(clusters_reference)),
    clusters_mapped = unname(as.list(clusters_mapped)),
    warnings = character(0)
  ))
  invisible(NULL)
}

main()
