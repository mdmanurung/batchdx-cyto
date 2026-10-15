## CytoNorm correction, non-interactive.
##
## Route: CytoNorm.train on the training split's measured reference aliquots, then
## CytoNorm.normalize over the rest.  Fit and apply are genuinely separable in the
## pinned API (config/upstream_lock.yaml records the commit), which is the whole
## reason this backend exists.
##
## CytoNorm is GPL (>= 2) and is not redistributed with this package.  It is called
## from a separate R process at arm's length, exactly as
## docs/decisions/0003-cytobatchflagr-integration.md requires; nothing from its
## source is copied here.
##
## FILES, NOT DATA FRAMES.  Both entry points take FCS file paths.  The Python
## adapter therefore writes the baseline out as per-sample FCS first
## (batchdx.cyto.io.fcs), hands this script the resulting paths, and reads the
## normalised files back afterwards.  The intermediate matrix exchange format in
## lib/binmat.R is deliberately not used here: it would fit CytoNorm to the wrong
## interface.

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

## Warnings are part of the result, not console noise.  CytoNorm falls back to the
## identity function for any cluster whose per-batch quantiles are degenerate
## ("Not enough cells for <batch>"), which means those cells come back unchanged.
## Reporting that as a successful correction without saying it would be exactly the
## failure this package is built to avoid, so every warning is captured and handed
## back to the adapter.
collect_warnings <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  ## The declared baseline is ALREADY transformed (asinh or logicle, named in the
  ## processing manifest).  CytoNorm's own vignettes pass a flowCore transformList
  ## because they start from raw fcs files.  Applying one here would transform the
  ## values a second time, so the transform list is the identity and this script
  ## refuses any other declaration rather than quietly ignoring it.
  transform <- as.character(p$transform_list)
  if (!identical(transform, "identity")) {
    stop(sprintf(
      paste0("refusing transform_list '%s': this adapter writes the already-transformed ",
             "preprocessed baseline as FCS, so the transform list handed to CytoNorm is ",
             "the identity. Passing a flowCore transform would transform the values twice."),
      transform))
  }
  if (!requireNamespace("CytoNorm", quietly = TRUE)) {
    stop("the CytoNorm R package is not installed in this environment")
  }
  ## Namespaced, not attached: attaching would mask any other package that exports
  ## a function of the same name in this session.
  stopifnot(is.function(CytoNorm::CytoNorm.train), is.function(CytoNorm::CytoNorm.normalize))

  channels <- as.character(p$channels)
  if (!length(channels)) stop("no channel names were supplied")

  train_files <- as.character(p$train_files)
  train_labels <- as.character(p$train_labels)
  apply_files <- as.character(p$apply_files)
  apply_labels <- as.character(p$apply_labels)
  if (length(train_files) != length(train_labels)) {
    stop("train_files and train_labels have different lengths")
  }
  if (length(apply_files) != length(apply_labels)) {
    stop("apply_files and apply_labels have different lengths")
  }
  if (length(unique(train_labels)) < 2) {
    stop(sprintf(
      paste0("the training split holds a measured reference in %d batch(es); CytoNorm fits ",
             "a per-batch quantile goal and needs at least two. This is a coverage failure ",
             "of the split, not a defect in CytoNorm, and it is not papered over by fitting ",
             "on a single batch."),
      length(unique(train_labels))))
  }
  for (f in c(train_files, apply_files)) {
    if (!file.exists(f)) stop(sprintf("input FCS file not found: %s", f))
  }

  train_dir <- as.character(p$train_dir)
  output_dir <- as.character(p$output_dir)
  dir.create(train_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  flowsom_params <- list(
    nCells = as.integer(p$nCells),
    xdim = as.integer(p$xdim),
    ydim = as.integer(p$ydim),
    nClus = as.integer(p$nClus),
    scale = isTRUE(p$scale)
  )
  norm_params <- list(nQ = as.integer(p$nQ))
  prefix <- if (is.null(p$prefix)) "Norm_" else as.character(p$prefix)

  set.seed(as.integer(p$seed))

  trained <- collect_warnings(
    CytoNorm::CytoNorm.train(
      files = train_files,
      labels = train_labels,
      channels = channels,
      transformList = NULL,
      outputDir = train_dir,
      FlowSOM.params = flowsom_params,
      normParams = norm_params,
      seed = as.integer(p$seed),
      clean = TRUE,
      plot = FALSE,
      verbose = FALSE,
      recompute = FALSE
    )
  )
  model <- trained$value
  metaclusters <- sort(unique(model$fsom$metaclustering))
  if (is.null(model$fsom) || is.null(model$clusterRes)) {
    stop("CytoNorm.train returned a model with no FlowSOM fit or no cluster splines")
  }

  ## CytoNorm.normalize returns a flowSet when given a flowSet and NULL when given
  ## file paths; on the file route the results exist only on disk.  That is why the
  ## adapter reads the written files back rather than trusting a return value.
  ## On the file route CytoNorm.normalize returns NULL and the results exist only on
  ## disk, so the warnings are what this script has to hand back here.
  normalised <- collect_warnings(
    CytoNorm::CytoNorm.normalize(
      model = model,
      files = apply_files,
      labels = apply_labels,
      transformList = NULL,
      transformList.reverse = NULL,
      outputDir = output_dir,
      prefix = prefix,
      clean = TRUE,
      verbose = FALSE,
      write = TRUE
    )
  )

  expected <- file.path(output_dir, paste0(prefix, basename(apply_files)))
  missing <- expected[!file.exists(expected)]
  if (length(missing)) {
    stop(sprintf(
      "CytoNorm.normalize reported success but wrote no output for %d of %d file(s); first: %s",
      length(missing), length(apply_files), missing[[1]]))
  }

  ## WHICH BATCHES WERE ACTUALLY FITTED.  CytoNorm fits one quantile spline per
  ## cluster per training label.  When it is asked to apply a label it has no model
  ## for it skips that file, and the final combine step then writes it back with no
  ## cluster replaced -- i.e. the output is the input, byte for byte.  A caller that
  ## only checks "did it write a file" would report that as a successful correction.
  ## So the fitted labels are reported explicitly and the unfitted ones named, and
  ## the Python adapter treats an unchanged file as a partial correction.
  ## clusterRes[[cluster]] is a list of channels / splines / quantiles / ... and the
  ## per-label models live in its $splines element, keyed by label.  Reading
  ## names(clusterRes[[cluster]]) instead returns those element names, which are the
  ## same for every cluster and say nothing about which batches were fitted.
  fitted_labels <- sort(unique(unlist(
    lapply(model$clusterRes, function(cluster_model) {
      names(cluster_model$splines)
    })
  )))
  unfitted_labels <- setdiff(sort(unique(apply_labels)), fitted_labels)
  if (length(unfitted_labels)) {
    warning(sprintf(
      paste0("no quantile model was fitted for label(s) %s: the training split held no ",
             "measured reference there, so every file in those batches was written back ",
             "unchanged. This is a partial correction."),
      paste(unfitted_labels, collapse = ", ")))
  }

  ## Housekeeping noise from upstream's own cleanup: file.remove() is called
  ## unconditionally on per-cluster files that may never have been written, and R
  ## warns for each.  Counted, not repeated eighteen times in the record, and not
  ## mixed in with the scientific warnings above.
  housekeeping <- grepl("cannot remove file", normalised$warnings)
  scientific <- normalised$warnings[!housekeeping]
  identity_warnings <- grepl("identity function will be used", scientific)

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    n_train_files = length(train_files),
    n_train_batches = length(unique(train_labels)),
    n_apply_files = length(apply_files),
    n_metaclusters = length(metaclusters),
    metaclusters = as.character(metaclusters),
    outputs = as.character(expected),
    ## as.list() so a single label still serialises as a JSON array.  With
    ## auto_unbox = TRUE a length-one character vector collapses to a bare string, and
    ## a reader that iterates it would then iterate its *characters* -- turning the
    ## batch name "B3" into ['3', 'B'].
    fitted_labels = as.list(as.character(fitted_labels)),
    unfitted_labels = as.list(as.character(unfitted_labels)),
    n_fitted_batches = length(fitted_labels),
    ## Clusters left on the identity function are clusters whose values came back
    ## unchanged.  The count is reported so a reader can tell a partial correction
    ## from a complete one instead of inferring it from a changed-cell count.
    n_identity_function_clusters = sum(identity_warnings),
    warnings = c(trained$warnings, scientific),
    n_housekeeping_warnings = sum(housekeeping),
    cytonorm_version = as.character(utils::packageVersion("CytoNorm")),
    flowsom_version = as.character(utils::packageVersion("FlowSOM")),
    flowcore_version = as.character(utils::packageVersion("flowCore")),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
