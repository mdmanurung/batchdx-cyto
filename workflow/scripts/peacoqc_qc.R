## PeacoQC event-level quality control, non-interactive.
##
## Entry point: PeacoQC::PeacoQC, from the Bioconductor package installed in this
## environment (bioconductor-peacoqc 1.16.0, the version the pinned SPiSCy reference
## resolves in .upstream/spiscy/docker/rbase.yaml:20).
##
## WHY THIS SCRIPT EXISTS WHEN .upstream/spiscy/workflow/scripts/QC.R ALREADY CALLS IT
## ------------------------------------------------------------------------------------
## QC.R calls PeacoQC and then throws the removal index away: it renames the report file
## (QC.R:150-160) and the gated flowFrame is simply a smaller file. That is adequate for
## producing a matrix and indefensible for a package whose whole premise is that nothing
## is dropped without being accounted for. This script keeps the index and hands it back.
##
## THERE IS NO $BadCell
## --------------------
## PeacoQC 1.16.0's return value has no BadCell element. The per-cell decision lives in
## GoodCells, a logical vector aligned to the input event order, so the removed set is
## !GoodCells. The element name is CHECKED rather than assumed, and a result without it
## is a hard error: a version that dropped GoodCells would otherwise make this script
## report "nothing was removed", which is the one answer this package must never give
## when it does not know.
##
## INSUFFICIENT DATA IS A STATUS, NOT A ZERO
## ------------------------------------------
## When PeacoQC declines to run the isolation tree -- "not enough bins for a robust
## isolation tree analysis" -- it sets ITPercentage to NA and then removes on MAD alone.
## determine_good_cells = "all" therefore completes, removes events, and has performed
## no IT analysis at all. That is reported as it_available = FALSE with the reason, so a
## MAD-only pass is never read as the "all" that was requested.

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

## PeacoQC's warnings are part of the result. "Not enough bins for a robust isolation
## tree analysis" is the one that decides whether this run did what was asked, so it is
## captured rather than left in the log where nothing reads it.
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

IT_REFUSAL <- "not enough bins for a robust isolation tree analysis"

main <- function() {
  parsed <- batchdx_args()
  p <- parsed$params
  started <- Sys.time()

  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop("the PeacoQC R package is not installed in this environment")
  }
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("the flowCore R package is not installed in this environment")
  }
  ## Namespaced, never attached: attaching would mask a same-named export elsewhere.
  stopifnot(is.function(PeacoQC::PeacoQC))

  files <- as.character(p$files)
  if (!length(files)) stop("no input files were supplied")
  absent <- files[!file.exists(files)]
  if (length(absent)) {
    stop(sprintf("input FCS file not found: %s", absent[[1]]))
  }
  ## Results are keyed by SAMPLE ID, not by file name. The reference implementation
  ## works in file names throughout, but this script hands the index back to a Python
  ## adapter that knows samples by the manifest's sample_id; keying by basename would
  ## make the two disagree by an extension, silently, on the one field the removal log
  ## is built from. Two parallel character vectors rather than a list of objects,
  ## because jsonlite's simplifyVector collapses a list of records into a data frame.
  sample_ids <- as.character(p$sample_ids)
  if (length(sample_ids) != length(files)) {
    stop("sample_ids and files have different lengths")
  }
  if (anyDuplicated(sample_ids)) {
    stop(sprintf("duplicate sample_id(s): %s",
                 paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")))
  }
  names(files) <- sample_ids

  channels <- as.character(p$channels)
  if (!length(channels)) {
    stop("no channel names were supplied; PeacoQC evaluates quality per channel and the list is never inferred")
  }
  criteria <- as.character(p$criteria)
  if (!criteria %in% c("all", "IT", "MAD")) {
    stop(sprintf("determine_good_cells='%s' is not one of all/IT/MAD", criteria))
  }
  if (criteria %in% c("all", "MAD") && is.null(p$mad)) {
    stop("determine_good_cells requires a non-null MAD value")
  }
  if (criteria %in% c("all", "IT") && is.null(p$it_limit)) {
    stop("determine_good_cells requires a non-null IT_limit value")
  }

  output_dir <- as.character(p$output_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ## The batch that matters. Only arguments the adapter actually varies are assembled
  ## here; every other formal keeps its upstream default, and those defaults are
  ## recorded in src/batchdx/cyto/preprocess/peacoqc.py::VERIFIED_FORMALS so a change
  ## upstream shows up in a diff of this repository rather than in a changed result.
  call <- list(
    channels = channels,
    determine_good_cells = criteria,
    ## save_fcs and plot differ from the upstream defaults (TRUE and 20) and from what
    ## the reference call site passes. Stated here rather than inherited: this adapter
    ## wants the index, not a second copy of the events.
    save_fcs = FALSE,
    plot = p$plot,
    report = isTRUE(p$report),
    output_directory = output_dir,
    name_directory = "PeacoQC_results",
    min_cells = as.integer(p$min_cells),
    max_bins = as.integer(p$max_bins),
    step = as.integer(p$step),
    consecutive_bins = as.integer(p$consecutive_bins),
    remove_zeros = isTRUE(p$remove_zeros),
    force_IT = as.integer(p$force_IT),
    peak_removal = as.numeric(p$peak_removal),
    min_nr_bins_peakdetection = as.integer(p$min_nr_bins_peakdetection)
  )
  if (!is.null(p$mad)) call$MAD <- as.numeric(p$mad)
  if (!is.null(p$it_limit)) call$IT_limit <- as.numeric(p$it_limit)
  ## time_channel_parameter defaults to "Time". Left unset when the caller declares no
  ## time channel, because naming a channel the file does not have would make PeacoQC
  ## fail on a parameter this adapter never needed.
  if (!is.null(as_opt(p$time_channel_parameter))) {
    call$time_channel_parameter <- as.character(p$time_channel_parameter)
  }

  results <- list()
  warnings_all <- character(0)
  ## PeacoQC writes its report into file.path(output_directory, name_directory) -- ONE
  ## directory for the whole call, appended to per sample, not one per sample. Verified
  ## from the installed package's own body rather than assumed, because this script
  ## reports the artifact it produced and a path it invented would be a path that does
  ## not exist.
  storing_directory <- file.path(output_dir, "PeacoQC_results")
  for (sample_id in names(files)) {
    f <- files[[sample_id]]
    label <- sample_id
    ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
    ## PeacoQC is chatty; everything it says is captured rather than streamed, so a
    ## diagnostic that decided the outcome cannot be lost in the noise.
    run <- collect_warnings(do.call(PeacoQC::PeacoQC, c(list(ff = ff), call)))
    res <- run$value
    warnings_all <- c(warnings_all, paste0(label, ": ", run$warnings))

    result_keys <- names(res)
    if (!"GoodCells" %in% result_keys) {
      stop(sprintf(
        paste0("PeacoQC returned no GoodCells element for %s; it returned [%s]. The ",
               "removal index is !GoodCells, so nothing here can say which events this ",
               "sample lost."),
        label, paste(result_keys, collapse = ", ")))
    }
    good <- as.logical(res$GoodCells)
    if (length(good) != nrow(ff)) {
      stop(sprintf("PeacoQC returned %d per-cell decisions for %d events in %s",
                   length(good), nrow(ff), label))
    }

    it_fraction <- as_opt(res$ITPercentage)
    mad_fraction <- as_opt(res$MADPercentage)
    it_available <- !is.null(it_fraction) && is.finite(as.numeric(it_fraction))
    it_reason <- if (it_available) NA_character_ else {
      hit <- grep(IT_REFUSAL, run$warnings, fixed = TRUE, value = TRUE)
      if (length(hit)) {
        paste(hit, collapse = "; ")
      } else {
        "the isolation tree reported no IT percentage; the warnings carried no reason"
      }
    }

    report_dir <- storing_directory
    written <- if (dir.exists(report_dir)) {
      list.files(report_dir, recursive = TRUE, full.names = TRUE)
    } else {
      character(0)
    }
    results[[sample_id]] <- list(
      n_input = as.integer(nrow(ff)),
      n_removed = as.integer(sum(!good)),
      ## as.list() so a single TRUE still serialises as a JSON array; with
      ## auto_unbox = TRUE a length-one logical vector collapses to a bare scalar and a
      ## reader that iterates it would iterate nothing at all.
      good_cells = as.list(as.logical(good)),
      percentage_removed = as_opt(res$PercentageRemoved),
      it_fraction = it_fraction,
      mad_fraction = mad_fraction,
      it_available = isTRUE(it_available),
      it_reason = as_opt(it_reason),
      result_keys = as.list(result_keys),
      source_file = f,
      report_path = report_dir,
      ## Reported as counted, not as requested: PeacoQC writes no report when
      ## determine_good_cells = "MAD", so report = TRUE is not a promise it keeps on
      ## every route. Asserting that here is what keeps "report = TRUE" from becoming a
      ## setting that is accepted and quietly does nothing.
      output_files = as.list(written),
      output_file_count = length(written),
      warnings = as.list(run$warnings)
    )
  }

  batchdx_write_result(parsed$result_path, list(
    status = "ok",
    samples = results,
    criteria = criteria,
    channels = as.list(channels),
    peacoqc_version = as.character(utils::packageVersion("PeacoQC")),
    flowcore_version = as.character(utils::packageVersion("flowCore")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    ## unname(): lapply over the named `results` list attaches the sample_id as a
    ## name, and a named character vector serialises as a JSON object rather than
    ## the array of strings the Python validator requires.
    warnings = unname(c(warnings_all, unlist(lapply(results, function(r) as.character(r$warnings))))),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  ))
}

main()
