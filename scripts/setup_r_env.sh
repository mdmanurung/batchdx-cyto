#!/usr/bin/env bash
# Build the isolated R environment that the batchdx-cyto R backends run against.
#
# Why an isolated environment rather than a personal R library: decision
# docs/decisions/0004-dependency-isolation.md requires the R dependencies of a run
# to be attributable to this package.  Installing into a shared user library (which
# on many machines already holds hundreds of unrelated packages) makes a run's R
# environment unstateable, because nothing records which of those packages the run
# actually needed.  The environment built here is self-contained and reproducible
# from a pinned list.
#
# This is the local-environment counterpart of containers/Dockerfile.r-backends.
# The container remains the answer for reproducible deployment; this is what lets
# the R adapters be executed and tested on a machine with no container runtime.
# Neither mechanism substitutes for the other.
#
# WHY A CONDA ENVIRONMENT AND NOT install.packages() INTO AN R LIBRARY
# ---------------------------------------------------------------
# Installing R packages from CRAN source on this machine fails, and the failure is
# a toolchain mismatch rather than anything about the packages:
#
#   * R 4.4.1 requires GLIBCXX_3.4.29; the default gcc here is 8.5.0, which
#     provides at most GLIBCXX_3.4.25.
#   * A conda-provided R 4.4.3 works, but CRAN's current Rcpp source requires
#     CXXABI_1.3.15, which only gcc 14 provides; the newest gcc on this node is
#     13.2.0 (CXXABI_1.3.14).
#   * Posit Package Manager has no prebuilt binary for centos8/R 4.4.
#
# A conda environment resolves the whole closure to prebuilt packages built
# against one consistent toolchain, so none of that applies. cyCombine and CytoNorm
# are installed from pinned local checkouts, without GitHub API credentials.
#
# Usage:
#   scripts/fetch_upstream.sh .upstream  # required pinned source checkouts
#   scripts/setup_r_env.sh            # build ./renv (the default prefix)
#   scripts/setup_r_env.sh <prefix>   # build elsewhere
#
# Re-running is safe: an existing environment at the prefix is reused and the
# source packages are only installed when absent.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${1:-$REPO_ROOT/renv}"
LOCK="$REPO_ROOT/config/upstream_lock.yaml"
MAMBA_BIN="${MAMBA_BIN:-micromamba}"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-$REPO_ROOT/.micromamba}"

# cyCombine's own commit is read from the lock rather than repeated here, so there
# is exactly one place a commit can be wrong.

CYCOMBINE_COMMIT="$(awk '
  /^  cycombine:/ {inside=1; next}
  inside && /^  [a-zA-Z]/ {inside=0}
  inside && /commit:/ {gsub(/[" ]/, "", $2); print $2; exit}
' "$LOCK")"
if [[ -z "$CYCOMBINE_COMMIT" ]]; then
  echo "could not read the cyCombine commit from $LOCK" >&2
  exit 2
fi

# Every Rscript invocation below passes --vanilla.  That is load-bearing, not
# stylistic.  A ~/.Renviron on many machines (including this one) sets
# R_LIBS_USER to a personal renv library, and Renviron entries override the
# environment variable, so without --vanilla an install silently writes into that
# personal library instead of into $PREFIX.  --vanilla implies --no-environ, so
# R_LIBS_USER from the environment is honoured.  The adapters already run under
# --vanilla (see batchdx.cyto.adapters.r_runner.run_r_script), which makes
# setup-time and run-time library resolution identical.

if ! command -v "$MAMBA_BIN" >/dev/null 2>&1; then
  cat >&2 <<'EOF'
error: micromamba not found on PATH.

This script builds the R environment with micromamba.  Fetch it with:

  curl -sSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp bin/micromamba
  export MAMBA_BIN=/tmp/bin/micromamba

Setting MAMBA_BIN to an existing micromamba is the supported override.
EOF
  exit 2
fi

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"

if [[ -d "$PREFIX" ]] && "$MAMBA_BIN" run -p "$PREFIX" Rscript --vanilla -e 'quit(status=0)' 2>/dev/null; then
  echo "== reusing existing environment at $PREFIX"
else
  echo "== creating R environment at $PREFIX"
  echo "   (r-base is pinned to 4.4 so the R version a run sees is a decided fact,"
  echo "    not whatever the machine happened to have)"
  # The list below is not cyCombine's DESCRIPTION Imports verbatim.  Three additions
  # are needed, and each was found by running the code rather than by reading the
  # manifest, so the reasons are recorded here:
  #
  #   pbmcapply  compute_emd() calls cyCombine:::check_package("pbmcapply") on its
  #              multicore path and aborts rather than degrading.  It is a Suggest,
  #              so dependency resolution does not install it, and without it the
  #              native detector fails.
  #   emdist     compute_emd() calls emdist::emd2d() unconditionally, so a Suggest
  #              is a hard requirement in practice.
  #   outliers   the same for the MAD evaluation path.
  #
  #   peacoqc    workflow/scripts/peacoqc_qc.R calls PeacoQC::PeacoQC, the
  #              event-level QC step of the preprocessing stage.  Pinned to 1.16.0
  #              because that is what the pinned SPiSCy reference resolves.
  #   flowvs     NOT called by any script.  It is here only so that
  #              tests/integration/test_peacoqc_backend.py can execute
  #              flowVS::transFlowVS and check batchdx's asinh against it, rather
  #              than against a reading of its documentation.
  #
  # The rest of cyCombine's Suggests are genuinely optional here: the plotting and
  # extra-format suggestions are not reached by the detection, correction or metric
  # routes this package uses.
  "$MAMBA_BIN" create -y -p "$PREFIX" -c conda-forge -c bioconda \
    'r-base=4.4' \
    bioconductor-sva bioconductor-flowcore bioconductor-flowsom \
    bioconductor-peacoqc bioconductor-flowvs \
    r-kohonen r-emdist r-outliers r-remotes r-pbmcapply \
    r-tibble r-readxl r-readr r-stringr r-dplyr r-magrittr r-jsonlite r-matrixstats r-testthat \
    r-anndata \
    r-ggridges r-ggplot2 r-cowplot r-matrix r-rcolorbrewer r-viridis r-plyr r-tidyr \
    r-uwot r-pheatmap r-ggpubr r-gridextra
fi

# CytoBatchFlagR's own dependency set, installed separately and idempotently.
#
# It is a separate step, not more lines in the create above, for two reasons.  The create
# only runs when the prefix is absent, so packages added later would never reach an
# environment that already existed; and upstream has no DESCRIPTION and no NAMESPACE, so
# the list cannot be derived from a manifest and has to be stated here explicitly.
#
# The binding set is the 29 UNCONDITIONAL library() calls at
# .upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R:6-34.  `library()` is not wrapped
# in require(), so a single absent package aborts the entire source() rather than
# degrading one function -- the list is all-or-nothing, and trimming it to work around a
# missing package would be the derivative work docs/decisions/0003 forbids.
#
# Of those 29, six were absent from this environment on 2026-09-30 and are installed
# here.  Channel placement was checked against the live indices, not from memory:
# bioconductor-complexheatmap is on bioconda; the other five are on conda-forge.  Each has
# an r44 build, so the pinned R 4.4 is satisfied.  r-umap pulls r-rcppannoy and r-uwot
# (uwot was already present); r-tidyverse pulls r-forcats and r-broom; and
# bioconductor-complexheatmap pulls r-circlize and r-clue.
#
# NOT installed, deliberately: randomcoloR, DT, shinyjs, shinycssloaders and shinyWidgets.
# Upstream's README and Rmd name them, but they are reached only through
# R/CytoBatchFlagR_workflow.R and the Shiny helpers, which this adapter deliberately never
# sources (see adapter_strategy in config/upstream_lock.yaml and ADR 0003).  Installing
# them would widen the dependency surface for code paths that cannot run here.  The
# verification block below asserts the absence of that reasoning being wrong: it derives
# the list from the pinned file at run time rather than from this comment.
echo "== ensuring the CytoBatchFlagR library() block can be satisfied"
"$MAMBA_BIN" install -y -p "$PREFIX" -c conda-forge -c bioconda \
  bioconductor-complexheatmap r-circlize r-laplacesdemon r-patchwork r-tidyverse r-umap r-testthat

echo "== R interpreter"
"$MAMBA_BIN" run -p "$PREFIX" Rscript --vanilla -e 'cat(R.version.string, "\n"); cat("libPaths:", .libPaths(), "\n")'

# CytoNorm's own commit, read from the lock for the same reason cyCombine's is.
CYTONORM_COMMIT="$(awk '
  /^  cytonorm:/ {inside=1; next}
  inside && /^  [a-zA-Z]/ {inside=0}
  inside && /commit:/ {gsub(/[" ]/, "", $2); print $2; exit}
' "$LOCK")"
if [[ -z "$CYTONORM_COMMIT" ]]; then
  echo "could not read the CytoNorm commit from $LOCK" >&2
  exit 2
fi

verify_checkout() {
  local name="$1" expected="$2" dir="$REPO_ROOT/.upstream/$1" observed
  if [[ ! -d "$dir/.git" ]]; then
    echo "missing pinned checkout $dir; run bash scripts/fetch_upstream.sh .upstream" >&2
    exit 2
  fi
  observed="$(git -C "$dir" rev-parse HEAD)"
  if [[ "$observed" != "$expected" ]]; then
    echo "$name checkout is $observed; expected $expected; rerun scripts/fetch_upstream.sh" >&2
    exit 2
  fi
}

verify_checkout cyCombine "$CYCOMBINE_COMMIT"
verify_checkout CytoNorm "$CYTONORM_COMMIT"

echo "== installing cyCombine from pinned local checkout at $CYCOMBINE_COMMIT"
if ! "$MAMBA_BIN" run -p "$PREFIX" Rscript --vanilla -e 'quit(status=if (requireNamespace("cyCombine", quietly=TRUE)) 0 else 1)'; then
  "$MAMBA_BIN" run -p "$PREFIX" R CMD INSTALL --no-multiarch "$REPO_ROOT/.upstream/cyCombine"
fi

echo "== installing CytoNorm from pinned local checkout at $CYTONORM_COMMIT"
# CytoNorm is GPL (>= 2).  It is installed and called, never copied into this package;
# see THIRD_PARTY_NOTICES.md and docs/decisions/0003-cytobatchflagr-integration.md for
# the arm's-length rule that applies to it. It is installed from the checkout,
# never copied into this package.
if ! "$MAMBA_BIN" run -p "$PREFIX" Rscript --vanilla -e 'quit(status=if (requireNamespace("CytoNorm", quietly=TRUE)) 0 else 1)'; then
  "$MAMBA_BIN" run -p "$PREFIX" R CMD INSTALL --no-multiarch "$REPO_ROOT/.upstream/CytoNorm"
fi

echo "== verification"
# The verification below resolves the pinned CytoBatchFlagR checkout relative to the
# repository rather than to the caller's working directory, so the result does not
# depend on where this script was invoked from.
BATCHDX_UPSTREAM_DIR="$REPO_ROOT/.upstream" \
"$MAMBA_BIN" run -p "$PREFIX" Rscript --vanilla -e '
  suppressPackageStartupMessages(library(cyCombine))
  libs <- c("cyCombine", "CytoNorm", "sva", "kohonen", "emdist", "outliers",
            "flowCore", "FlowSOM", "PeacoQC", "flowVS",
            "magrittr", "tibble", "dplyr", "readxl", "readr", "stringr", "jsonlite", "testthat")
  ok <- TRUE
  for (p in libs) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    if (is.na(v)) { ok <- FALSE; cat(sprintf("   MISSING  %s\n", p)) }
    else cat(sprintf("   ok       %-12s %s\n", p, v))
  }
  if (!ok) {
    cat("\nA required R package is absent.  The R-backed adapters will report\n")
    cat("not_available with a reason; they will not substitute another algorithm.\n")
    quit(status = 1)
  }
  cat("\ncyCombine location:", find.package("cyCombine"), "\n")
  cat("cyCombine exports the functions the adapters call:\n")
  need <- c("normalize", "create_som", "correct_data", "compute_emd",
            "detect_batch_effect", "detect_batch_effect_express",
            "evaluate_emd", "evaluate_mad", "get_markers")
  for (f in need) {
    present <- f %in% getNamespaceExports("cyCombine")
    cat(sprintf("   %-28s %s\n", f, if (present) "ok" else "ABSENT"))
    if (!present) quit(status = 1)
  }

  # CytoBatchFlagR, checked the only way it can be: against the library() block of the
  # file that will actually be sourced.  Upstream has no DESCRIPTION and no NAMESPACE, so
  # a list written down anywhere else -- including the comment above -- can drift from
  # the pinned checkout.  Deriving it here makes this a statement about the real binding
  # set, and it is all-or-nothing: one absence aborts the whole source().
  cbf <- file.path(Sys.getenv("BATCHDX_UPSTREAM_DIR", ".upstream"),
                   "CytoBatchFlagR", "R", "CytoBatchFlagR_functions.R")
  if (!file.exists(cbf)) {
    cat(sprintf("\nCytoBatchFlagR checkout absent at %s; its library() block cannot be\n", cbf))
    cat("verified here.  The adapter reports not_available in that case, which is correct,\n")
    cat("and the R-backed detector is simply not exercised.  Run scripts/fetch_upstream.sh\n")
    cat("to fetch the pinned commit and re-run this script.\n")
  } else {
    lines <- readLines(cbf, warn = FALSE)
    hits <- regmatches(lines, regexec("library[(][[:space:]]*([A-Za-z][A-Za-z0-9._]*)", lines))
    pkgs <- sort(unique(unlist(lapply(hits, function(m) if (length(m) > 1L) m[[2L]] else NULL))))
    cat(sprintf("\nCytoBatchFlagR library() block: %d packages\n", length(pkgs)))
    absent <- character(0)
    for (p in pkgs) {
      if (requireNamespace(p, quietly = TRUE)) {
        cat(sprintf("   ok       %-20s %s\n", p,
                    as.character(packageVersion(p))))
      } else {
        absent <- c(absent, p)
        cat(sprintf("   MISSING  %-20s\n", p))
      }
    }
    if (length(absent) > 0L) {
      cat(sprintf("\n%d package(s) of CytoBatchFlagR unconditional library() block are absent:\n",
                  length(absent)))
      cat(paste0("  ", paste(absent, collapse = ", "), "\n"))
      cat("library() is not wrapped in require(), so source() aborts entirely.  The R\n")
      cat("layer will report not_available; it will not fabricate a batch-effect verdict.\n")
      quit(status = 1)
    }
    # The four functions the adapter calls, confirmed present after a real source().
    ok_source <- tryCatch({
      suppressWarnings(suppressMessages(source(cbf, local = FALSE)))
      TRUE
    }, error = function(e) { cat("   source() FAILED:", conditionMessage(e), "\n"); FALSE })
    if (!ok_source) quit(status = 1)
    need <- c("neg_mfi_df_all", "pos_mfi_df_all", "percent_pos_df_all", "iqr_dataframe")
    for (f in need) {
      present <- exists(f, mode = "function")
      cat(sprintf("   %-20s %s\n", f, if (present) "ok" else "ABSENT"))
      if (!present) quit(status = 1)
    }
  }
'

cat <<EOF

R environment ready at: $PREFIX

Use it by putting its bin/ first on PATH in every shell that runs the R layer:

  export PATH="$PREFIX/bin:\$PATH"

To run an R script inside it directly:

  $MAMBA_BIN run -p $PREFIX Rscript --vanilla <script.R>

Nothing was written outside $PREFIX.
EOF
