# 0004 — Upstream R packages are pinned by commit and need a container

Status: accepted (Stage 0)

## Context

cyCombine and CytoNorm are distributed through GitHub, not CRAN or Bioconductor.
cyCombine imports `sva` (Bioconductor) and `kohonen` (CRAN); CytoNorm imports
`flowCore` and `FlowSOM` (Bioconductor). A conda environment file alone does not pin a
GitHub-sourced R package reproducibly.

## Decision

* `config/upstream_lock.yaml` pins a full 40-character commit for every resource, plus
  its licence, its adapter strategy, and the exact list of functions our adapters call.
* `scripts/fetch_upstream.sh` materialises those commits.
* Reproducible *execution* of the R backends requires the container definition in
  `containers/`, because the transitive CRAN and Bioconductor dependencies cannot be
  pinned by commit from R alone.
* Runner dependencies (Python, Snakemake) are separate from task dependencies (R and
  its packages).

## Consequences

* Without the container, `cycombine` reports `not_available` with the reason, and the
  Stage 6 gate for it is recorded as blocked rather than passed.
* The audit scripts run in two modes: against an installed library when one is present,
  and against the pinned source otherwise. Both *read* a signature; neither infers one.
