# 0003 — CytoBatchFlagR is invoked at arm's length, never vendored

Status: accepted (Stage 0). **This decision is a licence decision and was recorded
before the adapter was written.**

## Context

CytoBatchFlagR (published 2026 in Cytometry Part A; first released as `cytoFlagR`,
bioRxiv 2025) is distributed as **R scripts**, not as an installable package: the
checkout at the pinned commit has no `DESCRIPTION` and no `NAMESPACE`. Its entry point
`R/CytoBatchFlagR_workflow.R` is interactive — it calls `rstudioapi::selectDirectory`
and several Shiny helpers.

The repository carries **GPL-3.0**. batchdx-cyto is MIT.

## Decision

Integrate at arm's length:

1. The user supplies a checkout at the commit pinned in `config/upstream_lock.yaml`
   (`scripts/fetch_upstream.sh` fetches it).
2. `workflow/scripts/cytobatchflagr_detect.R` `source()`s `R/CytoBatchFlagR_functions.R`
   **unmodified**, in a separate R process.
3. The interactive workflow script is deliberately not sourced. The analysis functions
   (`neg_mfi_df_all`, `pos_mfi_df_all`, `percent_pos_df_all`, `iqr_dataframe`) are
   called directly with values taken from configuration: directories, marker
   selections, reference identities and thresholds.
4. **No upstream source is copied into this package and none is refactored.**

Refactoring upstream into a package — the obvious engineering improvement — would
create a derivative work and is therefore out of scope for an MIT distribution. If that
becomes necessary, the correct route is a separate GPL-3 package, not a change here.

## Consequences

* The adapter reports `not_available` with an explicit reason when no checkout is
  present, rather than failing obscurely or fabricating flags.
* The three IQR metrics require a control sample in every batch; that constraint, and
  the arithmetic fact that the 1.5 IQR fence cannot fire for n <= 3 batches, are
  encoded in the Stage 3 eligibility manifest with their origin stated.
* Thresholds must be supplied. The publication reports five experienced analysts
  choosing five different thresholds for the same markers, so the adapter refuses to
  pick them itself.
