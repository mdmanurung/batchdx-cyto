# 0006 — Success, ineligible and failed are three different things

Status: accepted (Stage 0)

## Context

The cheapest way to make a pipeline look complete is to let a missing dependency, an
invalid comparison or a crashed method fall through as "no result", which a reader then
interprets as "no problem".

## Decision

* `CorrectionStatus` distinguishes `success`, `ineligible`, `failed` and
  `not_available`; `MetricStatus` distinguishes `ok`, `ineligible`,
  `insufficient_data`, `undefined`, `failed` and `not_implemented`; `DecisionState`
  is separate again from both.
* A failed backend never returns the unchanged input. The identity candidate is the
  only thing allowed to return the baseline, and it is a genuine candidate that goes
  through the same validation, export and scoring path as any other.
* Missing values are never zero. A near-zero denominator produces `undefined`, not a
  score stabilised with an epsilon.
* Downstream components read status from records, never from the existence of a
  directory. `stage_status` builds the roster from the validated configuration and
  reports one entry per planned unit.
* Relative reductions use the same eligible comparisons in numerator and denominator,
  so dropping an unavailable comparison cannot improve a candidate's standing.

## Consequences

* The report has a "failures and limitations" section that is populated from records,
  and failed candidates stay visible in the comparison table.
* A partial run is reportable as a partial run.
* R exit zero requires a valid result envelope. Missing/malformed results and
  timeouts raise the checked boundary's existing failure type with logs; callers
  report the actual reason. Timeouts have no reported exit code. R execution uses
  a four-hour default and version probes a 30-second limit. Declared string arrays
  retain their list shape across the R/JSON/Python boundary.
* Expected candidate execution/output errors (including invalid values, missing
  files, malformed JSON, timeouts and export I/O failures) produce a failure record
  and a traceback in `logs/correct_<candidate>.log`. A failed software-version
  probe cannot replace the original error. Invalid shared inputs, cancellation,
  unexpected programming defects and inability to persist the failure remain
  whole-run errors.
* Before a correction retry, its old result, corrected H5AD and metric result are
  invalidated. Failed/partial exports are removed; only a successful correction
  exports a matrix. The old file cannot supply scores after a new attempt fails.
* Backend names are checked during configuration validation. Candidate IDs must
  be non-empty directory names without path separators, so retry invalidation
  cannot traverse into another candidate's directory.
