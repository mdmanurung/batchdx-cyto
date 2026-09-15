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
