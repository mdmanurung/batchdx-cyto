# 0005 — Evaluation partitions are frozen once, before any candidate runs

Status: accepted (Stage 0)

## Context

If a correction is allowed to re-derive the labels it is scored against, it is scored
against a partition it created. cyCombine's own benchmarking vignette describes exactly
this: correction-specific clustering changes the *uncorrected* EMD baseline.

## Decision

Two namespaces, kept apart:

* **native** — upstream-compatible mode. Upstream's own partitioning behaviour is
  preserved, including the moving uncorrected baseline. Results carry the `native`
  namespace and are never mixed with common numbers.
* **common** — one frozen partition shared by every candidate: the same evaluation
  cells, the same labels, the same marker subset, the same comparison pairs.

Evaluation labels come from one of three tiers, and the tier and its confidence are
carried on *every* metric record: user annotation (high), threshold gating in the
declared coordinate system (medium, the default when no annotation exists), or
clustering of the uncorrected baseline (low, flagged as potentially batch-driven and
never sufficient on its own to support a recommendation).

Quantitative neighbourhood metrics are pinned to four settings recorded on every
record: marker space, the evaluation marker mask, the declared per-marker scaling,
Euclidean distance, and a declared `k`. UMAP is visualisation and is never the scoring
coordinate system.

## Consequences

* Candidate order cannot change labels or comparisons (tested).
* A metric computed in native mode is not comparable with the same-named common metric,
  and the namespaces make that visible rather than implicit.
