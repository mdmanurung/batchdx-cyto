# 0001 — The internal data contract is AnnData plus seven typed records

Status: accepted (Stage 0)

## Context

Cytometry arrives as FCS or as AnnData/H5AD, and every downstream stage needs to agree
on what a cell is, what a marker is, and what has already been done to the values.

## Decision

AnnData is the canonical *logical* interchange contract. The expression values live in
a named layer that the configuration declares; `X` is never used implicitly.

Five manifests describe the data (cell index, sample, marker, processing, candidate),
and seven versioned records cross every stage boundary: `DatasetManifest`,
`DesignAudit`, `SamplingManifest`, `DetectorResult`, `CorrectionResult`,
`MetricResult`, `DecisionResult`. They are pydantic models with `extra="forbid"` and a
`schema_version`, and their JSON Schemas are exported for the R side
(`schemas/contracts.py::export_json_schemas`).

The preprocessed baseline is immutable. Each correction candidate writes its own
AnnData file; nothing accumulates candidates into one object.

## Consequences

* A record written by a newer build fails loudly instead of being half-read.
* Adding a field is a contract change and bumps `CONTRACT_SCHEMA_VERSION`.
* Memory cost: one baseline plus one candidate matrix at a time, not N candidates.
