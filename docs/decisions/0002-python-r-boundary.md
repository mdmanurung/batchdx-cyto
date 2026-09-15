# 0002 — R runs in subprocesses, never embedded in the Python process

Status: accepted (Stage 0)

## Context

cyCombine, CytoNorm and CytoBatchFlagR are R. The obvious alternative is an embedded
interpreter (rpy2), which puts a second garbage-collected runtime inside the analysis
process and turns an upstream segfault into a lost run.

## Decision

Every R call is `Rscript --vanilla workflow/scripts/<name>.R params.json result.json`
in a subprocess. Inputs and outputs cross as JSON plus a typed binary matrix
(`io/binary.py` and `workflow/scripts/lib/binmat.R`: raw little-endian float64 with a
JSON sidecar carrying shape, row ids and column names).

The binary format uses base R's `readBin` only, so the exchange adds no R package
dependency of its own beyond `jsonlite`.

## Consequences

* An R failure is a non-zero exit status plus a captured log, never a corrupted
  in-process state; `RScriptFailed` carries the log into the `CorrectionResult`.
* No large delimited-text intermediates.
* Cost: one process spawn and one matrix round trip per call. Acceptable, because the
  R side is doing minutes of work, not milliseconds.
