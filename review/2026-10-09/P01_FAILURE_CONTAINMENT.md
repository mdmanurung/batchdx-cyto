# P01 — Candidate failure containment

Date: 2026-10-10. Implementation authorized by the user's “next” after the
prototype task checklist and explicit next-ticket handoff. Status: implemented and
verified for the candidate-stage acceptance checks; broader suite has one existing
failure reproduced on an untouched HEAD snapshot.

## Problem and change

The first regression used the real output validator to reject a candidate matrix
containing infinity. `stage_correct` let `OutputValidationError` escape, stopping
the remaining candidates and report. The pre-fix test failed with that exception;
the receipt preserves its JSON report.

The correction boundary now contains expected backend, value/JSON, file/I/O and
subprocess-timeout errors, including H5AD export validation. It writes explicit
`failed` or typed `not_available` results and retains tracebacks in the existing
candidate workflow log. An availability-probe error preserves the original reason
and becomes a separate warning. Shared loading/config/design errors remain outside
the protected block; cancellation and unexpected programming defects still raise.

Before an attempt, the candidate's old correction result, corrected H5AD and metric
result are invalidated. A partially written export is removed after failure. A
returned non-success result does not export a matrix. Configuration validates
backend names and candidate directory IDs before any input manifests are written.

## Reach and validation

- Source: `src/batchdx/cyto/pipeline.py` and `config.py`; the shared stage is called
  by CLI and Snakemake. No backend algorithm or metric formula changed.
- New stage regressions: invalid output, missing files, malformed JSON, timeouts,
  typed availability, stale retry, partial export, export validation, software-probe
  failure, non-success return, cancellation, shared-input failure and invalid config.
- The public stage sequence exercises identity, an offset-removal test fixture and
  a broken candidate through scoring, decision, status and offline report.
- One Snakemake regression uses the actual cyCombine adapter with an interpreter
  fixture that answers version probes but omits output. This validates failure
  propagation through the DAG; it is not R execution or real correction evidence.
- Affected existing unit/scientific, backend and report checks plus scoped final
  config/candidate checks are recorded separately from the full package suite.

## Boundaries and continuation

The [verification receipt](prototype-evidence/2026-10-10-p01-verification.json)
links all six test reports and hashes the source. Final configuration/candidate
checks: 39 passed, including 23 new stage/config cases. The new Snakemake check:
1 passed, 14 deselected. Related unit/scientific/backend/report checks: 307 passed,
1 skipped, 1 failed, 5 deselected. The sole failure is
`test_cytonorm_is_registered_as_a_fit_apply_backend`: CytoNorm is unavailable on
the test PATH while its unchanged capability declares `fit_apply_verified=True`.
The identical assertion fails on a clean HEAD source snapshot under `/tmp`.
The related suite started before the final ID guard; the final scoped config and
candidate suite verified that addition. Ruff and whitespace checks passed.

The P01 acceptance checks pass; the broader validation run is not green. P02
still needs strict R result-payload validation, list cardinality and declared
timeouts. P03 still needs broader DAG input/output freshness; this change does not
declare those requirements complete. The full nine-score metric suite and ranked
benchmark visuals remain required and unimplemented by this ticket.

No dependency installation, real R/backend validation, measured-data correction,
full integration/interrupt suite, build, remote CI, scheduler work, commit, push or
publication is included. Existing mixed worktree changes are preserved.
