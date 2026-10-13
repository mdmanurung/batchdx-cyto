# P02 — Checked R result boundary

Date: 2026-10-10. Implementation authorized by the user's “next” and “resume”.
Status: implemented; scoped acceptance checks passed. Broader package validation
remains non-green. See the [source-hashed receipt](prototype-evidence/2026-10-10-p02-verification.json).

## Change

A successful subprocess exit now requires a readable JSON object in `result.json`
with `status` equal to `ok` or `not_available`. Missing files, malformed JSON,
non-object envelopes, unsupported status values, non-JSON NaN/Infinity constants
and malformed declared string arrays raise the existing `RScriptFailed` exception.
The captured stdout/stderr and validation reason are persisted to the existing
log. Old result JSON is removed before execution, as before. Domain matrices,
tables and per-sample results retain their existing caller-specific checks.
cyCombine additionally refuses a non-success envelope before reading any matrix.

The shared R writer preserves arrays for warnings, SOM node IDs, fitted/unfitted
batch labels, missing packages, metaclusters and output paths, including zero and
one item. This fixes the root serialization behavior across all six R scripts.
Python checks these declared string arrays rather than iterating scalar characters.

Execution defaults to 14,400 seconds, matching the correction rule's declared
240-minute budget. Existing explicit positive finite timeouts remain supported;
`None` now selects the bounded default. Version probes have a 30-second limit and
record a timeout as unavailable software with an explicit reason. Execution
timeouts retain partial stdout/stderr and timeout duration in the log and raise
`RScriptFailed` with `returncode=None`; no exit status is invented. All six callers
report the exception's reason instead of describing every failure as a non-zero
exit. No algorithm, new dependency or configuration option was added.

## Verification and boundaries

Six initial payload regressions failed before the change. Focused tests exercise
real interpreter-fixture subprocesses, stale result removal, list cardinality,
non-zero exit, partial timeout logs, bounded defaults, invalid timeouts and
explicit unavailability. Correction, detection and metric consumers exercise
the shared subprocess boundary with missing/malformed JSON and timeouts; failures
produce errors/failed records with no metric value or detector assessments.
These interpreter fixtures do not execute correction algorithms.

The R exchange suite uses R4_51 (R 4.5.1) and the installed `jsonlite`/`testthat`,
with actual R-written empty/single/multiple warnings and node IDs parsed by Python.
The first isolated new R test had an unavailable helper function; it was corrected
to use the existing shared helper, then the full exchange suite passed.
`testthat` prints a startup warning that it was built under R 4.5.2; that warning
is retained in the receipt separately from the suite's assertion summary.

Related checks on the default PATH expose three CytoNorm failures: the verified
capability assertion and two guards expecting split/reference errors before an
unavailable package. All three reproduce on an untouched HEAD archive; they are
not waived, and these checks do not establish a green full package validation.
The initial workflow attempt skipped because Snakemake was not on PATH; it was
rerun with the existing project executable. Its simulated interpreter now writes
a valid result envelope while still omitting the matrix, preserving the P01 test's
purpose under the stricter boundary.

Final focused Python checks: **40 passed** (30 boundary and 10 caller cases).
Full R exchange suite: **110 assertions passed**, no failures/errors/skips.
Real Snakemake missing-matrix regression: **1 passed**, 14 deselected, 414 seconds.
Related Python run: **93 passed, 3 failed**, 25 deselected; it preceded the final
non-JSON-number and correction-status guards, which the final focused suite covers.
Baseline reproduction: the same **3 failures**, 14 deselected. Ruff, whitespace,
document-link and strict Mycelium structure checks passed. The raw Python reports
are copied beside the receipt; R counts and its exact command are recorded from
observed tool output, without a separate raw R console file.

P03 workflow input/output freshness remains next. Full nine-score OpenProblems
integration and ranked benchmark visuals remain required. No real backend or
measured-data validation, full Python package gate, install/build, remote CI,
scheduler job, commit, push or publication is included.
