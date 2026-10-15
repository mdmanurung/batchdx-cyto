# batchdx-cyto — remaining work, tranche by tranche

Date: 2026-09-30
Last updated: 2026-09-30 (after Tranche 1)
Project: /exports/para-lipg-hpc/mdmanurung/batchdx-cyto
Status: in progress — Tranche 1 complete

## Objective

Close every remaining gap in the package in verifiable tranches: repair the report so
nothing is silently dropped, execute the one R backend that does not run, implement
preprocessing, make the deployment artifacts real, and reach a defensible release state.

Success for this plan is a suite that passes with no gate BLOCKED, and a status document
whose every line corresponds to a test that actually ran. It is not "feature complete".

## The correction that reopens Tranche 1

The previous plan recorded Step 2.1 as **delivered**, with six specific fixes claimed as
evidence. All six were verified against the code on 2026-09-30 and **all six are absent**:

| Recorded as delivered | Verified actual state (2026-09-30) |
|---|---|
| `pipeline.py:236` bare `pass` replaced | still there, now at `src/batchdx/cyto/pipeline.py:242` |
| `DecisionResult.overrides` given a producer | no caller of `overrides=` in `pipeline.py` |
| `report.html.j2:200` native-detector filter removed | still `{% if d.detector in ('acquisition_drift', 'composition') %}` |
| truncation caps declared rather than silent | `html.py:43` `MAX_EXPLORER_MARKERS = 8`, `:168` `markers[:3]`, both unchanged |
| `n_donors` on all metric families | one family only, `metrics/common.py:105` |
| `Redactor.path()` wired in | no callers anywhere in the repository |

The claim was written; the work was not done. The four parallel agents dispatched for the
prior tranche were scoped to metrics, splits, CytoNorm and OpenProblems, and Step 2.1 was
never assigned to one of them.

One of the six inverts the direction of the error: `docs/usage.md:234-238` says "nothing
calls `Redactor.path()`", and that documentation is **correct**. The plan was the wrong
document, not the code. `docs/usage.md` is therefore not a defect for this item.

## Second structural finding: 13 untracked files

`git status` shows 13 untracked paths, including `scripts/setup_r_env.sh` — the script
that builds the entire R environment, referenced by eleven tracked files — and
`workflow/scripts/cycombine_evaluate.R` and `tests/unit/test_config_keys.py`, which a
review agent initially reported as absent because they are in no commit. All work from
the previous tranche is therefore uncommitted and unbuildable by a fresh clone. Tranche 5
commits them.

## Progress Summary

- [x] Tranche 1 — report correctness and release blockers — **complete, 292 passed**
- [ ] Tranche 2 — CytoBatchFlagR happy path
- [ ] Tranche 3 — SPiSCy-style preprocessing
- [ ] Tranche 4 — workflow, CI, deployment
- [ ] Tranche 5 — release readiness

## Review Findings

| Finding | Resolution |
|---|---|
| The prior plan's Step 2.1 evidence was false | Step 2.1 is reopened as Tranche 1 items 1.1–1.8, and the false claim is recorded above rather than deleted |
| `docs/usage.md` §8 was suspected stale | Verified accurate; the plan was wrong. §8 is corrected only for the *example* it cites (see 1.5) |
| `workflow/profiles/slurm/config.yaml:28` is a dangling comment, not `Snakefile:28` | The earlier plan misattributed it to the Snakefile. Corrected in 4.7 |
| Seed 473 is in 8 places, not 5 | Verified 8 sites; corrected in 5.2 |
| `cycombine_evaluate.R` and `test_config_keys.py` reported "absent" | They exist on disk and are untracked. Recorded as a finding, not a defect |
| `docs/blocked.md:301-302` says three preprocessing keys are in `DECLARATIVE_ONLY` | They are in `config.UNSUPPORTED_KEYS`; corrected in 3.5 |
| Parallel agents in one checkout corrupted each other's test results (previous tranche) | Tranches run sequentially; the test suite is the only arbiter of success |

## Implementation Steps

### [x] Tranche 1 — report correctness and release blockers (complete)

Pure correctness and metadata. No R, no network, no external data. Run this tranche first:
a report that silently omits a successful native detector is the same failure class this
package exists to prevent, so it outranks new capability work.

#### [x] 1.1 The report stops dropping successful native detectors

`report.html.j2:200` gates the diagnostics section on
`{% if d.detector in ('acquisition_drift', 'composition') %}`. All five detector results
are fetched (`html.py:106-109`), put in context, iterated at `:199`, and three are then
skipped: `cycombine_express`, `cycombine_cluster`, `cytobatchflagr`. A **successful**
native detector currently appears nowhere; only a non-`ok` one surfaces, via the failures
table (`html.py:202-205`).

`test_stage10_report.py:189-190` asserts only the *skipped* path, which is why this went
unseen.

Fix: render every detector. Also render the fields that are fetched and never used —
`DetectorResult` (`contracts.py:292-309`) carries `tables`, `plots`, `eligibility`,
`requirements`, `parameters`, `software`, `detector_version`, and the template references
only `detector`, `status`, `reason`, `warnings`, `assessments`. So even the two rendered
detectors drop nine fields.

Validation: a test that runs the native detector successfully and asserts its name and
status appear in the HTML.

#### [x] 1.2 No assessment is silently discarded inside the diagnostics table

`report.html.j2:204` does `selectattr('flag')`, which drops every assessment whose `flag`
is `None`. That silently discards:

- `max_block_step_over_mad` (`diagnostics/drift.py:108-114`) — the module's own headline
  discontinuity statistic, emitted without a `flag` argument;
- `population_frequency` (`diagnostics/composition.py:67-74`) — every per-sample,
  per-population frequency row.

Render unflagged assessments in a separate labelled table. Same for the two silent
truncations: `flagged[:250]` (`:209`) and `joint[:200]` (`:308`) must state "showing N of
M".

#### [x] 1.3 The two hard-coded caps are declared

`html.py:43` `MAX_EXPLORER_MARKERS = 8` (applied `:141`) and `html.py:168` `markers[:3]`.
Either lift them or state the omitted count on the page. Surface
`sampling.evaluation_markers` length so a reader knows how many markers were not shown.

#### [x] 1.4 `_metric_rows` uses its `redactor` parameter, or loses it

`html.py:55` declares `redactor: Redactor` and the body (`:56-85`) never references it.
It is passed at `:134` and discarded. Meanwhile `MetricRecord.comparison_id` is projected
raw at `:66`, and `DetectorAssessment.sample_id` is projected raw at template `:210`.
Wire redaction through, or delete the parameter and say why it is not needed.

#### [x] 1.5 Decide `Redactor.path()` and fix the §8 example

Either call it, or delete it and simplify `usage.md` §8. Separately, §8's *example* is
wrong: it says the dataset manifest's provenance paths are printed, but the template
renders only `dataset_id` from that manifest. The real leak is `config | tojson(indent=2)`
(template `:342`), which carries `input.path` verbatim. Correct §8 to name the real one.

#### [x] 1.6 `n_donors` on every metric family

`metrics/common.py:105` sets it for EMD only, so the report's `n_donors` column
(`html.py:73`) is blank for MAD, frequency and joint-structure.

#### [x] 1.7 An exploratory override becomes visible

`pipeline.py:242` is the repository's only bare `pass`: with
`allow_exploratory_override` and an overridable ineligibility, the backend runs and
nothing is recorded. `DecisionResult.overrides` (`contracts.py:432`) has no producer,
because `stage_decide` never passes `overrides=`. Record the override and pass it through.

#### [x] 1.8 Stage 8 uncertainty reaches the report

`MetricRecord.uncertainty_low/high/kind` (`contracts.py:381-383`) are populated by every
preservation record and projected by nothing — not the table, not the CSV.

#### [x] 1.9 Add the MIT `LICENSE` file

`pyproject.toml:11` and `CITATION.cff:16` both declare MIT; **no licence text ships**.
Hard blocker for any distribution. Then switch to PEP 639
`license = "MIT"` + `license-files = ["LICENSE"]`.

#### [x] 1.10 Version and classifiers

Three files move together: `pyproject.toml:7`, `src/batchdx/cyto/_version.py:3`,
`CITATION.cff:15` — all `0.1.0.dev0`. Add `Programming Language :: Python :: 3.10/3.11/
3.12/3.13` and `:: R` (the package orchestrates five R scripts), and set a
`Development Status` that matches reality.

#### [x] 1.11 Correct the three stale tracked documents

The CI staleness gate (`.github/workflows/batchdx-cyto-tests.yaml:57-71`) covers only
`docs/implementation_status.md`, so these cannot drift silently:

- `THIRD_PARTY_NOTICES.md:30` — "no implementation was extracted" is false; six
  OpenProblems components are `compatibility` / `upstream_verified=True`. Its B-007
  pointer is also stale (B-007 is RESOLVED). Add the `metrics` extra's dependencies to
  the Python list at `:33-38`.
- `workflow/envs/README.md:3-4` — "One environment file lives here"; there are two.
- `config/upstream_lock.yaml:120-122` — `dependency_risk` names six packages; the enforced
  set is the 29 `library()` calls at `.upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R:6-34`.

#### [x] 1.12 The `WorkflowError` reference is now an explicit import

**The plan's premise here was wrong, and this item is corrected rather than executed as
written.** There was no `NameError`: Snakemake injects its own module globals into the
Snakefile namespace (`snakemake/workflow.py:234` `_globals = globals()`), so
`WorkflowError` was available. Two proofs: a Snakefile containing only
`raise WorkflowError(...)` works, and the existing test passes and prints the real
message. `test_correct_mode_without_a_candidate_fails_at_parse_time:136` **already
asserts the message**, not merely the exit code.

The explicit import was added anyway (`workflow/Snakefile:19`), labelled in place as a
dependency-instead-of-undocumented-injection change. Both `raise` branches were verified
live. Kept because an undeclared reliance on framework injection is a real fragility, not
because the original claim was true.

**Tranche 1 validation:** `pytest -q` green, `ruff check src tests scripts` clean, both
audit scripts exit 0, `python scripts/generate_status.py` regenerated and diff-clean
against the committed document.

---

### [ ] Tranche 2 — CytoBatchFlagR happy path (B-006)

#### [ ] 2.1 Install the 11 missing R packages

The binding dependency set is the 29 `library()` calls at
`.upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R:6-34`. **`library()` is
unconditional**, so one missing package aborts the entire `source()`. Missing from
`renv/lib/R/library`: `circlize`, `ComplexHeatmap`, `DT`, `LaplacesDemon`, `patchwork`,
`randomcoloR`, `tidyverse`, `umap`, `forcats`, `shinycssloaders`, `shinyjs`,
`shinyWidgets`. `ComplexHeatmap` is Bioconductor and pulls `circlize` + `DT`.

There is **no DESCRIPTION and no machine-readable dependency list anywhere** — the
`library()` block is the only authority, which is why ADR 0003's "unmodified" rule cannot
be worked around by installing a trimmed file.

#### [ ] 2.2 Extend the preflight from three packages to the full set

`cytobatchflagr_detect.R:40-44` checks only `dplyr`/`tidyr`/`matrixStats`. A missing
`circlize` currently produces a raw `source()` abort surfacing as `status="failed"`
(`cytobatchflagr.py:139-145`) rather than a stated `not_available`. This is original code,
not upstream, so ADR 0003 does not forbid the change.

#### [ ] 2.3 Guard the silent `<100`-cell drop

`neg_mfi_df_all` / `pos_mfi_df_all` default `min_cells=TRUE` and drop any control with
fewer than 100 cells for a marker (`functions.R:1041,1084,1090-1102`), `print()`ing a
summary. At `cells_per_sample=400` this can empty a batch — and a batch below four is
exactly the regime where the Tukey fence cannot fire, so a test would pass vacuously.
Either pass `min_cells=FALSE` and enforce coverage in the adapter, or assert per-batch
cell counts in the test.

#### [ ] 2.4 Write the happy-path test

It needs its own fixture with **≥4 batches and a control in every batch**
(`simulate.py:217` gives three, so the shared `frozen` fixture cannot be used), a real
`threshold_path`, and `@pytest.mark.requires_cytobatchflagr` — declared at
`pyproject.toml:72` and used by **zero** tests, which is itself the honest signal that the
path is unexercised. Both existing tests exit at the `min_batches=4` guard before reaching
their own assertions, so the threshold-refusal path (`cytobatchflagr.py:94-106`) has never
executed.

Include a negative control (no injected shift → no flags), or "flags something" is
indistinguishable from "flags anything". Assert the flagged batch is the shifted one, not
merely *a* batch. Assert the recorded commit equals `ad21e3b1` from
`config/upstream_lock.yaml:94` — that makes the arm's-length pin testable.

#### [ ] 2.5 Close out

Add the packages to `scripts/setup_r_env.sh` and `containers/Dockerfile.r-backends`,
update `dependency_risk`, add the test to `scripts/gates.yaml` Stage 5, adjust the
`partial_reason`, close B-006 in `docs/blocked.md`.

**Tranche 2 validation:** the new test executes rather than skips, `pytest -q` green,
status regenerated and diff-clean.

---

### [ ] Tranche 3 — SPiSCy-style preprocessing (B-009)

The largest design decision left, with a trap at 3.2.

#### [ ] 3.1 Decide the gating model before writing code

SPiSCy gates on 2-D geometry: `rectangleGate` / `polygonGate` over `mk1`/`mk2`
(`.upstream/spiscy/workflow/scripts/prelim_gating.R:81-112`, applied at `:215-216`).
Our existing `GatingRule` (`sampling/partitions.py:46-54`) is 1-D per-marker thresholds.
**These are not the same shape** and `GatingRule` cannot express a polygon. Add a 2-D gate
table, or state that only 1-D gating is supported.

#### [ ] 3.2 Check the transform convention before adopting FlowVS

`io/export.py:36` inverts as `sinh(x) * cofactor`, i.e. it assumes
`x_in = asinh(x_raw / cofactor) * cofactor`. SPiSCy uses
`flowVS::transFlowVS` (`.upstream/spiscy/workflow/scripts/transformation.R:159`). If the
conventions differ, the exporter's exact-inverse claim (`export.py:47-52`) becomes false.

Note also that `grep -i logicle` over the SPiSCy checkout returns **zero** matches. There
is no logicle in the reference implementation, so `config.py:186`'s
`Literal["asinh","logicle","none"]` is a schema choice SPiSCy does not motivate, and
`test_uninvertible_transformation_keeps_transformed_values_and_says_so` already documents
that logicle has no tested inverse.

#### [ ] 3.3 PeacoQC adapter

`bioconductor-peacoqc` is absent from `renv/` and from `scripts/setup_r_env.sh:107-113`;
SPiSCy pins 1.16.0 in `.upstream/spiscy/docker/rbase.yaml:20`. Entry point:
`PeacoQC(ff=, channels=, determine_good_cells=, save_fcs=, plot=, report=,
name_directory=, IT_limit=, MAD=)`.

**SPiSCy discards the removed-event index** — it only renames the report file — so our
`removed_events` writer must read `peacoqc_res$BadCell` or equivalent, which SPiSCy never
does. Use `rscript="Rscript-does-not-exist"` for negative tests, per the pattern that made
the Stage 5 unavailable path testable in every environment.

#### [ ] 3.4 Write the `removed_events` writer

The field exists at `io/dataset.py:54-56` with exactly the right columns and **has no
writer and no reader anywhere**; `load_dataset:261-269` does not pass it. Thread it into
the manifest and the report.

#### [ ] 3.5 Config key discipline

The three keys sit in `UNSUPPORTED_KEYS` (`config.py:39-57`) and must be moved out before
re-adding them as fields. Any new key must be read outside `config.py` / `examples.py` /
`schemas/`, **or** listed in `DECLARATIVE_ONLY` with an >80-character justification **and**
pinned in `SUPPORTED_VALUES` — enforced by `tests/unit/test_config_keys.py`.

Also correct `docs/blocked.md:301-302`, which claims the three keys are in
`DECLARATIVE_ONLY`. They are in `UNSUPPORTED_KEYS`; the two sets are asserted disjoint.

#### [ ] 3.6 Any R call must be added to `upstream_lock.yaml` in the same change

`docs/usage.md:252-257`: adding an R call means adding it to `adapter_calls`, or the
Stage 0 test fails by design. That test **skips** when no R is on `PATH`, so it will not
catch a mistake locally.

---

### [ ] Tranche 4 — workflow, CI, deployment

#### [ ] 4.1 Make `envs/` real, or say they are recipes

`python.yaml` and `openproblems.yaml` exist and **nothing consumes them**: no `conda:`
directive in the Snakefile, no `--use-conda` in `cli.py:50-65` or in the test harness, no
CI step. Either wire them into all nine rules, or state in the README that they are
human-run recipes.

#### [ ] 4.2 `rule status` is the only rule without `benchmark:`

`Snakefile:220-235`. Add it, then raise the assertion from `>=7` to `>=9` and add
`"status"` to the rule-name tuple (`test_stage11_workflow.py:166-172`). The current
threshold has exactly enough slack to hide the omission.

#### [ ] 4.3 Add a resume test

No test interrupts a run. The nearest is post-hoc deletion (`:69-89`). Add a
`SIGTERM`-mid-run test using `Popen` + `signal`, then assert byte-identical outputs on
rerun. `rerun-incomplete: true` is already set in both profiles, so the semantics are
configured but unverified.

#### [ ] 4.4 Add `snakemake --lint` to CI

Needs a valid config on the command line, because `Snakefile:28` loads config at parse
time.

#### [ ] 4.5 Container — record what cannot be verified

No apptainer, singularity, docker or podman exists on this node, so the image cannot be
built, run or digest-pinned here. Three concrete defects are fixable without building:
the base image is a tag not a digest (`:12`); it installs **none** of the 11 packages
CytoBatchFlagR needs; and `BATCHDX_UPSTREAM_DIR` (`:37`) is read by nothing. Update the
`NOT BUILT` banner, whose CRAN-unreachable premise is now stale.

#### [ ] 4.6 Run the R layer in CI

`scripts/setup_r_env.sh` is untracked, so there is nothing for CI to call. Replace
`r-base-core r-cran-jsonlite` (`batchdx-cyto-tests.yaml:31`) with the real toolchain, run
the script before `pytest`, and expect the staleness gate at `:57-71` to fail until the
status document is regenerated with new skip counts.

**Do not relax the skip guards** to make CI green. Gating a test on whether a package is
installed is what made the Stage 5 unavailable path untested in the first place.

#### [ ] 4.7 SLURM profile

`workflow/profiles/slurm/config.yaml:15-16` still has `REPLACE_ME` twice, and `:28` is a
standalone comment with no directive after it. The earlier plan attributed that comment
to `Snakefile:28`; it is in this file. `local/config.yaml:3` claims the integration suite
runs that profile; no test passes `--profile`.

---

### [ ] Tranche 5 — release readiness

#### [ ] 5.1 Commit the 13 untracked files

Everything that makes the R layer work exists only on disk, including the script eleven
tracked files reference.

#### [ ] 5.2 Multi-seed

Seed 473 is hard-coded in **8** places, not 5: `config.py:161`, `diagnostics/base.py:26`,
`examples.py:94`, `examples.py:154`, `simulate.py:60,209,237,254`. Thread one value from
`ProjectConfig.seed` and parametrise the scientific suite over at least five seeds, adding
the seed to the fingerprint.

#### [ ] 5.3 `tests/testthat/` for the R scripts

Five tracked R files; only `scripts/audit_upstream.R` is executed by CI.
`cytobatchflagr_detect.R` has never been executed at all. Highest-value target first:
`workflow/scripts/lib/binmat.R`, which is hand-mirrored against
`src/batchdx/cyto/io/binary.py` by its own admission — a byte-order or
row-major/column-major divergence there would be silent and would corrupt every adapter.

#### [ ] 5.4 One Python support story

Four different answers today: `requires-python = ">=3.10"`, CI pins 3.11,
`envs/python.yaml:29` pins 3.12, and there are no version classifiers. Reconcile, then
write `docs/compatibility_matrix.md` and a CI matrix.

#### [ ] 5.5 Dependency bounds

All twelve core dependencies have lower bounds only. `anndata>=0.10` already produced a
real user-visible bug: anndata 0.13 exposes a `None` key in `.layers`, which crashed
`sorted(adata.layers)` at `io/dataset.py:272`.

#### [ ] 5.6 `CHANGELOG.md` and `RELEASE.md`

Neither exists. Also add `docs/plans/` to the `docs/README.md` index, and decide whether
the roadmap itself is tracked.

#### [ ] 5.7 Ship synthetic-only, stated

Make the decision explicit in `README.md`, `docs/usage.md` and Stage 12's
`partial_reason`, and remove "a permissibly usable real cytometry example with provenance
and checksums" as a not-done item. **B-008 closes as a stated scope boundary, not as a
capability.** Every scientific result in this package rests on deterministic synthetic
fixtures; that must be on the face of the documentation, not only in `docs/blocked.md`.

#### [ ] 5.8 Populate Stage 12's `tests:` list

`scripts/gates.yaml:382` is `tests: []`, and `generate_status.py:73-74` maps an empty list
to DEFERRED unconditionally. **Stage 12 cannot leave DEFERRED by running tests**; its
list must be populated first.

---

## Final Validation

Run after the last tranche, not incrementally:

1. `pytest -q` — 0 failures, and no unexplained skips.
2. `pytest --json-report --json-report-file=build/report.json` then
   `python scripts/generate_status.py --report build/report.json`; the committed
   `docs/implementation_status.md` must diff clean, which CI already enforces.
3. `ruff check src tests scripts`.
4. Both audit scripts exit 0: `python scripts/audit_upstream.py` and
   `Rscript scripts/audit_upstream.R`.
5. Every declared gate test exists on disk and was executed.
6. A deliberate mutation per gate still reports failure — pointing a gate at a
   nonexistent test must produce `NOT RUN`, and `generate_status.py:88-90` must refuse to
   report PASSED while any declared test is missing, failed or skipped.
7. `git status` shows no untracked source file.
8. Every entry in `docs/blocked.md` either cites a command that currently fails or is
   marked RESOLVED.

## What cannot be verified on this machine

- **The container.** No apptainer, singularity, docker or podman on PATH. It cannot be
  built, run or digest-pinned here. Tranche 4 records the fixable defects and states the
  limit.
- **A real cytometry dataset.** None is obtainable under the current egress policy.
  Tranche 5.7 converts this into a stated scope boundary rather than pretending it is
  closable.
- **CI behaviour.** Every CI claim in this plan is read from YAML, never executed.
- **PeacoQC and the 11 CytoBatchFlagR packages.** Installability needs a solve that has
  not been run.
- **`snakemake --lint`.** Unexecuted, and its result depends on the Snakefile parsing with
  a valid config at parse time.

## Decision and Evidence Log

- 2026-09-30 — Plan reviewed. Found that the prior plan's Step 2.1 evidence was false in
  six particulars, and that 13 files including the R environment setup script are
  untracked. Both recorded as findings above rather than deleted.
- 2026-09-30 — Decided to ship synthetic-only and state it, rather than leave B-008 as an
  open blocker that never converges.
- 2026-09-30 — Decided to execute tranche by tranche, verifying after each, rather than
  in one pass. The previous tranche's parallel agents corrupted each other's test
  results in a shared checkout; the test suite is the only acceptable arbiter.
- 2026-09-30 — Implementation not started.

## Risks and Assumptions

- **Tranche 2 may not converge.** If the 11 packages cannot be resolved for this platform,
  or upstream's unconditional `library()` cannot be satisfied, B-006 stays open. An
  honest `not_available` is the correct outcome; a fabricated pass is not.
- **Tranche 3's transform convention is a genuine trap.** If FlowVS and the exporter's
  inverse disagree, adopting FlowVS silently invalidates the exact-inverse claim at
  `export.py:47-52`. Read both before choosing.
- **`anndata` compatibility is a live class of bug, not a one-off.** `pyproject.toml`
  pins `>=0.10` with no upper bound; 0.13 changed `.layers` semantics. The 0.13 defect is
  fixed, but only a version matrix catches the next one.
- **A `~/.Renviron` on the build host can silently redirect R installs.** It overrode
  `R_LIBS_USER` and wrote into a personal library. Mitigated by `--vanilla` in
  `setup_r_env.sh` and in `run_r_script`; a CI runner with the same file and a script that
  omits `--vanilla` would reintroduce it.
- **Single seed still gates the scientific claims** until 5.2 lands. The Stage 8
  resampling is explicitly seeded and a test asserts reproducibility from that seed, but a
  one-seed run is not a sweep.
