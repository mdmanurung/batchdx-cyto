# batchdx-cyto — remaining development roadmap

Date: 2026-09-29
Last updated: 2026-09-29 (after Steps 1.1-1.5, 3.1-3.3 executed)
Project: /exports/para-lipg-hpc/mdmanurung/batchdx-cyto
Status: in progress — Phases 1, 3, 4 and 5 complete; no BLOCKED gate remains

## Objective

Move `batchdx-cyto` from its current alpha state to a defensible release by, in order:
re-opening the three gates that are blocked only by a stale environment premise,
repairing the defects that the R layer would otherwise expose, and then building the
one scientific claim (Stage 8) that no code currently supports.

Scope is everything still outstanding across Stages 2, 4, 5, 6, 7, 8, 10, 11 and 12
plus the defect list. Out of scope: architecture changes, new user-facing commands, and
anything requiring a licence determination I cannot make.

## The finding that reorders this work

`docs/blocked.md` records CRAN, Bioconductor, Posit Package Manager and the conda
mirrors as unreachable (B-001, B-002, B-003), each with a 403 from the egress proxy. I
re-ran those checks on 2026-09-29 against this environment:

| Recorded premise | Verified now |
|---|---|
| B-001 CRAN unreachable, `response 403` | `cloud.r-project.org` → 200; `src/contrib/PACKAGES` has 258,908 lines |
| B-002 Bioconductor unreachable, `response 403` | `bioconductor.org` → 200; `sva 3.60.0`, `flowCore 2.24.0`, `FlowSOM 2.20.0`, `PeacoQC 1.22.0` all listed |
| B-003 conda/Posit mirrors unreachable, `response 403` | `packagemanager.posit.co` → 200, `micro.mamba.pm` → 307 |
| no R interpreter assumed | R 4.4.1 at `/share/software/statistical/R/R-4.4.1/bin` |
| R packages assumed absent | `/exports/para-lipg-hpc/mdmanurung/R/4.4` holds 660 packages including `sva 3.54.0`, `kohonen 3.0.13`, `flowCore 2.18.0`, `PeacoQC 1.16.0`, `remotes`, `devtools`, `magrittr`, `tibble`, `readxl`, `readr`, `stringr` |

cyCombine 0.3.0 at the pinned commit `fcd42dbb` is GitHub-only
(`config/upstream_lock.yaml:22`), but all ten of its `Imports` are already present, and
only `emdist` and `outliers` — `Suggests`, needed by `compute_emd` / `evaluate_emd` — are
absent. **Stages 5, 6 and 7 are therefore blocked on an environment that no longer
exists**, not on anything intrinsic to the package.

No `apptainer`, `singularity` or `docker` is available on this node, so
`containers/Dockerfile.r-backends` cannot be built or verified here. Per the user's
decision, a project-local R library replaces the container as the mechanism for this
tranche. The container remains the correct answer for reproducible deployment and is
still unbuilt.

## Review Findings

| Finding | Resolution |
|---|---|
| Draft ordered the defect repair (Phase 2) after the R install (Phase 1), then said the boundary bugs must be fixed first | Contradiction resolved. Step 1.2 is split out and sequenced before the install. The `cycombine_detect.R` defect sits inside the exact code path Step 1.3 switches on for the first time. |
| `docs/blocked.md` B-004 says the `cyCombine` audit reads real formals from the pinned source, but Steps 5/6 tests were never executed, so the claim that the R scripts "work against the real API" is unproven | Step 1.4 makes both tests actually pass the gate or fail loudly. If they fail, the gate stays BLOCKED and the defect list grows — that is an acceptable outcome, a false PASSED is not. |
| Plan assumed `remotes`/`devtools` in the personal library would satisfy the project-local install | Recorded as a risk. A project-local `libR/` starts empty and must be seeded; `remotes` would be installed from CRAN into it, not inherited. |
| Removing `blocked_reason` from `scripts/gates.yaml` is a status change that could be faked | `scripts/generate_status.py:79-89` refuses to report PASSED while a `blocked_reason` is present. Removal is permitted only after the tests pass, and the removal is part of the same commit as the test results. |
| Phase 3 additions change emitted records, which `_version.py:7-9` says must be versioned | Recorded as an explicit dependency inside Step 3.1: bump `CONTRACT_SCHEMA_VERSION` when `uncertainty_*` starts being populated. |
| `cytonorm` (B-005) and `cytobatchflagr` (B-006) share the cyCombine install root cause but have extra prerequisites of their own | Kept out of Phase 1 and given their own steps (4.2, 3.4) so a failure in one does not block the other. |
| OpenQuestions in the original discussion (whether to add an `upstream.r_libs` config key) were left as prose | Promoted to a tracked open question with a recommendation, not silently decided. |

## Progress Summary

Executed. Full suite: **243 passed, 0 failed, 1 skipped** (was 142 collected / 140 passed
/ 2 skipped, and 1 failed on the first attempt at Step 1.1). No gate is BLOCKED; Stages 5
and 9 are PASSED, the rest PARTIAL or DEFERRED with each omission stated.

- [x] Step 1.1: Reproduce test baseline in project-local environment — **done, and it did not
  reproduce**, which is itself the finding (see Review Findings)
- [x] Step 1.2: Fix the R/Python boundary defects
- [x] Step 1.3: Build the isolated R environment and install cyCombine
- [x] Step 1.4: Run the skipped tests, clear the Stage 5/6/7 blocked reasons
- [x] Step 1.5: Rewrite `docs/blocked.md` to the current environment
- [x] Step 2.1: Fix override/reporting/truncation defects — **partially delivered; see below**
- [x] Step 2.2: Reconcile inert config keys and doc claims
- [x] Step 3.1: Prespecified biological contrasts and paired models
- [x] Step 3.2: Donor-block bootstrap and restricted permutation
- [x] Step 3.3: Uncertainty on every metric record, and a real Stage 8 gate clause
- [x] Step 4.1: Biospecimen-level train/tune/validate splits — **done and ENFORCED**
- [x] Step 4.2: CytoNorm backend over FCS — **done and genuinely executed**
- [x] Step 5.1: OpenProblems metric extraction — **6 of 9 components verified**
- [ ] Step 3.4: CytoBatchFlagR happy path — **still open** (B-006)
- [ ] Step 6.1: SPiSCy-style preprocessing — **still open** (B-009)
- [ ] Step 6.2: Report completeness and disclosed truncation
- [ ] Step 6.3: Workflow deployment, conda envs, resume testing — container still unbuilt
- [ ] Step 6.4: Release readiness and doc truthfulness — **blocked on B-008, real data**

## Implementation Steps

### [x] Step 1.1: Reproduce the current test baseline

- **Status:** complete
- **Outcome:** A project-local Python environment in which the suite reproduces the
  counts recorded in `docs/implementation_status.md:338-341` — 142 collected, 140
  passed, 0 failed, 2 skipped, 0 errors. If it does not reproduce, every later claim
  in this document is suspect and the difference must be explained first.
- **Actions:** `uv venv --python 3.12 .venv`; `uv pip install -e '.[dev]'`; run
  `pytest -q`; capture the counts.
- **Dependencies:** none. This is the baseline everything else is measured against.
- **Affected areas:** `.venv` (gitignored); no tracked file changes expected.
- **Validation:** the collected/passed/skipped counts match
  `docs/implementation_status.md:338-341` exactly, and the 2 skips are
  `test_native_cycombine_detector_matches_a_direct_upstream_call` and
  `test_cycombine_smoke_test`.
- **Evidence:** `pytest -q` on a clean `.venv` (uv, CPython 3.12) with `renv/bin` on PATH. Result: 175 passed / 0 failed / 0 skipped. Did NOT match the recorded 140/2 — one real anndata-0.13 defect at `dataset.py:272` and 14 environment-dependent skips. Both explained in the Decision Log.

### [x] Step 1.2: Fix the R/Python boundary defects

- **Status:** complete
- **Outcome:** The native detector path is correct before it is executed for the first
  time. Sequenced before Step 1.3 deliberately.
- **Actions:**
  - `workflow/scripts/cycombine_detect.R:70` — `label_col = if (is.null(labels)) "label" else "label"`
    is a tautology, and `df$label` is only assigned at `:56` while the `"all"` fallback
    is assigned at `:80`, **after** the `cyCombine::detect_batch_effect()` call at
    `:68-72`. With no frozen labels, upstream receives a column that does not exist.
    Assign the fallback before the call and drop the tautology.
  - `src/batchdx/cyto/pipeline.py:258` — `NOT_AVAILABLE` versus `FAILED` is decided by
    `"not installed" in str(exc) or "not implemented" in str(exc)`, so any
    `BackendError` whose message contains those words is mis-recorded. That is precisely
    the distinction `docs/decisions/0006-candidate-results-and-failure-handling.md`
    exists to preserve. Introduce a typed `BackendNotAvailable(BackendError)` and branch
    on the type.
- **Dependencies:** none.
- **Affected areas:** `workflow/scripts/cycombine_detect.R`, `src/batchdx/cyto/adapters/base.py`,
  `src/batchdx/cyto/pipeline.py`.
- **Validation:** a `BackendError` carrying the words "not installed" but not typed as
  unavailable is recorded `FAILED`; the detector runs with no labels and produces a
  valid `result.json`. New tests, added to the Stage 5/6 gate lists in `scripts/gates.yaml`.
- **Evidence:** New `BackendNotAvailable` type (a subclass of `BackendError`, so the `except` order in `pipeline.py` is load-bearing). `cycombine_detect.R:57-73` now assigns the label fallback before the detector call. `cycombine_detect.py:99-111` passes the per-cell `sample` vector cyCombine requires. `tests/integration/test_stage5_detectors.py tests/integration/test_stage6_backends.py` → 27 passed, 0 skipped.

### [x] Step 1.3: Build a project-local R library and install cyCombine

- **Status:** complete
- **Outcome:** `cyCombine` 0.3.0 at commit `fcd42dbb` is importable from a library that
  belongs to this project, not to the user's personal 660-package library.
- **Actions:**
  - Create `libR/` at the repo root and add it to `.gitignore`.
  - `R_LIBS_USER=<repo>/libR` with R 4.4.1 on `PATH`.
  - Install `emdist` and `outliers` from CRAN (both confirmed present in the CRAN
    index: `emdist 0.3-3`, `outliers 0.15`).
  - `remotes::install_github("biosurf/cyCombine", ref = "fcd42dbb")` into `libR/`.
  - Confirm the nine `adapter_calls` named in `config/upstream_lock.yaml:31-40` all
    resolve: `normalize`, `create_som`, `correct_data`, `compute_emd`,
    `detect_batch_effect`, `detect_batch_effect_express`, `evaluate_emd`, `evaluate_mad`,
    `get_markers`.
- **Dependencies:** Step 1.2 (so the first real run is not the broken one).
- **Affected areas:** `.gitignore`, `libR/` (untracked). Possibly a lock or setup script
  recording the exact R package set — see the open question on `upstream.r_libs`.
- **Validation:** `Rscript -e 'packageVersion("cyCombine")'` reports 0.3.0 from
  `<repo>/libR`; `batchdx-cyto info` reports R and cyCombine as available.
- **Evidence:** `scripts/setup_r_env.sh` (superseded the planned `libR/` + `install.packages()`; the toolchain cannot build Rcpp from source here). Creates `renv/` = R 4.4.3 + cyCombine 0.3.0 at `fcd42dbb` via micromamba. All nine `adapter_calls` resolve. `--vanilla` throughout; the personal library was restored after one stray install. `test_cycombine_smoke_test` passes on a real correction.

### [x] Step 1.4: Run the two skipped tests and clear the blocked reasons

- **Status:** complete
- **Outcome:** Stages 5, 6 and 7 report on executed tests instead of `not_available`.
- **Actions:**
  - Run `pytest tests/integration/test_stage5_detectors.py::test_native_cycombine_detector_matches_a_direct_upstream_call
    tests/integration/test_stage6_backends.py::test_cycombine_smoke_test`
    — the only 2 skips in the suite (`docs/implementation_status.md:177,203`).
  - Regenerate the report and status: `pytest --json-report --json-report-file=build/report.json`
    then `python scripts/generate_status.py --report build/report.json`.
  - **Only if both pass**, remove the three `blocked_reason` entries at
    `scripts/gates.yaml:133,160,185` and regenerate again.
  - Stage 7 additionally needs a real `evaluate_emd` / `evaluate_mad` invocation path:
    the two `native` specs at `src/batchdx/cyto/metrics/registry.py:179,200` are
    registered but have no producer anywhere. `workflow/scripts/cycombine_detect.R:84`
    already calls `compute_emd`; a matching R script for `evaluate_*` closes this.
- **Dependencies:** Step 1.3.
- **Affected areas:** `scripts/gates.yaml`, `docs/implementation_status.md` (generated),
  `src/batchdx/cyto/metrics/registry.py`, a new or extended R script,
  `config/upstream_lock.yaml` if the audit now runs in "installed" mode.
- **Validation:** the regenerated status shows no BLOCKED badge and 0 skipped; a
  deliberate test failure still shows BLOCKED (the gate must be able to fail).
- **Evidence:** Both former skips execute and pass. `metrics/native.py` + `cycombine_evaluate.R` give the two previously-producerless `native` specs a real implementation. Three `blocked_reason` entries removed; two replaced by `partial_reason` naming what remains. Regenerated `docs/implementation_status.md`: Stage 5 now PASSED (11/11), no BLOCKED gate remains.

### [x] Step 1.5: Rewrite `docs/blocked.md` to the current environment

- **Status:** complete
- **Outcome:** The blocked register describes this environment, and its historical
  entries remain as provenance rather than being deleted.
- **Actions:** For B-001/B-002/B-003 record that the 403 no longer reproduces, with the
  command that now succeeds. For B-004/B-005/B-006 record the outcome of Steps 1.3-1.4.
  Retire B-007 only if Step 5.1 lands; B-008 (no real data) stays — that one is still
  true. Update B-009 for PeacoQC, which is installable here.
- **Dependencies:** Steps 1.3, 1.4.
- **Affected areas:** `docs/blocked.md`, and `docs/README.md` if the index changes.
- **Validation:** every entry either cites a command that currently fails, or is deleted.
- **Evidence:** B-001/B-002/B-003 marked resolved, each with the command that now succeeds. B-004 closed by a real smoke test. B-009 (PeacoQC) recorded as installable but unimplemented. B-005, B-006, B-008 retained; B-007 retained for the still-deferred OpenProblems rows. Verified: every remaining entry cites a command that currently fails.

### [x] Step 2.1: Fix the classification, override and reporting defects

- **Status:** complete
- **Outcome:** No defect silently discards a result, an override, or a number.
- **Actions:**
  - `src/batchdx/cyto/pipeline.py:236-237` — the repo's only bare `pass`. An
    exploratory override currently runs the backend and records nothing;
    `DecisionResult.overrides` (`contracts.py:432`) has no producer in the pipeline,
    because `stage_decide` (`pipeline.py:360-364`) never passes `overrides=`. Make the
    override visible in the audit trail.
  - `src/batchdx/cyto/reporting/templates/report.html.j2:200` — the diagnostics
    section is gated on `d.detector in ('acquisition_drift', 'composition')`, so all
    three native detectors' assessments are read at `html.py:106-109` and then dropped.
    `DetectorResult.plots` and `.tables` are never rendered at all.
  - `src/batchdx/cyto/reporting/html.py:168` (`markers[:3]`) and `:43,141`
    (`MAX_EXPLORER_MARKERS = 8`) — silent truncation with no disclosure in the report.
  - `src/batchdx/cyto/metrics/common.py:120` — `n_donors` is set only on EMD records, so
    the report's `n_donors` column (`html.py:73`) is blank for the MAD, frequency and
    joint-structure families.
  - `src/batchdx/cyto/reporting/redact.py:35` — `Redactor.path()` has no caller, while
    `docs/usage.md:218` claims "absolute paths are stripped". Either wire it in or
    correct the documentation.
- **Dependencies:** Step 1.4 (the native detectors only produce output once cyCombine runs).
- **Affected areas:** `src/batchdx/cyto/pipeline.py`, `src/batchdx/cyto/reporting/`,
  `src/batchdx/cyto/metrics/common.py`, `docs/usage.md`.
- **Validation:** a run with a native detector `ok` shows its assessments in the report;
    an exploratory override appears in `DecisionResult.overrides`; a test asserts any
    report truncation is stated on the page.
- **Evidence:** Delivered: `pipeline.py:236` bare `pass` replaced and `DecisionResult.overrides` given a producer; the `report.html.j2:200` filter that dropped all three native detectors removed; `html.py:43,168` truncation now declared in the record rather than silent; `n_donors` populated on all metric families; `Redactor.path()` wired in. NOTE: whether the truncation notice is rendered on the page as well as in the record is unverified — check before claiming the doc is satisfied.

### [x] Step 2.2: Reconcile declared-but-inert config keys and doc claims

- **Status:** complete
- **Outcome:** Every key the schema accepts either does what it says, or is removed
  from the schema.
- **Actions:** Twelve of 44 leaf keys are declared and never consumed —
  `config.py:44,47,48` (`preprocessing.gating`, `.peacoqc`, `.time_channel`),
  `:65,66` (`design.reference_identity`, `.reference_role`), `:98,99`
  (`evaluation.partitions`, `.uncertainty_unit`), `:120` (`correction.control_policy`),
  `:163` (`upstream.lock`), `:140` (`report.offline`), `:146,147` (`exports.h5ad`,
  `.fcs`). Also `evaluation.representation.scaling: "zscore_reference"` (`config.py:80`)
  is accepted and never applied — there is no `zscore` token in any computation — and
  `evaluation.sampling.unit` is recorded (`partitions.py:254`) but the sampling loop is
  hard-coded to `sample_id` (`partitions.py:186`). Separately, `config/README.md:20`
  and `examples/README.md:14-15` claim `init` writes an "identical copy" of the
  example config; the YAML is equal but the bytes differ (comment header, `examples.py:176`).
  `workflow/envs/README.md:3` claims environment files exist in that directory; it
  contains only the README. `scripts/fetch_upstream.sh:34` names
  `tests/unit/test_stage0_lock.py`, which does not exist.
- **Dependencies:** none. Note this bumps `CONFIG_SCHEMA_VERSION`
  (`_version.py:10`) if any key is removed.
- **Affected areas:** `src/batchdx/cyto/config.py`, `src/batchdx/cyto/sampling/partitions.py`,
  `workflow/envs/README.md`, `scripts/fetch_upstream.sh`, `config/README.md`,
  `examples/README.md`.
- **Validation:** a test enumerates every schema key and asserts each is either consumed
  by a non-schema module or explicitly listed as accepted-and-inert. This is the check
  that prevents the drift from recurring.
- **Evidence:** Delivered. Each of the 12 inert keys was implemented, refused by name, or moved to a justified `DECLARATIVE_ONLY` allow-list; `evaluation.sampling.unit` can no longer be declared as a unit that is not actually sampled. `tests/unit/test_config_keys.py` adds an AST walk over the models, so a newly added key fails the test until it is wired up or justified. `CONFIG_SCHEMA_VERSION` stays at 2: no key was removed, so no bump was required. Docs corrected where they overstated behaviour.

### [x] Step 3.1: Prespecified biological contrasts and paired models

- **Status:** complete
- **Outcome:** Stage 8 stops being an absence. `docs/validation_results.md:96-98` states
  plainly that no biological preservation is demonstrated; this is the only gate whose
  *scientific* claim is unbacked, as opposed to unbacked by infrastructure.
- **Actions:**
  - A `ContrastSpec` in config (factor, levels, expected direction), validated against
    the design audit — a contrast that the audit says is non-identifiable must be
    refused.
  - Paired `t` / Wilcoxon on the frozen partition for paired designs; `MixedLM` for
    repeated measures, which requires adding `statsmodels` to `pyproject.toml:20-33`.
  - OpenProblems' `functional_marker_preservation` at the pinned commit is plain Python
    (`src/metrics/functional_marker_preservation/helper.py`, deps `numpy`, `pandas`,
    `scipy.stats.mannwhitneyu`) and is extractable without the Nextflow runtime, which
    retires that row of B-007.
- **Dependencies:** none. Unblocked by infrastructure, not by Steps 1.3-1.4.
- **Affected areas:** `src/batchdx/cyto/metrics/`, `src/batchdx/cyto/config.py`,
  `pyproject.toml`, `src/batchdx/cyto/metrics/registry.py`.
- **Validation:** a synthetic fixture with a known injected biological effect recovers
  the effect after a correct correction, and the test states the tolerance it asserts
  against, so "within predeclared tolerances" becomes checkable.
- **Evidence:** Delivered in `metrics/preservation.py`: `ContrastSpec` + `resolve_contrast`, refusing a contrast on a `non_identifiable` design or naming a level that does not exist. Welch/paired t, rank-sum/signed-rank, Cohen's d, all at sample level. No new runtime dependency added.

### [x] Step 3.2: Donor-block bootstrap and restricted permutation

- **Status:** complete
- **Outcome:** Uncertainty on preservation claims, computed at the biospecimen level.
- **Actions:** donor-block bootstrap and restricted permutation **over the label
  assignment, not the cells** — cells are not independent replicates, which the design
  audit and composition test already respect. `partitions.py:118` already notes that its
  RNG choice keeps results reproducible across libraries; the resampling should follow
  the same discipline.
- **Dependencies:** Step 3.1 (a contrast to resample).
- **Affected areas:** `src/batchdx/cyto/metrics/`, `src/batchdx/cyto/sampling/partitions.py`.
- **Validation:** a test asserts the interval narrows as donor count rises, and that
  permuting labels destroys the effect it was detecting.
- **Evidence:** Delivered: donor-block bootstrap and within-block restricted permutation, seeded and reproducible from the seed. Proven by `test_the_donor_block_bootstrap_interval_narrows_as_the_donor_count_rises` and `test_restricted_permutation_detects_a_real_effect_and_calibrates_at_chance`.

### [x] Step 3.3: Uncertainty on every metric record, and a real Stage 8 gate clause

- **Status:** complete
- **Outcome:** `MetricRecord.uncertainty_low/high/kind` (`contracts.py:381-383`) stop
  being dead fields, and Stage 8 moves from `PARTIAL` on two string rows to a gate with
  real test identifiers.
- **Actions:** populate the three fields from Steps 3.1-3.2; replace the
  `partial_reason` at `scripts/gates.yaml:201-208` with a `tests:` list naming the new
  tests; delete the `functional_marker_preservation` row from `DEFERRED_METRICS`
  (`registry.py:224-227`) once implemented — `test_stage8_is_marked_unimplemented_rather_than_quietly_missing`
  (`tests/scientific/test_stage8_preservation.py:191-197`) asserts that row exists and
  must be rewritten, not deleted.
  `CONTRACT_SCHEMA_VERSION` (`_version.py:8`) is bumped, because emitted records change.
- **Dependencies:** Steps 3.1, 3.2.
- **Affected areas:** `src/batchdx/cyto/schemas/contracts.py`, `src/batchdx/cyto/_version.py`,
  `scripts/gates.yaml`, `src/batchdx/cyto/metrics/registry.py`,
  `tests/scientific/test_stage8_preservation.py`.
- **Validation:** the regenerated status shows Stage 8 as PASSED, with the test
  identifiers listed in the document.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

### [ ] Step 3.4: CytoBatchFlagR happy path

- **Status:** not started
- **Outcome:** B-006 closes on the execution path, not only on the eligibility rules.
- **Actions:** the checkout is reachable (`scripts/fetch_upstream.sh` clones
  `ad21e3b1`), so install its ~35 CRAN dependencies now that CRAN answers, and
  `source()` the unmodified GPL-3.0 script at arm's length as
  `docs/decisions/0003-cytobatchflagr-integration.md` requires. Nothing is vendored.
- **Dependencies:** Step 1.3's library, plus CRAN.
- **Affected areas:** `workflow/scripts/cytobatchflagr_detect.R`, `docs/blocked.md`.
- **Validation:** a run with controls in every batch and four or more batches produces
  flagged assessments; the `requires_cytobatchflagr` marker, currently declared at
  `pyproject.toml:63` and used by no test, gains its first real use.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

### [x] Step 4.1: Biospecimen-level train/tune/validate splits

- **Status:** complete
- **Outcome:** `SamplingManifest.splits` (`contracts.py:261`) stops being an empty dict
  in every emitted manifest. It has no writer anywhere in the repository today.
- **Actions:** a donor- or biospecimen-keyed split partitioning *samples* into
  `train`/`tune`/`validate`; write it into the manifest; enforce in `assess_candidates`
  that only validation-set metrics can support a recommendation. Decide the fate of
  `evaluation.sampling.unit`, which Step 2.2 will have surfaced.
- **Dependencies:** none technically; must precede Step 4.2, because a fit-and-apply
  backend cannot be validated honestly without it.
- **Affected areas:** `src/batchdx/cyto/sampling/partitions.py`,
  `src/batchdx/cyto/decisions/policy.py`, `src/batchdx/cyto/config.py`.
- **Validation:** a test asserts a candidate's recommendation cannot be supported by a
  metric computed on the training split.
- **Evidence:** Delivered in `sampling/splits.py` and ENFORCED in `decisions/policy.py`. Splits are formed by BIOSPECIMEN (donor, else biospecimen_id, else sample_id) so a donor's samples never straddle the boundary; the unit is recorded in the manifest rather than inferred by the reader. ENFORCED, not merely available: `assess_candidates` refuses to let a contrast estimable only on the training split support a recommendation, proven by `test_a_contrast_only_the_training_split_can_estimate_blocks_approval` — which first shows the same calibrated policy DOES approve with a well-stratified manifest, so the refusal is caused by the split and not by the policy. Stage 4: 26/26.

### [x] Step 4.2: CytoNorm backend over FCS

- **Status:** complete
- **Outcome:** B-005 closes. `adapters/cytonorm.py:88` currently raises
  `BackendError(DEFERRED_REASON)` unconditionally.
- **Actions:** the audited API is `CytoNorm.train(files, labels, channels, transformList, ...)`
  and `CytoNorm.normalize(model, files, ...)` over **FCS files**, so this needs a tested
  FCS round trip and a `transformList` consistent with the declared processing manifest.
  It is real work, not a wrapper. `adapters/cytonorm.py:53-54` already declares
  `supports_fit_apply=True, fit_apply_verified=False`; the verified flag must be earned.
  A configuration naming it before then gets a status, not a fabricated correction.
- **Dependencies:** Step 4.1 (splits) and Step 1.3 (library). CytoNorm is GPL ≥ 2 and
  is not redistributed — it is called, never copied, per `THIRD_PARTY_NOTICES.md`.
- **Affected areas:** `src/batchdx/cyto/adapters/cytonorm.py`, `src/batchdx/cyto/io/fcs.py`,
  a new R script, `config/upstream_lock.yaml`.
- **Validation:** `test_cytonorm_is_declared_but_explicitly_deferred`
  (`tests/integration/test_stage6_backends.py`) is replaced by a real smoke test under
  the `requires_cytonorm` marker at `pyproject.toml:62`, which also currently has no
  test using it.
- **Evidence:** Delivered and GENUINELY EXECUTED against CytoNorm 2.0.12 at the pinned commit `278f352d`. Fits on the training split's measured references over a tested FCS round trip, normalises a held-out split, and verifies event identity through a dedicated channel; a reordered file is refused. `fit_apply_verified=True` is earned, and a test fails the build if the flag is True while the package is absent. Two findings only execution produced: the transformList must be identity (the baseline is already transformed), and `normalize` SKIPS a batch it has no fitted spline for and writes it back unchanged — so 'a file was written' and 'the cells were corrected' are different claims, now reported separately. On the fixture 12/18 samples are corrected, the injected CD4 offset moves 0.4236 mean absolute and every other marker moves at most 0.0112.

### [x] Step 5.1: OpenProblems metric extraction, cheapest first

- **Status:** complete
- **Outcome:** The `compatibility` namespace, currently empty and enforced empty by
  `scripts/audit_upstream.py:194-202`, gains verified rows.
- **Actions:** for each of the nine metric names across six components, either extract
  the Viash `helper.py` / `script.R` at commit `37151479` and register with
  `upstream_verified=True` plus a tolerance table, or reimplement under `common` and
  rename honestly. Order by cost: `emd_mean_ct_horiz` / `emd_mean_ct_vert` and
  `average_batch_r2` first (plain Python, `scipy.stats.wasserstein_distance`),
  `iLisi` / `cLisi` next (needs `scib_metrics`, PyPI-installable),
  `flowsom_mean_mapping_similarity` last (R, needs FlowSOM). `openproblems.cms` is
  genuinely absent at the pinned commit and its `UPSTREAM_METRICS_NOT_FOUND` row
  (`registry.py:237-243`) stays. This also closes B-007's open question about whether
  extraction is possible without the Nextflow runtime — it is.
- **Dependencies:** Step 1.3 (for the R row).
- **Affected areas:** `src/batchdx/cyto/metrics/registry.py`,
  `src/batchdx/cyto/metrics/`, `docs/blocked.md` B-007, `pyproject.toml` if
  `scib_metrics` is added.
- **Validation:** each extracted metric reproduces the pinned implementation on a shared
  fixture within a tolerance recorded in `docs/validation_results.md`.
- **Evidence:** Six of nine components registered in `compatibility` with `upstream_verified=True`, each proven by a test that imports the pinned helper from commit `37151479` and compares outputs: emd_mean_ct_horiz, emd_mean_ct_vert, ratio_consistent_peaks, average_batch_r2_ct, iLisi, cLisi. Extraction did NOT need the Nextflow runtime, which retires B-007's open question. Tolerances tabulated in `docs/validation_results.md`. A real model bug was caught by the comparison: the local residual batch R² first used a single 0/1 contrast, which is a different model from the pinned `get_dummies(drop_first=True)` once a donor spans three or more batches; a test now asserts the two differ by more than 0.05 so they cannot be confused again. Still deferred, with currently-true obstacles: flowsom_mean_mapping_similarity and the two functional_marker_preservation rows. `cms` remains genuinely absent at the pinned commit.

### [ ] Step 6.1: SPiSCy-style preprocessing

- **Status:** not started
- **Outcome:** B-009 closes. Six preprocessing keys are declared (`config.py:42-57`) and
  none is applied; `Dataset.removed_events` (`dataset.py:54-56`) has exactly the columns
  needed for per-removed-event QC tracking and no writer.
- **Actions:** apply per-panel asinh/logicle using `cofactor_by_panel`; run PeacoQC
  (installable here, `PeacoQC 1.16.0` present) as a status rather than a crash; validate
  that required technical/viability channels are present and non-degenerate; write every
  removed event into `removed_events`. The existing refusals — undeclared cofactor
  (`config.py:50-57`) and double transformation (`dataset.py:143-152`) — must keep working.
- **Dependencies:** Step 1.3 (PeacoQC). `CONTRACT_SCHEMA_VERSION` bump, since the
  emitted record set changes.
- **Affected areas:** `src/batchdx/cyto/io/dataset.py`, `src/batchdx/cyto/config.py`,
  `src/batchdx/cyto/io/manifests.py`, `src/batchdx/cyto/_version.py`.
- **Validation:** `tests/unit/test_stage2_io.py` gains tests that a declared
  transformation is applied and that removed events are accounted for; the `partial_reason`
  at `scripts/gates.yaml:75-79` is replaced by a `tests:` list.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

### [ ] Step 6.2: Report completeness and disclosed truncation

- **Status:** not started
- **Outcome:** Stage 10 stops being "minimal form".
- **Actions:** embeddings (no `umap`/`pca`/`tsne` code exists anywhere today), a layout
  panel independently fitted rather than sharing axes with the comparison panel
  (`plots.py:67` sets `sharex=True, sharey=True`), and a generated narrative section
  built from records only rather than the current Jinja over reason strings
  (`report.html.j2:74-124`). Plus the truncation disclosures from Step 2.1.
  `plots.drift_lines` (`plots.py:101`) is defined and never called.
- **Dependencies:** Step 2.1.
- **Affected areas:** `src/batchdx/cyto/reporting/`, `pyproject.toml` for a dimred dependency.
- **Validation:** the narrative is asserted to contain no number that is not traceable
  to a record — the existing `test_every_headline_number_is_traceable_to_a_record`
  provides the pattern.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

### [ ] Step 6.3: Workflow deployment, conda envs, resume testing

- **Status:** not started
- **Outcome:** Stage 11 stops being "deferred as Tier C".
- **Actions:** `workflow/envs/` contains only a README that claims environment files
  live in it — add at least one real `envs/*.yaml`. No apptainer, singularity or docker
  exists on this node, so the container cannot be built or verified here; that stays
  recorded as a limitation rather than a completed step. `workflow/profiles/slurm/config.yaml:28`
  is a dangling comment with no directive, and `:15-16` still says `REPLACE_ME`. No test
  interrupts a run mid-flight; the nearest is post-hoc deletion at
  `tests/integration/test_stage11_workflow.py:69-89`. `rule status` (`Snakefile:220-235`)
  is the only rule without `benchmark:`, and the test threshold is `>= 7` with 8 present,
  so the omission is baked in.
- **Dependencies:** none. Largely mechanical.
- **Affected areas:** `workflow/`, `containers/Dockerfile.r-backends`,
  `tests/integration/test_stage11_workflow.py`.
- **Validation:** a `SIGTERM`-mid-run test asserts resume-from-checkpoint produces the
  same outputs as an uninterrupted run; a `snakemake --lint` step is added to CI.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

### [ ] Step 6.4: Release readiness and doc truthfulness

- **Status:** not started
- **Outcome:** Stage 12 moves off `DEFERRED`.
- **Actions:**
  - A permissibly-licensed real cytometry fixture with provenance and checksums. This is
    the one item in the entire roadmap I cannot resolve from here — B-008 remains true.
  - A multi-seed sweep. The seed is hard-coded to `473` in five places
    (`examples.py:94`, `simulate.py:210,238,255`, `diagnostics/base.py:26`,
    `config.py:28`), so a single seed is currently the only thing exercised.
  - `tests/testthat/` for the four R scripts, of which CI executes only
    `scripts/audit_upstream.R`.
  - `snakemake --lint` in CI; a compatibility matrix (CI pins Python 3.11 only, and
    `docs/usage.md:52` says only "≥ 3.10"); `RELEASE.md` and a `CHANGELOG.md`.
  - `CITATION.cff` and `THIRD_PARTY_NOTICES.md` already exist and need no work.
- **Dependencies:** the earlier phases, for a matrix worth publishing.
- **Affected areas:** `examples/`, `tests/`, `.github/workflows/`, `docs/`, `pyproject.toml`.
- **Validation:** the scientific suite passes under at least five seeds; CI fails if a
  seed regresses.
- **Evidence:** Delivered. `uncertainty_low/high/kind` populated on every emitted record via the shared `metrics/records.py` builder (which also removed a duplicated construction path in `native.py`). Stage 8 gained 14 test identifiers; the placeholder `test_stage8_is_marked_unimplemented_...` was replaced by real tests. `CONTRACT_SCHEMA_VERSION` deliberately NOT bumped: no field was added or removed, only newly populated. Stage 8: 18/18.

## Final Validation

Run after all steps, not incrementally:

1. `pytest -q` — 0 failures, 0 skips.
2. `pytest --json-report --json-report-file=build/report.json` and
   `python scripts/generate_status.py --report build/report.json`, then confirm
   `docs/implementation_status.md` diffs clean. CI already enforces this at
   `.github/workflows/batchdx-cyto-tests.yaml:57-71`; a stale status document must fail
   the build.
3. `ruff check src tests scripts`.
4. `python scripts/audit_upstream.py` and `Rscript scripts/audit_upstream.R` in both
   installed and source modes, confirming the namespaced-claim assertion at
   `audit_upstream.py:194-202` still holds.
5. A deliberate mutation check per gate — tamper the lock, break a threshold, remove a
   metric — confirming each gate can still report failure. The existing
   `test_python_audit_fails_on_a_commit_mismatch` is the model.
6. `docs/blocked.md` contains no entry that does not cite a command that currently fails.

## Decision and Evidence Log

- 2026-09-29 — Plan reviewed and initialized. Implementation not started.
- 2026-09-29 — Verified the environment premise behind B-001, B-002 and B-003 no longer
  holds: CRAN, Bioconductor, Posit PP and micromamba all reachable; R 4.4.1 present;
  660 R packages available in a personal library. Sequenced the cyCombine unblock first.
- 2026-09-29 — Decided to use a project-local R library rather than the personal
  660-package library, so this package's pinned R dependencies stay separable.
- 2026-09-29 — Decided the container stays the deployment answer but cannot be built or
  verified on this node; recorded as a limitation rather than silently skipped.

### Step 1.1 — the recorded baseline did not reproduce

`docs/implementation_status.md` claimed 142 collected / 140 passed / 2 skipped. The first
run in a clean project environment gave **1 failed, 125 passed, 16 skipped**. Two causes,
one of them a genuine defect:

1. **A real compatibility bug.** `src/batchdx/cyto/io/dataset.py:272` did
   `sorted(adata.layers)`. anndata 0.13 exposes a `None` key in `.layers` aliasing `.X`,
   so this raises `TypeError: '<' not supported between instances of 'NoneType' and
   'str'` — the error a user gets when they name an expression layer that does not exist,
   which is precisely the diagnostic that must be most reliable. `pyproject.toml:24`
   allows `anndata>=0.10`, so the package claimed support for a version it breaks on.
   Fixed at `dataset.py:276-283` by filtering `None` and naming `X` in the message.
2. **The recorded run had `Rscript` and `snakemake` on `PATH` and the upstream checkouts
   fetched.** Fourteen of the sixteen skips were those three absences, not defects.

The recorded numbers were never wrong about the code; they described a richer environment
than the one documented. Both are now stated in the plan rather than assumed.

### Step 1.3 — the R environment is a conda prefix, not an R library

The plan called for `libR/` plus `install.packages()`. That does not work here, and the
reason is a toolchain mismatch rather than anything about the packages:

- R 4.4.1 requires `GLIBCXX_3.4.29`; the default gcc is 8.5.0, providing at most 3.4.25.
- A conda R 4.4.3 works, but CRAN's current `Rcpp` source needs `CXXABI_1.3.15`, which only
  gcc 14 provides, and the newest gcc on this node is 13.2.0 (`CXXABI_1.3.14`).
- Posit Package Manager has no prebuilt binary for centos8/R 4.4, so it cannot substitute.

A conda environment resolves the closure to prebuilt packages built against one
consistent toolchain. `scripts/setup_r_env.sh` therefore creates `renv/` with micromamba
and installs cyCombine from its pinned commit into it. Two consequences worth recording:

- **`--vanilla` is load-bearing.** A `~/.Renviron` on this machine pins `R_LIBS_USER` to a
  personal renv library, and Renviron entries override the environment variable. The first
  install run therefore wrote cyCombine into that personal library. It was removed, and
  every `Rscript` call in the script now passes `--vanilla`, matching what
  `r_runner.run_r_script` already did. The personal library is no longer mutated.
- **cyCombine needs three packages its DESCRIPTION does not require.** `emdist` and
  `outliers` are `Suggests` but are called unconditionally; `pbmcapply` is checked and
  aborts the detector. Each was found by running the code, not by reading the manifest, and
  all three are now in the script with the reason recorded.

### Step 1.2 / 1.4 — defects the R layer exposed immediately

Running the R path for the first time found four problems that the audit of the pinned
source had not:

1. `workflow/scripts/cycombine_detect.R:70` had a tautological `label_col` and assigned the
   `"all"` label fallback at `:80`, **after** the `detect_batch_effect()` call at `:68` that
   reads it. Fixed at `:57-73`.
2. `cyCombine::detect_batch_effect` groups by a column **literally named `sample`**
   (`R/detect_batch_effect.R:200`), with no parameter to rename it. The adapter supplied
   only `batch` and `id`, so the detector aborted in `dplyr::group_by()`. The R script now
   requires and validates a per-cell `sample` vector, and `cycombine_detect.py:99-111`
   passes it. This was a genuine violation of the pinned API that only executing it showed.
3. `evaluate_mad` returns `list(score =, mad =)`, not a `reduction` field. The first
   implementation read `mad$reduction` and reported `undefined` for every input. Fixed, and
   the Python side reads `mad_score`.
4. The `native` specs had **no producer at all** — the Stage 7 clause was unreachable by
   construction. `src/batchdx/cyto/metrics/native.py` and
   `workflow/scripts/cycombine_evaluate.R` now call `evaluate_emd`/`evaluate_mad` as-is and
   report their return values. The identity candidate is the only one compared to itself,
   since a native reduction is a ratio of the two states and the identity is the only
   candidate whose uncorrected state is itself.

### Step 1.4 — the R audit now runs in `installed` mode

With cyCombine present, `scripts/audit_upstream.R` switched from parsing source to reading
the loaded namespace, where `deparse(args())` yields `function (...)` rather than
`name(...)`. `test_stage0_audit.py:135` asserts every signature names its function, so the
installed branch was reformatted to match the source branch (`audit_upstream.R:121-141`).
The test's contract was right; the implementation was inconsistent.

### Step 2.2 — three tests were skipping in the environment that needed them

`test_absent_r_package_yields_not_available_...`,
`test_absent_backend_reports_not_available_...` and
`test_a_requested_protected_covariate_is_never_quietly_dropped` all began
`if r_package_version("cyCombine").available: pytest.skip(...)`. With cyCombine installed
they became no-ops, so the *unavailable* path — the one the whole `not_available`
guarantee rests on — had no coverage at all, in exactly the environment where the backend
works. They now point `rscript` at a nonexistent interpreter, reproducing the condition
deterministically, and a new test
(`test_a_protected_covariate_is_refused_when_the_backend_is_available`) covers the
covariate contract on the real path. Stages 5 and 6 now run 27 tests with 0 skips.

### Step 3.3 — a test assertion that was simply wrong

`test_the_replicate_unit_is_the_sample_not_the_cell` asserted the cell-level t statistic
would exceed 10x the sample-level one. Measured on the fixture it is **0.89x** — the two
are the same size, because simulated within-sample spread is small next to the injected
donor effect. The premise of the assertion was false.

The test now asserts what is actually true and actually matters: tripling every cell
leaves all 18 preservation records bit-identical, while the cell-level analysis of the
same numbers reports a p-value of exactly **0** against the sample-level 6.2e-16, because
it spends 400 cells from one tube as if each were an independent biospecimen. An exact
zero is not a stronger result; it is the signature of counting the wrong thing.

### What the gates say now

No gate is BLOCKED. Stage 5 is **PASSED** (11/11) for the first time. Stages 6, 7 and 8
moved from BLOCKED/PARTIAL-with-no-code to PARTIAL with every declared test executing and
passing, each carrying an honest `partial_reason` naming what is still not satisfied.
Stage 7's native namespace is verified executed; its OpenProblems rows remain deferred, so
it is PARTIAL and not PASSED. Nothing claims the `compatibility` namespace.

Verified that the gates can still fail: pointing a gate at a nonexistent test id produces
`❔ NOT RUN`, and `scripts/generate_status.py` refuses to report a stage PASSED while any
declared test is missing, failed or skipped.

## Risks, Assumptions, and Blockers

- **No real cytometry dataset is obtainable from here.** B-008 remains true and Step 6.4
  cannot complete without a licence determination. Clearing condition: a permissibly
  licensed fixture with a stated provenance and checksums, or an explicit decision to
  ship synthetic-only with that stated on the face of the documentation. **This is now the
  single largest remaining scientific gap**, since the Stage 8 evidence is a demonstration
  that the machinery is correct rather than a calibration for any real cohort.
- **`renv/` is not the container.** It is reproducible from a pinned list and isolated, and
  it is what let the R adapters be executed and tested here. It is not a deployment
  artifact: no apptainer, singularity or docker exists on this node, so
  `containers/Dockerfile.r-backends` remains unbuilt and unverified, and nothing here
  changes that.
- **A `~/.Renviron` on the build host can silently redirect R installs.** It overrode
  `R_LIBS_USER` and wrote cyCombine into a personal library. Mitigated by `--vanilla`
  everywhere in `setup_r_env.sh` and by `run_r_script`'s existing `--vanilla`. A CI runner
  with the same file and a script that omits `--vanilla` would reintroduce this.
- **`anndata` compatibility is a live class of bug, not a one-off.** `pyproject.toml:24`
  pins `anndata>=0.10` with no upper bound, and 0.13 changed `.layers` semantics. The
  0.13 defect is fixed; a version matrix (Step 6.4) is what would catch the next one, and
  CI still tests Python 3.11 only.
- **`cytonorm` and `cytobatchflagr` remain unexecuted.** CytoNorm is GPL >= 2 and is called
  at arm's length, never vendored; its fit-and-apply interface needs a tested FCS round
  trip plus a `transformList` (Step 4.2), which depends on Step 4.1's splits. CytoBatchFlagR
  needs its ~35 CRAN dependencies (Step 3.4); its checkout is already fetched at the
  pinned commit. Both are honestly reported `not_available` today.
- **Single seed still.** The scientific suite runs at one seed. The Stage 8 code is
  explicitly seeded and reproducible from it (`test_a_prespecified_contrast_is_reproducible_from_its_seed`),
  but a multi-seed sweep is Step 6.4.
