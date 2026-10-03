# Changelog

All notable changes to `batchdx-cyto` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**How to read the `Unreleased` section.** This repository's history is short and
unusual: the five tracked commits are all dated 2026-09-14/15, and every change made since
— the report-correctness tranche, preprocessing, the R environment, the R test suite and
this release-readiness tranche — is **still uncommitted**. `Unreleased` below is therefore
written from the working tree through 2026-10-03, not from commit messages, and the
commit that eventually carries it will be the first since `9850c33`. Where a claim here
could be checked, it names the command that checks it; `docs/compatibility_matrix.md` is
the authority on which claims are `verified` and which are `declared-only`.

The project version is `0.1.0`, mirrored in `pyproject.toml:14`,
`src/batchdx/cyto/_version.py:3` and `CITATION.cff:15`. All three move together.

## [Unreleased] — 2026-10-03

The first release candidate. `Development Status :: 3 - Alpha` is accurate: several
gates are `PARTIAL`, and real-cohort calibration and deployment remain unverified.

### Added

- **SPiSCy-style preprocessing** (`src/batchdx/cyto/preprocess/`, new): 2-D rectangular
  and polygon gating over `mk1`/`mk2`, the per-panel asinh transform, a PeacoQC 1.16.0
  pass, and a `removed_events` writer. Covered by `tests/unit/test_preprocessing.py` (30)
  and `tests/integration/test_peacoqc_backend.py` (17). PeacoQC executes for real: a
  fixture with a known contaminated block has every injected event removed and
  attributed, on both the isolation-tree and MAD paths.
- **OpenProblems metric components** (`src/batchdx/cyto/metrics/openproblems.py`, new):
  six of the nine upstream components registered in the `compatibility` namespace with
  `upstream_verified=True`, each proven by a test that imports the **pinned** helper file
  and compares outputs. `scripts/audit_upstream.py` fails the build if a `compatibility`
  row ever appears without a recorded verification.
- **Train/tune/validate splits** (`src/batchdx/cyto/sampling/splits.py`, new): a
  three-way split over samples, with held-out batch-pair checking so a validation split
  cannot be credited with a contrast it cannot supply.
- **R-side test suite** (`tests/testthat/`, new): `test-binmat.R` holds
  `workflow/scripts/lib/binmat.R` and `src/batchdx/cyto/io/binary.py` to a contract
  stated by neither of them, so a row-major/column-major or byte-order divergence cannot
  hide behind a matching divergence. Driven by `tests/run_r_tests.R`.
- **Multi-seed checks** (`tests/scientific/test_multiseed.py`, new): 25 tests over five
  seeds, covering the determinism and stability claims rather than re-running the whole
  suite per seed. Caught a real defect — see *Fixed*.
- **Conda environment recipes** (`workflow/envs/python.yaml`,
  `workflow/envs/openproblems.yaml`): the versions this build was developed against,
  recorded as numbers rather than as an opinion.
- **MIT `LICENSE`** and PEP 639 `license = "MIT"` / `license-files = ["LICENSE"]`. The
  project declared MIT in `pyproject.toml` and `CITATION.cff` with no licence text
  shipping, which blocks any distribution. The metadata form needs `hatchling>=1.27`;
  built with 1.26.3 the wheel silently degrades to the deprecated `License:` field and
  drops `License-File`.
- **`CHANGELOG.md`** and **`RELEASE.md`** (this file and the next), neither of which
  existed.
- **Three R adapter scripts** (`workflow/scripts/{cycombine_evaluate,cytonorm_correct,peacoqc_qc}.R`)
  and **`scripts/setup_r_env.sh`**, the script that builds the isolated R environment and
  which eleven tracked files reference.
- The `docs/plans/` roadmap and tranche records.

### Fixed

- **The report silently dropped successful native detectors.**
  `report.html.j2:200` gated the diagnostics section on
  `{% if d.detector in ('acquisition_drift', 'composition') %}`, so a *successful* run of
  `cycombine_express`, `cycombine_cluster` or `cytobatchflagr` appeared nowhere; only a
  non-`ok` one surfaced, through the failures table. Every detector is now rendered, and
  so are the nine `DetectorResult` fields that were fetched and never used.
- **Unflagged assessments were discarded.** `report.html.j2:204` did `selectattr('flag')`,
  which dropped every assessment with `flag=None` — including
  `max_block_step_over_mad`, the drift module's own headline discontinuity statistic, and
  every per-sample `population_frequency` row. They now render in a separate labelled
  table, and the two silent truncations state "showing N of M".
- **The two hard-coded report caps are declared** rather than silent:
  `MAX_EXPLORER_MARKERS = 8` and `markers[:3]`.
- **`_metric_rows` accepted a `redactor` and discarded it**, projecting
  `MetricRecord.comparison_id` and `DetectorAssessment.sample_id` raw.
- **`n_donors` was set on the EMD family only**, so the report's `n_donors` column was
  blank for MAD, frequency and joint-structure.
- **An exploratory override was invisible.** `pipeline.py:242` was the repository's only
  bare `pass`: with `allow_exploratory_override` and an overridable ineligibility the
  backend ran and nothing was recorded. The override is now recorded and passed through
  to `DecisionResult.overrides`, which had no producer.
- **Stage 8 uncertainty never reached the report.**
  `MetricRecord.uncertainty_low/high/kind` are populated by every preservation record and
  projected by nothing, not the table and not the CSV.
- **Seed was threaded from one source.** `materialise_example(..., seed=7)` generated
  data with seed 7 and wrote a configuration declaring `seed: 473`; every stage reads
  `config.project.seed`, so the shipped run record described a run that was not the run.
  The default now lives once, as `config.DEFAULT_SEED`, and `ProjectConfig.seed`,
  `SimulationSpec.seed`, the three named scenarios, `materialise_example` and
  `DetectorContext.seed` all derive from it.
- **`io/dataset.py:272` crashed on anndata 0.13**, which exposes a `None` key in
  `.layers`. Fixed at `:361-364` by filtering `if name is not None`.
- **Three stale tracked documents**: `THIRD_PARTY_NOTICES.md:30` claimed no
  implementation was extracted (six components are), `workflow/envs/README.md` said one
  environment file where there are two, and `config/upstream_lock.yaml:120-122` named six
  `dependency_risk` entries for a 29-package enforced set.
- **`workflow/scripts/lib/binmat.R` was swallowed by `.gitignore`**, whose inherited
  Python template carries a bare `lib/`. A fresh clone got a workflow whose R adapters
  could not load. Fixed by an explicit negation.
- **One-element R sidecar identifiers unboxed to JSON strings.** `write_binmat()` now
  preserves `row_ids` and `col_names` as arrays; the R exchange suite passes all 101
  assertions on R 4.5.1.
- **A second AnnData 0.13 `None` layer alias path raised `TypeError`.** The missing-layer
  error now filters that alias and has a regression check.

### Changed

- `preprocessing.transformation` no longer accepts `"logicle"`. It was a declaration with
  no capability behind it: this build implements no logicle, `io/export.py` has no
  tested inverse for one, and the pinned reference implementation contains none either
  (`grep -i logicle` over `.upstream/spiscy` returns zero matches). Accepting it validated
  a value and bought nothing.
- `preprocessing.gating` and `preprocessing.peacoqc` moved out of `UNSUPPORTED_KEYS`
  into real fields, read as *paths* to a gate table and a QC settings file, on the same
  change that gave them readers. `preprocessing.time_channel` **remains refused**,
  because nothing makes the drift detector read a configured channel.
- Version `0.1.0.dev0` → `0.1.0`; `Development Status` and the `Programming Language ::
  Python :: 3.10/3.11/3.12/3.13` and `:: R` classifiers added.
- An explicit `import WorkflowError` in `workflow/Snakefile:19`. Snakemake injects its
  own module globals, so the previous code worked; the import replaces an undeclared
  reliance on that injection and is labelled in place as such.

### Known issues, stated rather than fixed

These are open on purpose. `docs/blocked.md` carries each with the exact command that
still fails.

- **No real cytometry data, and this is a scope decision.** No dataset with reference or
  replicate structure was obtainable under the current egress policy, so this package
  ships synthetic fixtures only. B-008 closes as a stated boundary, not a capability.
  Every `verified` row in `docs/compatibility_matrix.md` was measured against synthetic
  input.
- **CytoBatchFlagR flags remain uncalibrated.** The detector now executes, but its
  1.5-IQR fence has a high per-test false-positive rate at realistic batch counts.
- **All sixteen dependency floors are untested**, and no dependency carries an upper
  bound. The recommended ceilings are tabulated in `docs/compatibility_matrix.md`.
- **Python 3.10, 3.11 and 3.13 are declared-only.** 3.12 is the tested version. The
  source parses under all four grammars, which is the only claim that can be made without
  installing a stack.
- **The container has never been built.** No docker, podman, apptainer or singularity on
  this node.
- **The SLURM profile has never run on a cluster**, and still carries two `REPLACE_ME`
  values.

## [0.1.0.dev0] — 2026-09-15

The five commits that exist.

- `fc06be5` 2026-09-14 — Initial commit (the cookiecutter template this repository also
  carries; 53 files).
- `3749f94` 2026-09-15 — **Add batchdx-cyto**: cytometry batch diagnosis, correction and
  evaluation. 102 files, +14,740. The whole package: the 12-stage pipeline, the strict
  configuration schema, the three shipped examples, the Snakemake workflow, the adapters
  for cyCombine, CytoNorm and CytoBatchFlagR, the OpenProblems metric integrations, the
  offline HTML report, and the gate/audit scripts.
- `0668462` 2026-09-15 — Add CI for the package and scope the template's test task.
- `08776c8` 2026-09-15 — Fix the template test job: guard the shared conftest, rename the
  new CI job.
- `89bddbf` 2026-09-15 — Ship `workflow/scripts/lib/binmat.R`, which gitignore was
  swallowing.
- `9850c33` 2026-09-15 — Merge PR #1.

## A note on what is *not* in this changelog

Stage statuses. `docs/implementation_status.md` is generated from an executed test report
and is the authority on which gates passed, partially passed or were deferred; this file
does not restate it, because a hand-maintained copy of a generated document drifts.
Compatibility claims live in `docs/compatibility_matrix.md`, and every blocked item lives
in `docs/blocked.md` with the command that still fails it.

[Unreleased]: https://github.com/mdmanurung/batchdx-cyto/compare/v0.1.0...HEAD
[0.1.0.dev0]: https://github.com/mdmanurung/batchdx-cyto/releases/tag/v0.1.0.dev0
