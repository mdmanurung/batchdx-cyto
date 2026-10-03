# Release checklist — batchdx-cyto 0.1.0

The procedure for cutting a release, and the state of this one. Every box states the
command that establishes it. **A box that cannot be checked on the machine doing the
release is not ticked**; it is left open with the reason, which is the only thing that
keeps this document worth reading.

Target: `0.1.0` · Development status: **3 - Alpha** · Date of this checklist: 2026-10-03

Alpha remains the correct status. The correction backends and CytoBatchFlagR detector
have executed, but real-cohort calibration, container build and cluster execution have
not been verified.

---

## 1. Version is consistent in all three places

They move together; `tests/unit/test_release_metadata.py` fails if they drift.

- [x] `pyproject.toml:14` — `version = "0.1.0"`
- [x] `src/batchdx/cyto/_version.py:3` — `__version__ = "0.1.0"`
- [x] `CITATION.cff:15` — `version: 0.1.0`
- [x] `CHANGELOG.md` `[Unreleased]` section written from the actual tree

```bash
grep -n 'version' pyproject.toml | head -3
grep -n '__version__' src/batchdx/cyto/_version.py
grep -n 'version:' CITATION.cff | head -2
```

## 2. The tree is clean and buildable from a fresh clone

- [x] `LICENSE` exists and is the MIT text. It was **untracked** until this release; a
      package that declares MIT in two files and ships no licence text cannot be
      distributed at all.
- [x] `pyproject.toml` uses PEP 639 (`license = "MIT"`, `license-files = ["LICENSE"]`).
- [x] `requires = ["hatchling>=1.27"]`. Load-bearing: built with 1.26.3 the wheel succeeds
      but emits the deprecated `License:` field and drops `License-File`, so the
      published metadata would depend on whichever hatchling the resolver picked.
- [x] `THIRD_PARTY_NOTICES.md` present and accurate.
- [x] `CITATION.cff` present and valid.
- [x] Required source and the R setup script are committed in the A1 commit series.
- [x] Disposable fresh clone at `/tmp/batchdx-fresh-20261003` built wheel and sdist,
      installed `.[dev]` into a new Python 3.12 virtual environment, imported version
      `0.1.0`, fetched all five pinned upstream checkouts, and passed the Python audit.

```bash
git status --porcelain | grep '^??' | wc -l    # must become 0
```

## 3. The test suite is green, and the skips are explained

- [ ] `pytest -q` — a complete current run is still needed. The non-slow run on
      2026-10-03 had 355 passed, 1 conditional skip, 23 slow tests deselected.
- [x] No unexplained skips. The skips that exist are gated on an absent package, and each
      names the package.
- [x] `ruff check src tests scripts` — clean.
- [x] Both audit scripts exit 0: `python scripts/audit_upstream.py` and
      `Rscript --vanilla scripts/audit_upstream.R`.
- [x] `python scripts/generate_status.py` regenerated `docs/implementation_status.md`
      from the completed non-slow report; Stages 1 and 11 show `NOT RUN` because their
      workflow tests were deselected. CI staleness checking remains unverified here.
- [x] **The R exchange suite is green locally.** R 4.4.3 ran 101 assertions with
      0 failures, 0 errors and 0 skips on 2026-10-03. The same suite passed in
      draft PR CI run 37102250467.

```bash
Rscript --vanilla tests/run_r_tests.R
.venv/bin/pytest -q -m 'not slow'                         # 355 passed, 1 skipped, 23 deselected
.venv/bin/ruff check src tests scripts
```

## 4. Every gate test exists and actually ran

- [x] `docs/implementation_status.md` is generated from a real report, not written by
      hand.
- [x] `scripts/generate_status.py:88-90` refuses to report `PASSED` while any declared
      test is missing, failed or skipped. Verified by pointing a gate at a nonexistent
      test: the status becomes `NOT RUN`, never `PASSED`.
- [x] **Stage 12 has executed test identifiers.** `scripts/gates.yaml` now registers
      release metadata and seed checks; its remaining limits are stated as `PARTIAL`.
- [x] **Every declared gate test exists.** `pytest --collect-only -q` found every
      identifier in `scripts/gates.yaml`; execution of Stages 1 and 11 remains open.

## 5. The release claims match what was measured

- [x] `docs/compatibility_matrix.md` reconciles tested against declared. **3.12 is the
      only `verified` Python**; 3.10, 3.11 and 3.13 are `declared-only`, and the source
      parsing under their grammars is labelled as the narrow claim it is.
- [x] No row anywhere claims a version was exercised when it was not.
- [x] `docs/validation_results.md` states what the scientific suite does **not**
      establish, including that the multi-seed sweep covers the determinism and stability
      claims but **not** a per-seed re-run of the whole suite.
- [ ] **Complete CI is not yet green.** Draft PR #2 CI run 37102250467 reached the
      full suite: 375 passed, 2 skipped, 2 failed. Both failures were interrupted-run
      byte comparisons caused by elapsed walltime embedded in result records. The
      focused fix is under CI rerun; its status-document comparison has not passed.
- [ ] **The container is unbuilt.** No docker, podman, apptainer or singularity on this
      node, so the image cannot be built, run or digest-pinned. `containers/` is shipped
      with a `NOT BUILT` banner and that is the honest state.

## 6. Known defects, stated rather than silently shipped

Each is on the record in `CHANGELOG.md` and `docs/blocked.md`. A release is allowed to
ship a known defect; it is not allowed to ship an unrecorded one.

- [x] **`binmat.R:49` — length-1 sidecar identifier fields now remain arrays.** A one-row
      or one-column matrix previously unboxed an identifier to a JSON string.
      `I(as.character(...))` now preserves `row_ids` and `col_names` as arrays.
      `tests/testthat/test-binmat.R` covers degenerate shapes and Python round trips.

- [x] **`io/dataset.py:255` — the anndata 0.13 `None`-key defect, second instance.**
      Both error paths exclude the `None` alias before sorting names.
      `test_missing_expression_layer_with_none_alias_names_layers` covers the path.

- [ ] **No dependency carries an upper bound**, and all sixteen lower bounds are
      untested. The recommended ceilings are tabulated in
      `docs/compatibility_matrix.md` § *Dependency bounds*. `anndata` is the one that has
      already bitten.

## 7. Documentation is consistent with the code

- [x] `docs/README.md` indexes every document, including `docs/plans/`.
- [x] `CHANGELOG.md` and `RELEASE.md` exist. Neither did.
- [x] The synthetic-only decision is on the face of `examples/README.md`,
      `docs/validation_results.md` and `docs/compatibility_matrix.md`, not only in
      `docs/blocked.md`. `README.md` is the one remaining surface and needs the text in
      §8.
- [x] The four tracked documents corrected for staleness are regenerated and
      diff-clean: `THIRD_PARTY_NOTICES.md`, `workflow/envs/README.md`,
      `config/upstream_lock.yaml`, and `docs/implementation_status.md`.
- [x] **`README.md` states the synthetic-only boundary** beside the package description.

## 8. Central changes required before this release is honest

Items 2–5 are complete locally. Item 1 remains necessary for a fresh clone.

1. **Done locally:** commit the required source and release files, including
   `scripts/setup_r_env.sh` and `workflow/envs/{python,openproblems}.yaml`.
2. **Done:** preserve R sidecar identifier arrays in `workflow/scripts/lib/binmat.R`.
3. **Done:** filter the AnnData `None` layer alias in `src/batchdx/cyto/io/dataset.py`.
4. **Done:** register Stage 12 tests in `scripts/gates.yaml`.
5. **Done:** state the synthetic-only scope in `README.md`.

6. **Configured, awaiting CI execution:** the R suite uses the job's Python
   interpreter and the isolated R setup installs `testthat`.
7. **Configured, awaiting CI execution:** the primary CI job and workflow environment
   both use Python 3.12. Python 3.10, 3.11 and 3.13 remain declared-only.
8. **Set upper bounds on the sixteen dependencies** per the table in
   `docs/compatibility_matrix.md`, and add a lowest-supported-dependency CI job so the
   lower bounds become measurements rather than promises.

## 9. Signing off

Cutting this release is defensible **with** §8 items 1–5 done. Items 6–8 are follow-up
work that does not block a `0.1.0` marked Alpha, provided each stays on the record —
which §6 and `docs/blocked.md` exist to guarantee.

Remaining Alpha release blocks are a complete current suite and a green CI run. The tag
and publication require final approval.
