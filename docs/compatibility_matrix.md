# Compatibility matrix

What is actually known about the versions and backends this package runs on, and — for
every row — the command that establishes it.

**The rule this file follows.** A cell says `verified` only if a command was run on this
checkout and produced the stated result. It says `declared-only` when a tracked file
asserts it but nothing here executed it. It says `untested` when there is neither. No row
is `verified` because a version number looks plausible. Where a number was resolved from a
remote index over the network, that is stated in the row, because "resolved" and "ran"
are different claims.

Measured on 2026-09-30 at commit `9850c33` plus uncommitted work in the tree.

**Environment of record.** Python 3.12.13, R 4.4.3 (`renv/bin/R`), `.venv` for Python and
`renv/lib/R/library` for R. Every `verified` cell below was measured here.

---

## Python

### Tested versus declared, reconciled

Four tracked files made four different statements about which Python this package runs
on, and the disagreement was itself the defect:

| Source | Said |
|---|---|
| `pyproject.toml` `requires-python` | `>=3.10` — a floor, i.e. everything at or above 3.10 is *permitted* |
| `pyproject.toml` classifiers | 3.10, 3.11, 3.12, 3.13 — four versions *advertised* |
| `.github/workflows/batchdx-cyto-tests.yaml:29` | `python-version: "3.12"` — configured for CI, with no completed runner job yet |
| `workflow/envs/python.yaml:29` | `python=3.12` — one version *pinned for the workflow* |

A floor, a classifier list, a CI pin and an environment pin are four different claims,
and only two of them are ever executed. The reconciliation is deliberately boring:
**3.12 is the tested version; the other three are declared.** Nothing is widened to make
the table look better, and no version is marked `verified` because a number next to it
looks plausible.

| Version | Status | What was actually run | How to verify |
|---|---|---|---|
| **3.12** | `verified` | the suite in two runs: `.venv/bin/python -V` → `Python 3.12.13`; `pytest -q tests/unit` → 123 passed in 25s; `pytest -q tests/scientific` → 95 passed / 1 skipped in 1m47s; together 243 passed / 1 skipped in 2m13s, of which 25 are `tests/scientific/test_multiseed.py` | `python -V`; `.venv/bin/pytest -q tests/unit tests/scientific` |
| 3.11 | `declared-only` | **nothing.** No environment was built and no test was run. A grammar-level parse of all 85 `.py` files under 3.11 succeeds, which rules out a syntax regression and establishes nothing else | `python3.11 -m venv .venv311 && .venv311/bin/pip install -e '.[dev]' && .venv311/bin/pytest -q` |
| 3.10 | `declared-only` | **nothing.** Same as 3.11: grammar-level parse only. This is the `requires-python` floor, so a 3.10 incompatibility is a *packaging* bug, not an untested feature | as above, with `python3.10` |
| 3.13 | `declared-only` | **nothing.** Grammar-level parse only. A trove classifier with no environment behind it is a claim; it is labelled as one | as above, with `python3.13` |

**The grammar-level check, and exactly what it is.** All 85 `.py` files under `src/`,
`tests/` and `scripts/` were parsed with `ast.parse(..., feature_version=(3, N))` for
N in 10, 11, 12, 13; all 85 parse under all four. That establishes that the source uses
no syntax newer than 3.10, which is a real constraint on the `requires-python` floor and
is the only thing about 3.10 that can be checked without installing a stack. It
establishes **nothing** about the standard library surface, C-extension ABI, or
dependency resolution, all of which differ across these versions.

```bash
python - <<'PY'
import ast, pathlib
files = list(pathlib.Path("src").rglob("*.py")) + list(pathlib.Path("tests").rglob("*.py"))
for minor in (10, 11, 12, 13):
    for p in files:
        ast.parse(p.read_text(), filename=str(p), feature_version=(3, minor))
print("all", len(files), "files parse under 3.10-3.13")
PY
```

**There is no CI matrix across versions.** The workflow names Python 3.12, matching the
workflow environment. Three of the four rows above remain claims, not measurements,
until a matrix runs. The Python 3.12 draft-PR job reached the complete suite on a
GitHub runner in [run 37102250467](https://github.com/mdmanurung/batchdx-cyto/actions/runs/37102250467):
375 passed, 2 skipped and 2 interrupted-run byte comparisons failed. This is executed
3.12 evidence, not yet a green compatibility result.

### Core dependency versions that were resolved and run

These are the versions in `.venv`, and they are exactly the pins in
`workflow/envs/python.yaml`, so the environment used here is the one that file describes.

| Package | Version | How to verify |
|---|---|---|
| numpy | 2.5.3 | `python -c "import importlib.metadata as m; print(m.version('numpy'))"` |
| pandas | 3.0.6 | as above, `pandas` |
| scipy | 1.18.1 | as above, `scipy` |
| anndata | 0.13.4 | as above, `anndata` |
| h5py | 3.16.0 | as above, `h5py` |
| pydantic | 2.13.5 | as above, `pydantic` |
| typer | 0.27.2 | as above, `typer` |
| jinja2 | 3.1.6 | as above, `jinja2` |
| PyYAML | 6.0.3 | as above, `PyYAML` |
| flowio | 1.4.0 | as above, `flowio` |
| pyarrow | 25.0.1 | as above, `pyarrow` |
| matplotlib | 3.11.2 | as above, `matplotlib` |
| snakemake | 9.27.0 | as above, `snakemake` |
| scikit-learn (`metrics` extra) | 1.9.1 | as above, `scikit-learn` |
| scib-metrics (`metrics` extra) | 0.6.1 | as above, `scib-metrics` |
| umap-learn (`metrics` extra) | 0.5.12 | as above, `umap-learn` |

The floors in `pyproject.toml` are deliberately much wider than this table, and
**every one of those lower bounds is untested** — see the next section, which is the
honest reading of that gap rather than a defence of it.

### Dependency bounds: lower bounds only, and what that costs

All twelve core dependencies and all three `metrics`-extra dependencies carry a lower
bound and no upper bound. That is a real gap, not a stylistic one, and the
`anndata>=0.10` line has already produced a user-visible failure:

* **What happened.** anndata 0.13 added a `None` key to `.layers` that aliases `.X`. Any
  code that iterates or sorts the layer view hits it. `io/dataset.py:361-364` filters
  `if name is not None` and is correct; **`io/dataset.py:255` does not**, and still does
  `sorted(adata.layers)` in an error message. With anndata 0.13.4 a file carrying both a
  `None` key and a named layer raises
  `TypeError: '<' not supported between instances of 'str' and 'NoneType'` where the
  intended `DatasetValidationError` should name the available layers. Reproduce it with:
  ```python
  import anndata as ad, numpy as np
  a = ad.AnnData(X=np.zeros((3, 2), dtype="float32")); a.layers["preprocessed"] = np.zeros((3, 2))
  sorted(a.layers)   # TypeError
  sorted(n for n in a.layers if n is not None)   # the fix already applied at :361-364
  ```
  A `TypeError` from a library's error path is the worst of both: the message is
  unreadable *and* the validation it was written to perform never runs.
* **Why one fix does not close the class.** A lower bound cannot express "and not the
  version that changed this". Only an upper bound, or a version matrix, stops the next
  one. `pandas 3.0` and `numpy 2.x` are the same shape of risk and have not yet broken
  anything here, which is luck rather than a bound.

**Recommended bounds.** These are recommendations, not edits — `pyproject.toml` is not
this tranche's to change. Each ceiling is the next major version above the one row in
the table above that was actually executed, so the bound is "everything up to and
including what we ran", which is the most that can honestly be promised.

| Package | Now | Recommended | Ceiling rationale |
|---|---|---|---|
| numpy | `>=1.24` | `>=1.24,<3` | 2.x changed the ABI and scalar promotion; 3 is untested |
| pandas | `>=2.0` | `>=2.0,<4` | 3.0 is a major release (copy-on-write, removals); 4 is untested |
| scipy | `>=1.10` | `>=1.10,<2` | 1.x is tested; 2 is untested |
| anndata | `>=0.10` | `>=0.10,<0.14` | **the live one**: 0.13 changed `.layers` semantics |
| h5py | `>=3.8` | `>=3.8,<4` | 3.x is tested; 4 is untested |
| pydantic | `>=2.5` | `>=2.5,<3` | 2.x is tested; 3 is untested |
| typer | `>=0.9` | `>=0.9,<1` | pre-1.0: the ceiling is 1.0, not 0.28 |
| jinja2 | `>=3.1` | `>=3.1,<4` | 3.x is tested; 4 is untested |
| pyyaml | `>=6.0` | `>=6.0,<7` | 6.x is tested; 7 is untested |
| flowio | `>=1.2` | `>=1.2,<2` | 1.x is tested; 2 is untested |
| pyarrow | `>=14.0` | `>=14.0,<26` | 25.0.1 is tested; 26 is untested |
| matplotlib | `>=3.7` | `>=3.7,<4` | 3.x is tested; 4 is untested |
| snakemake (`workflow`) | `>=8.0` | `>=8.0,<10` | 9.27.0 is tested; the 8→9 jump moved APIs |
| scikit-learn (`metrics`) | `>=1.3` | `>=1.3,<2` | 1.9.1 is tested; 2 is untested |
| scib-metrics (`metrics`) | `>=0.5` | `>=0.5,<0.7` | 0.6.1 is tested; pre-1.0 ceiling is the minor |
| umap-learn (`metrics`) | `>=0.5` | `>=0.5,<0.6` | 0.5.12 is tested; pre-1.0 ceiling is the minor |

**The lower bounds are the larger half of the problem and ceilings do not touch them.**
`numpy>=1.24` permits 1.24, which has never been run here; `pandas>=2.0` permits 2.0;
and so on down the table. A ceiling stops the *next* breaking release; only a
lowest-supported-dependency job — a second CI job that resolves every dependency to its
floor and runs the suite — turns those promises into measurements. Neither exists today,
and no row in this file claims otherwise.

### Build backend

| Item | Status | Evidence | How to verify |
|---|---|---|---|
| `hatchling` >= 1.27 required | `verified` | built with 1.26.3 and with 1.27.0. 1.27.0 emits `License-Expression: MIT` + `License-File: LICENSE`; 1.26.3 succeeds but silently emits the deprecated `License: MIT` and drops `License-File` | `python -m build --no-isolation` against each version; compare `METADATA` in the wheel |

`hatchling` is a build-time requirement only and is **not** installed in `.venv`, so
nothing in the test suite pins it. `requires = ["hatchling>=1.27"]` is what stops a build
silently degrading to the old metadata field.

---

## R

| Version | Status | Evidence | How to verify |
|---|---|---|---|
| 4.4.3 | `verified` | `renv/bin/R`; all 12 R-marked tests pass under it | `R --version`; `pytest -m "requires_r or requires_cycombine or requires_cytonorm"` |
| 4.4.1 | `declared-only`, and **never built** | the container base image. No docker, podman, apptainer, singularity, skopeo or crane exists on this node | build `containers/Dockerfile.r-backends` somewhere with a runtime |
| any other 4.x | `untested` | nothing declares or runs one | — |

R 4.4 is what makes Bioconductor 3.20 the relevant release; that is the release the
container's `BiocManager::install()` will resolve against.

### R-side test coverage

R is not an optional garnish here: **eight** R files under `scripts/` and
`workflow/scripts/` do work that `batchdx.cyto.pipeline` drives, which is why
`Programming Language :: R` is declared in `pyproject.toml`. Until this tranche, one of
them (`scripts/audit_upstream.R`) was the only one CI executed.

| R file | Unit tests | Driven end to end by pytest |
|---|---|---|
| `workflow/scripts/lib/binmat.R` | **yes** — `tests/testthat/test-binmat.R`, 101 assertions | yes, `tests/unit/test_contracts.py:124-165` |
| `workflow/scripts/cycombine_correct.R` | no | yes — cyCombine correction runs at its pinned commit |
| `workflow/scripts/cycombine_detect.R` | no | yes — Stage 5 native-detector comparison |
| `workflow/scripts/cytobatchflagr_detect.R` | no | **no — never executed**, see B-006 below |
| `workflow/scripts/cycombine_evaluate.R` | no | no |
| `workflow/scripts/cytonorm_correct.R` | no | yes — CytoNorm backend, 4/4 |
| `workflow/scripts/peacoqc_qc.R` | no | yes — `tests/integration/test_peacoqc_backend.py` |
| `scripts/audit_upstream.R` | no | yes — CI, and the Stage 0 gate |

`binmat.R` is the one covered directly, and it is covered directly on purpose: its own
header calls it a hand-mirror of `io/binary.py`, and it is the only file in the
package where a byte-order or row-major/column-major divergence would be silent *and*
would corrupt every adapter. `tests/testthat/test-binmat.R` holds both implementations to
a contract stated by neither of them — the expected bytes are produced by base R
`writeBin` against a `raw()` connection — so a divergence cannot hide behind a matching
divergence. The suite currently reports **94 passing and 7 failing**, all 7 from one
defect in `binmat.R:49` (`auto_unbox = TRUE` collapses a length-1 `row_ids`/`col_names`
to a JSON scalar, which Python's `list(...)` then iterates character by character). See
the tranche report for the one-line fix and the verified blast radius.

Run the R suite with:

```bash
Rscript --vanilla tests/run_r_tests.R          # exits non-zero on any failure
BATCHDX_PYTHON=.venv/bin/python Rscript --vanilla tests/run_r_tests.R   # force cross-language
```

---

## Native R backends

R is not an optional garnish here: eight tracked R files
(`scripts/audit_upstream.R`, `workflow/scripts/{cycombine_correct,cycombine_detect,cytobatchflagr_detect,cycombine_evaluate,cytonorm_correct,peacoqc_qc}.R`,
`workflow/scripts/lib/binmat.R`) do work that `batchdx.cyto.pipeline` drives, which is why
`Programming Language :: R` is declared in `pyproject.toml`. Per-file unit-test coverage is
tabulated under *R-side test coverage* above.

| Backend | Pinned commit | Library status | Status | How to verify |
|---|---|---|---|---|
| cyCombine | `fcd42dbb` | 0.3.0 installed, loads | `verified` — **executed**: 12 R-marked tests pass, none skipped | `pytest -m requires_cycombine`; `Rscript --vanilla -e 'library(cyCombine)'` |
| CytoNorm | `278f352d` | 2.0.12 installed, loads | `verified` — **executed**: `tests/integration/test_cytonorm_backend.py` passes 4/4 | `pytest tests/integration/test_cytonorm_backend.py`; `Rscript --vanilla -e 'library(CytoNorm)'` |
| CytoBatchFlagR | `ad21e3b1` | all 29 packages in its unconditional `library()` block installed | `verified` — detector executes on shifted and unshifted four-batch fixtures | `pytest tests/integration/test_cytobatchflagr_backend.py` |

### CytoBatchFlagR execution and limits

The enforced dependency set is the 29 unconditional `library()` calls at
`.upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R:6-34`. All 29 load in the
isolated R environment. The integration test calls the pinned source, observes the
injected flag, and checks the unshifted design. The earlier tests that stopped at an
eligibility guard remain useful refusal-path tests; the new integration tests reach R.

The 1.5-IQR fence has a high per-test false-positive rate at realistic batch counts,
so `flag=True` is an uncalibrated finding. The adapter also requires at least 100 cells
per control and marker to avoid treating a dropped control as a clean result. See
`docs/blocked.md` B-006 for the measured specificity limitation and rationale.

Check the executed path with:

```bash
.venv/bin/pytest -q tests/integration/test_cytobatchflagr_backend.py
```

---

## OpenProblems compatibility components

Six of the nine OpenProblems metric names are registered in this package's
`compatibility` namespace with `upstream_verified=True`. "Verified" here has one specific
meaning: a test imports the **pinned** Viash helper from the pinned checkout and compares
its output against ours under the tolerance below. Each of the six was executed on
2026-09-30 and passed.

| Component | Pinned helper | Tolerance asserted by the test | Status | Test |
|---|---|---|---|---|
| `emd_mean_ct_horiz` | `emd/helper.py` | `rel=1e-9` | `verified` | `test_the_pinned_horizontal_emd_matches_the_extracted_component` |
| `emd_mean_ct_vert` | `emd/helper.py` | `rel=1e-9` | `verified` | `test_the_pinned_vertical_emd_matches_the_extracted_component` |
| `ratio_consistent_peaks` | `ratio_consistent_peaks/helper.py` | `rtol=0, atol=1e-12` on density; exact array equality on the counts | `verified` | `test_the_pinned_peak_helpers_reproduce_the_extracted_ones` |
| `average_batch_r2_ct` | `average_batch_r2/helper.py` | `rel=1e-9, abs=1e-12` per marker and on the aggregate | `verified` | `test_the_pinned_batch_r2_helper_matches_the_extracted_component` |
| `iLisi` | `lisi/helper.py` | `rtol=0, atol=0` per cell; `rel=1e-12, abs=1e-15` on the normalisation | `verified` | `test_the_pinned_lisi_helpers_match_the_extracted_ones` |
| `cLisi` | `lisi/helper.py` | `rtol=0, atol=0` per cell; `rel=1e-12, abs=1e-15` on the normalisation | `verified` | `test_the_pinned_lisi_helpers_match_the_extracted_ones` |

All six node IDs above are in `tests/scientific/test_stage7_openproblems.py`. Three
supporting helpers are verified the same way and are worth knowing about because they are
what the deferred rows would build on: `test_the_pinned_bin_array_reproduces_the_extracted_one`,
`test_the_pinned_compute_emd_reproduces_the_extracted_distance`,
`test_the_pinned_wilcoxon_helper_matches_the_local_rank_statistic` and
`test_the_pinned_cohens_d_helper_matches_the_local_cohens_d`.

### The three that are *not* verified

| Component | Status | Obstacle | Where |
|---|---|---|---|
| `flowsom_mean_mapping_similarity` | `untested` | FlowSOM is Bioconductor-only and absent from `renv/`; the metric also needs three `uns['parameter_*']` benchmark keys this package does not define | `DEFERRED_METRICS` in `src/batchdx/cyto/metrics/registry.py`; `docs/blocked.md` B-007 |
| `functional_marker_preservation_wilcoxon` | `untested` | upstream reports a retention *ratio*, not the statistic this package computes; needs two integrated splits and an uncorrected baseline | as above |
| `functional_marker_preservation_cohens_d` | `untested` | upstream reports effect-size *drift* against an uncorrected baseline, not a per-pair d | as above |
| `cms` | `not present upstream` | named in the CHANGELOG, but `src/metrics/cms` does not exist at the pinned commit | `test_cms_stays_recorded_as_absent_at_the_pinned_commit` |

---

## Data: this package ships no real cytometry data, and that is the design

**Stated on the face of the documentation, not only in `docs/blocked.md`.** This package
ships **no real cytometry data and no real-data example**, because none is obtainable
under the current egress policy. That is a scope decision, taken deliberately, not a gap
waiting to be filled: no gate depends on real data, and no claim in this package is
about real data.

The consequence is a hard limit on everything above. Every `verified` row in this file
was measured against deterministic synthetic fixtures from `batchdx.cyto.simulate`. A
version being exercised here means *this code ran on this input*; it does not mean the
code has been exercised on anything a flow cytometrist would recognise. The specific
things real data would test and these fixtures do not are heavy tails, multimodality,
marker dropout, instrument artefacts and genuine spillover compensation.

| | |
|---|---|
| Real data included | **none** |
| Real-data example included | **none** |
| Origin of every fixture | `batchdx.cyto.simulate` — five populations, declared marker means, a shared within-population latent factor, and separately controllable `batch_shift` / `batch_spread` / `group_shift` / `drift` / frequency knobs |
| Obstruction | egress policy; no permissibly usable dataset with reference or replicate structure was obtainable, and none was written from memory |
| Where it is stated | `docs/blocked.md` B-008, `docs/validation_results.md` §4, `examples/README.md`, the `README.md` scope note, and here |
| What would change it | one permissibly licensed dataset with provenance and checksums, plus a run over it. Not a code change |

Anyone reading a `verified` cell above should read this section with it. The
compatibility of a *library version* is a real, transferable fact; the compatibility of
this package's *scientific conclusions* is a statement about synthetic fixtures only.

## Not verified here, and why

Naming these is part of the matrix. Each one is a claim someone could otherwise make by
omission.

* **The container.** `containers/Dockerfile.r-backends` has never been built, run or
  pulled. Its base image digest was resolved over the Docker Hub registry HTTP API on
  2026-09-30 — twice, same answer — and that is the only thing established about it. It is
  the digest of the multi-platform OCI *index*, not of any single platform image. Which
  repository each R package belongs to was checked against the live CRAN index and
  Bioconductor 3.20; whether the installs succeed was not.
* **CI.** Draft PR #2 has executed both Python and R suites. The first complete Python
  report has two resilience failures; R exchange tests, Ruff, Snakemake lint and both
  upstream audits passed. A corrected complete run is pending.
* **Python 3.10, 3.11 and 3.13.** Declared; the source parses under their grammars and
  nothing else about them was run.
* **Any R version other than 4.4.3.**
* **Real cytometry data.** See the section above: this is a stated scope boundary, and
  every scientific result in this package rests on synthetic fixtures.
* **The R-side test suite against every R script.** `tests/testthat/test-binmat.R` covers
  `workflow/scripts/lib/binmat.R`, the exchange library every adapter loads. The adapter
  scripts are covered by Python integration tests that drive them end to end; the R
  exchange suite is now configured in CI, but has not run on a GitHub runner.
