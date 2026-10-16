# Blocked work

Every entry records the exact command attempted, the error it produced, and which stage
it affects. Nothing below was worked around by inferring a function signature or
substituting a different algorithm, and no gate that depends on one of these is
reported as passed.

**Every remaining entry below was re-verified on 2026-09-29.** Entries whose obstacle
lifted are marked RESOLVED and kept, with the historical command and error, so the
provenance of a gate's status is never lost. A resolved entry is not a silent deletion:
`scripts/generate_status.py` reads the *executed test report*, not this file, so a status
change still requires a test that actually ran.

## The environment has changed since this file was written

The original build environment routed all outbound HTTPS through a policy-enforcing
egress proxy. PyPI, npm, crates.io and github.com were reachable; **CRAN, Bioconductor,
Posit Package Manager and the conda mirrors were not** (403 from the proxy). That single
policy was the root cause of every R-package entry here.

As of 2026-09-29 all of them answer:

```
$ for h in cloud.r-project.org bioconductor.org packagemanager.posit.co micro.mamba.pm; do
    curl -sS -o /dev/null -w "%$h %{http_code}\n" "https://$h/"; done
cloud.r-project.org        200
bioconductor.org           200
packagemanager.posit.co   307
micro.mamba.pm             307
```

The R layer is therefore built by `scripts/setup_r_env.sh` (a conda environment, not an
R library — see the toolchain notes in the roadmap plan's Decision Log), and `renv/`
carries cyCombine 0.3.0 and CytoNorm 2.0.12 at their pinned commits.

**The one obstacle this file can still not resolve is B-008: there is no permissibly
usable real cytometry dataset available here.** No gate depends on it, but every
scientific result in this package is established on deterministic synthetic fixtures,
and that is a material limit on what the numbers mean.

---

## B-001 — RESOLVED (2026-09-29): CRAN is reachable

**Stage affected:** 5, 6, 7 (native namespace), 12. **Status: closed.**

The command below is the original failure, kept for provenance:

```
$ curl -sS -o /dev/null -w "%{http_code}\n" https://cloud.r-project.org/
curl: (56) CONNECT tunnel failed, response 403
```

It now returns 200, and `src/contrib/PACKAGES` serves 258,908 lines. `kohonen`,
`emdist`, `outliers`, `dplyr`, `tidyr` and `matrixStats` are all installed in `renv/`.

**Note on *how* they are installed.** Installing R packages from CRAN **source** on this
host does not work, and the reason is a toolchain mismatch rather than the packages:
R 4.4.1 requires `GLIBCXX_3.4.29` while the default gcc 8.5.0 provides at most 3.4.25,
and CRAN's current `Rcpp` source requires `CXXABI_1.3.15` which only gcc 14 provides
(the newest gcc here is 13.2.0 → 1.3.14). Posit Package Manager has no prebuilt binary
for centos8/R 4.4, so it cannot substitute. A conda environment resolves the closure to
prebuilt packages built against one consistent toolchain, which is why
`scripts/setup_r_env.sh` uses micromamba.

---

## B-002 — RESOLVED (2026-09-29): Bioconductor is reachable

**Stage affected:** 5, 6. **Status: closed.**

Original failure, kept for provenance:

```
$ curl -sS -o /dev/null -w "%{http_code}\n" https://bioconductor.org/
curl: (56) CONNECT tunnel failed, response 403
```

Now 200. `sva` (3.54.0), `flowCore` (2.18.0) and `FlowSOM` (2.14.0) are installed in
`renv/` from `bioconductor-sva`, `bioconductor-flowcore` and `bioconductor-flowsom`.
Note that `bioconductor-cytonorm` and `r-flowsom` do not exist, so CytoNorm comes from a
GitHub install at its pinned commit.

---

## B-003 — PARTIALLY RESOLVED (2026-09-29): the mirrors are reachable; the container is still unbuilt

**Stage affected:** 11, 12. **Status: mirrors closed; container still open.**

Original failure, kept for provenance:

```
$ curl -sSL -o /dev/null -w "%{http_code}\n" https://micro.mamba.pm/api/micromamba/linux-64/latest
curl: (56) CONNECT tunnel failed, response 403
$ curl -sS -o /dev/null -w "%{http_code}\n" https://packagemanager.posit.co/cran/latest/
curl: (56) CONNECT tunnel failed, response 403
```

Both now resolve (307 / 200), and `scripts/setup_r_env.sh` uses micromamba against
conda-forge and bioconda to build `renv/` — which is how the whole R layer became
executable at all.

**What is still open, and is a tooling gap rather than a network one:** an apptainer
was later installed via conda (apptainer 1.5.4, 2026-10-16) and it can **build** a SIF
from a `docker://` source, but it **cannot execute one on this node**: every `apptainer
run` fails with a setgroups permission error (the sandbox requires user namespaces this
node does not grant). The container image in `containers/Dockerfile.r-backends`
therefore **still cannot be built from its own Dockerfile or verified by execution**,
and it remains labelled **NOT BUILT OR TESTED**. Its base image is also still a tag
(`rocker/r-ver:4.4.1`) rather than a digest pin, which is a real reproducibility gap for
a deployment artifact. The conda environment is the local substitute, not a replacement:
`renv/` is reproducible from a pinned list and isolated, but it is not a distributable
image and nothing here makes it one.

---

## B-004 — RESOLVED (2026-09-29): `cyCombine` executes and the smoke test runs

**Stage affected:** 6, 7, 5. **Status: closed.**

Original failure, kept for provenance:

```
$ Rscript -e "cat(as.character(utils::packageVersion('cyCombine')))"
  (empty; the package is not installed)
```

cyCombine 0.3.0 is now installed at the pinned commit `fcd42dbb` in `renv/`, and
`test_cycombine_smoke_test` **executes a real correction** through
`workflow/scripts/cycombine_correct.R` — asserting the output is finite, shape-stable,
preserves excluded markers, and is genuinely different from the baseline. The Stage 5
native detector comparison also runs and passes. Stage 5 is the gate that moved from
BLOCKED to PASSED.

The Stage 0 source audit still runs and is still load-bearing: `scripts/audit_upstream.R`
now reads the *installed* namespace rather than parsing source, and
`tests/integration/test_stage0_audit.py` fails the build if any of the nine
`adapter_calls` in `config/upstream_lock.yaml` is absent.

**What executing it found that the source audit could not.** These are recorded in the
code and were the substance of the old "what was done instead" list being insufficient:

* `cyCombine::detect_batch_effect` groups by a column **literally named `sample`**
  (`R/detect_batch_effect.R:200`) with no parameter to rename it. The adapter supplied
  only `batch` and `id`, so the detector aborted in `dplyr::group_by()`. A per-cell
  `sample` vector now crosses the boundary and is length-checked.
* `compute_emd` calls `emdist::emd2d()` unconditionally and `cyCombine:::check_package`
  aborts without `pbmcapply`, so three packages that are merely `Suggests` in the
  DESCRIPTION are hard requirements. All three are in `scripts/setup_r_env.sh` with the
  reason recorded.
* `evaluate_mad` returns `list(score =, mad =)`, **not** a `reduction` field. Reading the
  wrong key reported `undefined` for every input including real ones.

---

## B-005 — RESOLVED (2026-09-29): `CytoNorm` executes, and the backend is implemented

**Stage affected:** 6 (Tier B). **Status: closed.**

The blocker was availability, and it lifted with B-001/B-002. `scripts/setup_r_env.sh`
now builds an environment carrying `CytoNorm 2.0.12` at the pinned commit `278f352d`,
`FlowSOM 2.14.0` and `flowCore 2.18.0`. Note that neither `bioconductor-cytonorm` nor
`r-flowsom` exists, so the GitHub install at the pinned commit plus a conda-provided
FlowSOM is the only working combination.

The adapter is implemented, not wrapped: it writes per-sample FCS carrying event
identity in a dedicated `BATCHDX_EVENT_IDX` channel, fits on the **training split's**
measured references only, normalises a held-out split, reads the files back, and refuses
a reordered or dropped event rather than matching it. `fit_apply_verified=True` is
earned, and a test fails the build if the flag is ever True while the package is absent.

Two findings that only executing it could have produced, both now recorded in the code:

* The `transformList` must be pinned to the identity, because the baseline is already
  transformed. A declared transform list is rejected as a parameter, and the R script
  refuses anything the Python side cannot reach.
* `CytoNorm.normalize` **skips** a batch it has no fitted spline for and writes that file
  back unchanged. "A file was written" and "the cells were corrected" are different
  claims, so the result reports `n_fitted_batches`, `n_uncorrected_samples` and a
  `PARTIAL CORRECTION` warning naming the batches. On the test fixture, 12 of 18 samples
  are corrected; the injected CD4 offset moves 0.4236 mean absolute while every other
  marker moves at most 0.0112.

**What remains, and is not an availability problem:** a design-coverage refusal (no
control in the training split, a control in one batch, an arm-unbalanced split) is
currently recorded as `not_available` because no exception type maps to
`CorrectionStatus.INELIGIBLE`. `ineligible` is the honest state and is a pending change.
Also, the shipped `paired_reference` example cannot run the cytonorm candidate, because
its single reference donor cannot fill a training split; the adapter refuses correctly
and a two-reference-donor-per-batch layout is needed.

---

## B-006 — RESOLVED (2026-09-30): `CytoBatchFlagR` executes, with a specificity caveat

**Stage affected:** 5, 7 (flag-resolution row). **Status: closed as a blocker; the
remaining limit is statistical, not availability.**

The original obstacle was availability: the checkout was present but its dependencies were
not. The binding dependency set is not a DESCRIPTION — there is none — but the 29
unconditional `library()` calls at `.upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R:6-34`,
because `library()` is not wrapped in `require()` and a single missing package aborts the
whole `source()`. **All 29 are now installed**; six were added for this purpose.

`tests/integration/test_cytobatchflagr_backend.py` drives the real detector through
`run_r_script`. It needed its own fixture: the shared example has three batches, so both
earlier tests exited at the `min_batches=4` eligibility guard before reaching any of their
own assertions, and the threshold-refusal path had therefore never executed. The upstream
code fires correctly — the injected shift is flagged, and the same design unshifted is
not.

**Two findings that only execution could produce:**

1. **The 1.5-IQR fence is not specific at realistic batch counts.** Measured over 400,000
   draws it fires on four values with probability **0.207 per test** (0.0515 per value),
   and it does **not** improve with more batches (0.27 at six, 0.27 at eight, 0.30 at
   twenty). "No flags" is therefore not evidence of absence, and a reader of a report must
   not treat it as such. The clean negative result in the suite is a differential (the
   shifted batch flags, the unshifted one does not), which fence noise cannot produce;
   the module docstring records the measured rate so the clean case is not over-read.
2. **A silent eligibility constraint lives in an argument default.** `neg_mfi_df_all` and
   `pos_mfi_df_all` default to `min_cells=TRUE` and drop any control with fewer than 100
   cells for a marker, `print()`ing only a summary. `percent_pos_df_all` has no such
   argument at all. A panel too thin to evaluate would therefore have reported itself as
   *clean* — the fabricated-negative case. The wrapper now passes `min_cells=FALSE` and the
   adapter enforces `MIN_CELLS_PER_CONTROL_MARKER = 100` (upstream's number, adopted
   rather than invented), returning `ineligible` for a short panel.

Also worth recording: the README's own `package_installer()` list names `DT`,
`randomcoloR`, `shinyjs`, `shinycssloaders` and `shinyWidgets`, but those are reached only
through the interactive Shiny files the adapter deliberately never sources. They are
neither installed nor enforced, and widening the surface for code that cannot run would
cost without benefit.

**Remaining limitation:** a multiplicity-aware false-positive rate is not attached to a
`flag=True` assessment, so a report reader sees an uncalibrated flag. That is the honest
follow-up and is not forgotten.

---

## B-007 — RESOLVED (2026-09-29): extraction does NOT need the Nextflow runtime

**Stage affected:** 7 (OpenProblems rows), 8 (functional-marker component).
**Status: closed as to the open question; partially closed as to coverage.**

The open question was whether the metric code at commit `37151479` could be extracted
without the Nextflow runtime the components are built and run through. **It can.** Each
component's `helper.py` / `helper.R` is an ordinary module with plain imports; only the
`script.py` / `script.R` assembly layers read `resources_test/…` `.h5ad` files the
checkout does not carry, and those loops are transcribed and tested directly.

Seven components are now registered in the `compatibility` namespace with
`upstream_verified=True`, each proven by a test that **imports the pinned helper from the
pinned commit and compares outputs** against a recorded tolerance:

| Component | Tolerance | Agreement measured |
|---|---|---|
| `emd_mean_ct_horiz` | `rel=1e-9, abs=1e-12`; `array_equal` on `bin_array` | exact on the bin array |
| `emd_mean_ct_vert` | `rel=1e-9, abs=1e-12` | — |
| `ratio_consistent_peaks` | `atol=1e-12` on density | exact peak counts |
| `average_batch_r2_ct` | `rel=1e-9, abs=1e-12` | 2e-16 |
| `iLisi` | `rtol=0, atol=0` per cell; `rel=1e-12` on the normalisation | exact |
| `cLisi` | same, at the pinned seed 0 | exact |
| `flowsom_mean_mapping_similarity` | `abs=1e-8`, the pinned helper.R run on shared assignments | exact |

`scripts/audit_upstream.py` still fails the build if a `compatibility` row ever appears
without a recorded verification, so the claim cannot be made loosely.

One correctness finding from the comparison: the local residual-batch-R² first used a
single 0/1 contrast, which is a **different model** from the pinned
`pd.get_dummies(drop_first=True)` once a donor appears in three or more batches. The two
agree to 2e-16 only after the full dummy block is fitted, and a test now asserts the
single-contrast fit is more than 0.05 away, so the two cannot be confused again.

**Still deferred, each with a currently-true obstacle:**

* `functional_marker_preservation_wilcoxon` / `_cohens_d` — the previous deferral
  reasons were **wrong**: upstream reports a *retention ratio* and an effect-size
  *drift*, not the statistic Stage 8 computes. Both need two integrated splits, an
  uncorrected significance baseline and `var['marker_type'] == 'functional'`, which the
  one-candidate-per-frozen-partition contract does not supply. The shared per-pair
  statistics are verified against the pinned helpers.

**`cms` is still absent at the pinned commit** and its `UPSTREAM_METRICS_NOT_FOUND` row
is correct: the CHANGELOG records it being added in PR #79 and modified in PR #83, but
`src/metrics/cms` does not exist at `37151479`. A test now pins that row so it cannot be
quietly retired.

---

## B-008 — public data obtained; real-data validation remains open

**Stage affected:** 12.

All 74 files from the public CytoBatchFlagR spectral-control record were downloaded
outside Git and checksum verified on 2026-10-03. The FCS audit found six channel layouts;
the 24 files from batches 13–18 form one homogeneous four-control technical subset.
See `data/metadata/cytobatchflagr_spectral/DATASET.md` and its design audit. No raw-to-report
run or real-data correction endpoint has been executed, and no donor identity is
established. The scientific suite still runs on deterministic synthetic fixtures from
`batchdx.cyto.simulate`, which are not cohort evidence. The `0.1.0` Alpha keeps this
limit explicit; the `0.2.0` candidate needs measured real-data results.

---

## B-009 — PARTIALLY RESOLVED (2026-09-30): the machinery exists; it is not yet in the pipeline

**Stage affected:** 2 (optional preprocessing clause). **Status: partially closed.**

Two reasons this was open have both gone. Bioconductor is reachable, and the work itself
is now written: `src/batchdx/cyto/preprocess/` implements 2-D gating, the per-panel asinh
transform, the PeacoQC pass, and a `removed_events` writer, with
`tests/unit/test_preprocessing.py` (30 tests) and
`tests/integration/test_peacoqc_backend.py` (17). PeacoQC 1.16.0 executes for real: a
fixture with a known contaminated block has every injected event removed and attributed,
and both the isolation-tree and MAD paths run.

**What is still open, and it is the reason this is only partially closed: none of it is
wired into a run.** `preprocess` sits outside the pipeline by design, because it operates
on raw events *before* a baseline exists. A preprocessing stage, a Snakemake rule, and a
report section for the removal log are all still to be built, so today `init` and every
shipped example take the already-QCed bypass path. `docs/usage.md` §9 records the
extension point.

**Two schema corrections made along the way**, both removing an accepted-and-ignored
value:

* `preprocessing.transformation` offered `"logicle"`, which this build does not implement,
  has no tested inverse for, and which the pinned reference implementation does not
  contain either (`grep -i logicle` over `.upstream/spiscy` returns zero matches). It is
  removed from the `Literal`.
* `preprocessing.gating` and `preprocessing.peacoqc` were in `UNSUPPORTED_KEYS`; they are
  now real fields, read as *paths* to a gate table and a QC settings file.
  `preprocessing.time_channel` **remains refused**, because nothing makes the drift
  detector read a configured channel, and resurrecting it would recreate the very defect
  the key table exists to prevent.

**The transform convention was verified rather than assumed, and the plan's reading of it
was wrong.** `flowVS::transFlowVS` at the pinned version computes
`asinh(x / cofactor)` with no output scaling, which is exactly what `io/export.py` inverts
via `sinh(x) * cofactor` (round-trip error 4.7e-10 absolute). Adopting FlowVS is therefore
what makes the exporter's exact-inverse claim true, not a threat to it; the scaled
convention the plan recorded would have broken it by a relative error of 2.6e22.
