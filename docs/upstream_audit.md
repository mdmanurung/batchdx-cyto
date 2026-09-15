# Upstream audit (Stage 0)

This file **interprets** the audit artifacts; it does not replace them. The artifacts
are produced by two executable scripts and are the authority:

```
scripts/fetch_upstream.sh .upstream
python  scripts/audit_upstream.py --upstream-dir .upstream --out docs/audit/audit_python.json
Rscript scripts/audit_upstream.R docs/audit/audit_r.json .upstream config/upstream_lock.yaml
```

Both scripts exit non-zero when something does not reconcile: a checkout at the wrong
commit, a metric this package claims to reconcile against that is absent upstream, a
local metric claiming to be an upstream compatibility component, or an
adapter-referenced R function that does not exist in the audited source.
`config/upstream_lock.yaml` is the pinned manifest they read.

The R audit runs in one of two modes and records which:

* **installed** — the package is on the library path; formals come from the loaded
  function objects.
* **source** — the package is not installed; the checkout at the pinned commit is
  *parsed* (never evaluated) and formals come from the assignment expressions.

Neither mode infers a signature. Both read one. In the environment this build was
developed in, all three R resources were audited in **source** mode, because CRAN and
Bioconductor are unreachable (`docs/blocked.md`).

---

## cyCombine — `biosurf/cyCombine` @ `fcd42dbbe92bfa6bffa6094bf9b4de536e71c249`

Version 0.3.0, MIT (the short `LICENSE` stub points at `LICENSE.md`, which carries the
full text; the audit scans both).

Confirmed against the source, correcting or sharpening the plan's premises:

* **Exports.** All of `detect_batch_effect`, `detect_batch_effect_express`,
  `evaluate_emd`, `evaluate_mad`, `prepare_data`, `create_som`, `correct_data`,
  `normalize`, `transform_asinh`, `impute_across_panels` and `salvage_problematic` are
  exported, as the plan said. The plan did not mention `compute_emd`, `compute_mad`,
  `compute_lisi` and `compute_cilisi`, which are also exported and which matter: see
  below.
* **`ComBat_seq` is reachable.** `correct_data(df, label, markers, method = c("ComBat",
  "ComBat_seq"), ...)`, and the NAMESPACE imports both `sva::ComBat` and
  `sva::ComBat_seq`. The adapter therefore pins `method = "ComBat"` and *rejects* the
  parameter outright — a user cannot pass either value, because the route is not theirs
  to choose. `workflow/scripts/cycombine_correct.R` refuses to run if the caller sends
  anything else.
* **The 3-to-30-batch rule is real and is in the code, not only the docs.**
  `R/detect_batch_effect.R` line 336: `if (length(levels(df$batch)) >= 3 &
  length(levels(df$batch)) <= 30)`. Outside that range the function emits UMAPs only.
  Encoded as `ELIGIBILITY_RULES["cycombine_cluster"]` with `min_batches=3`,
  `max_batches=30` and the origin string quoting the source.
* **Neither detector returns anything structured.** Both functions communicate through
  `message()` calls and plot files. A structured table can therefore only come from
  calling the exported `compute_emd` directly, which is exactly what
  `workflow/scripts/cycombine_detect.R` does. Nothing scrapes console output and
  nothing reads a rendered image.
* **The EMD is binned and lives in bin units.** `compute_emd` builds bins as
  `seq(floor(min(df[,markers])) - binSize, ceiling(max(...)) + binSize, by = binSize)`
  using *global* limits across all markers, fills them with `hist(...)$density`, and
  calls `emdist::emd2d` on single-column matrices. That is a different quantity from an
  empirical Wasserstein distance in data units, which is why this package's own metric
  is named `common.emd_within_population` and says so in its registry notes.
* **The MAD score is unsigned.** `evaluate_mad` returns
  `median(|MAD_uncorrected - MAD_corrected|)`, so it cannot distinguish a correction
  that restores spread from one that destroys it. That is why
  `common.mad_signed_change` exists alongside it under a different name.
* **`evaluate_emd` recomputes both sides.** It computes EMD for the corrected *and* the
  uncorrected data, each with its own `label` column, then filters pairs below
  `filter_limit = 2` in both states. Under the benchmarking procedure where each
  correction re-clusters, the "uncorrected" baseline therefore moves. That behaviour is
  preserved in native compatibility mode only and is never mixed with common-partition
  numbers (`docs/decisions/0005-evaluation-partitions.md`).
* **Fit and apply are not separable.** `correct_data` returns a corrected table, not a
  reusable model, so the adapter declares `supports_fit_apply=False` and does not
  advertise frozen application to unseen batches.

---

## CytoNorm — `saeyslab/CytoNorm` @ `278f352d8edffd7a5d3e6b0b76c402d991a1e7a9`

Version 2.0.12, **GPL (>= 2)** (read from `DESCRIPTION`).

* **Fit and apply *are* separable.** `CytoNorm.train(files, labels, channels,
  transformList, outputDir, FlowSOM.params, normMethod.train, normParams, seed, ...)`
  and `CytoNorm.normalize(model, files, labels, transformList, transformList.reverse,
  outputDir, prefix, ...)` are distinct, and `getCytoNormQuantiles(model, type)` exposes
  the fitted quantiles. This resolves one of the plan's open questions in the
  affirmative.
* **They operate on FCS files, not data frames.** That is the reason the backend is
  deferred rather than wrapped: a correct adapter has to write per-sample FCS from the
  immutable baseline, hand CytoNorm a `transformList` consistent with the declared
  processing manifest, and read the normalised files back. None of that is implemented
  or tested here.
* **Behaviour under missing reference populations** — the plan's other open question —
  **was not resolved**, because it cannot be determined from signatures alone and the
  package could not be executed. It stays open.
* Licence consequence: GPL-2-or-later means no source may be copied into this MIT
  package. Arm's-length invocation only.

---

## CytoBatchFlagR — `AndorfLab/CytoBatchFlagR` @ `ad21e3b104126783044351eb11fa24622a0da2f7`

* **The rename is confirmed by GitHub itself.** Cloning `AndorfLab/cytoFlagR` and
  `AndorfLab/CytoBatchFlagR` yields the same HEAD commit, i.e. the old name redirects.
  Both names are recorded in the lock (`repository` and `former_repository`). Cite the
  published name; record the earlier `cytoFlagR` identifier.
* **It is not an installable package.** The checkout has `R/` and `LICENSE` and nothing
  else: no `DESCRIPTION`, no `NAMESPACE`. ~7900 lines of R across nine files.
* **It is GPL-3.0**, which is the decisive fact for integration. See
  `docs/decisions/0003-cytobatchflagr-integration.md`: the scripts are `source()`d
  unmodified from a user-supplied checkout, in a separate R process, and nothing is
  copied or refactored into this package.
* **Its entry point is interactive.** `R/CytoBatchFlagR_workflow.R` calls
  `rstudioapi::selectDirectory`; `R/editMarkerTable.R`, `R/selectMarkerList.R` and
  `R/referenceMarkerSelection.R` are Shiny apps. The adapter never sources those; it
  calls the analysis functions in `R/CytoBatchFlagR_functions.R` with configured values.
* **The IQR rule, read from the source:** `detect_outlier(x)` is
  `x > Q3 + 1.5*IQR(x) | x < Q1 - 1.5*IQR(x)`, applied per marker, grouped by control,
  across batches, for three metric types (`-MFI`, `+MFI`, `%pos`). Two constraints
  follow and are encoded in Stage 3:
  1. a control sample must exist in **every** batch;
  2. for n ≤ 3 batches the fence can never be exceeded, so the assessment is degenerate
     below four batches. This is arithmetic, not a guess, and is proven over random
     draws in `test_iqr_fence_cannot_flag_with_three_values`.
* **Thresholds are the analyst's, not the tool's.** The publication reports five
  experienced analysts setting different positive/negative thresholds for the same
  markers. The adapter refuses to choose them and requires
  `evaluation.gating_thresholds` in the declared coordinate system.

---

## OpenProblems cytometry batch integration — `openproblems-bio/task_cyto_batch_integration` @ `37151479aaa89edf37ae23068de8fbb41a1b7bb4`

MIT. Viash components executed through Nextflow.

**Nine metric names exist at the pinned commit, across six components** — reconciled
from the actual `src/metrics/*/config.vsh.yaml` files, not from a summary:

| Component | Metric names |
|---|---|
| `average_batch_r2` | `average_batch_r2_ct` |
| `emd` | `emd_mean_ct_horiz`, `emd_mean_ct_vert` |
| `flowsom_mapping_similarity` | `flowsom_mean_mapping_similarity` |
| `functional_marker_preservation` | `functional_marker_preservation_wilcoxon`, `functional_marker_preservation_cohens_d` |
| `lisi` | `iLisi`, `cLisi` |
| `ratio_consistent_peaks` | `ratio_consistent_peaks` |

**Correction to the plan's premise:** the plan states that "a CellMixS-based `cms`
metric is present in the repository". **It is not present at the pinned commit.** The
CHANGELOG records `Added metrics/cms (PR #79)` and a later entry about disabling its
parallelisation (PR #83), but `src/metrics/cms` does not exist in the tree. It is
recorded as unavailable rather than as something to reconcile against, and the audit
asserts it stays that way.

Every one of these metrics is **deferred** in this build (Tier B), and B-007 in
`docs/blocked.md` records that extraction from the Viash components was not attempted.

---

## SPiSCy — `emilier17/spiscy` @ `a6e55c8f722112bf01d5364996b1867844efbd02`

MIT (verified from the checkout, which the plan listed as unverified). A Snakemake
workflow with `workflow/{Snakefile,rules,scripts,profiles}`, a `docker/` directory and
shell entry points. No release and no associated publication was located.

Treated as reference reading for workflow layout only. **No Tier A capability depends on
it**, and nothing here is built against it as a contract.

---

## Environments where versions conflict

None were needed: the Python side has no conflict, and the R side is a single
environment per the container definition. Should cyCombine and CytoNorm ever require
incompatible versions of a shared dependency, the subprocess boundary
(`docs/decisions/0002-python-r-boundary.md`) already makes separate environments a
configuration change (`upstream.rscript`) rather than a redesign.
