# Blocked work

Every entry records the exact command attempted, the error it produced, and which stage
it affects. Nothing below was worked around by inferring a function signature or
substituting a different algorithm, and no gate that depends on one of these is
reported as passed.

The environment this build was developed in routes all outbound HTTPS through a
policy-enforcing egress proxy. PyPI, npm, crates.io and github.com are reachable;
**CRAN, Bioconductor, Posit Package Manager and the conda/mamba mirrors are not.**
That single policy is the root cause of every R-package entry here.

---

## B-001 — CRAN is unreachable, so no CRAN R package can be installed

**Stage affected:** 5, 6, 7 (native namespace), 12.

```
$ curl -sS -o /dev/null -w "%{http_code}\n" https://cloud.r-project.org/
curl: (56) CONNECT tunnel failed, response 403
```

Proxy status endpoint:

```
$ curl -sS "$HTTPS_PROXY/__agentproxy/status"
  "recentRelayFailures": [
    { "kind": "connect_rejected",
      "detail": "gateway answered 403 to CONNECT (policy denial or upstream failure)",
      "host": "cloud.r-project.org:443" }, ...
```

**Consequence:** `kohonen`, `emdist`, `outliers`, `dplyr`, `tidyr`, `matrixStats` and
the rest of the CRAN dependency surface of cyCombine and CytoBatchFlagR cannot be
installed.

---

## B-002 — Bioconductor is unreachable

**Stage affected:** 5, 6.

```
$ curl -sS -o /dev/null -w "%{http_code}\n" https://bioconductor.org/
curl: (56) CONNECT tunnel failed, response 403
```

**Consequence:** `sva` (needed by cyCombine for ComBat) and `flowCore` / `FlowSOM`
(needed by CytoNorm) cannot be installed.

---

## B-003 — conda / mamba mirrors are unreachable

**Stage affected:** 11, 12.

```
$ curl -sSL -o /dev/null -w "%{http_code}\n" https://micro.mamba.pm/api/micromamba/linux-64/latest
curl: (56) CONNECT tunnel failed, response 403
$ curl -sS -o /dev/null -w "%{http_code}\n" https://packagemanager.posit.co/cran/latest/
curl: (56) CONNECT tunnel failed, response 403
```

**Consequence:** the per-rule conda environments described in `workflow/envs/` cannot be
created or solved here, and the container image in `containers/Dockerfile.r-backends`
cannot be built (its `BiocManager::install` and `install.packages` layers need CRAN and
Bioconductor). The Dockerfile is shipped labelled **NOT BUILT OR TESTED**.

---

## B-004 — `cyCombine` cannot be executed, so Stage 6 `cycombine` has no smoke test

**Stage affected:** 6 (gate clause "all implemented backends execute real smoke tests"),
7 (native EMD and MAD), 5 (both cyCombine detectors).

```
$ Rscript -e "cat(as.character(utils::packageVersion('cyCombine')))"
  (empty; the package is not installed)
```

Installing it requires B-001 and B-002.

**What was done instead:**

* The adapter and its R script are written against the pinned API and the **source at
  the pinned commit was audited**, not guessed: `scripts/audit_upstream.R` parses
  `fcd42dbb` and prints the real formals of every function the adapter calls
  (`normalize`, `create_som`, `correct_data`, `compute_emd`,
  `detect_batch_effect[_express]`, `evaluate_emd`, `evaluate_mad`, `get_markers`).
  `tests/integration/test_stage0_audit.py` fails the build if any of them is absent.
* The adapter reports `not_available` with an actionable reason and never substitutes
  another algorithm.
* `tests/integration/test_stage6_backends.py::test_cycombine_smoke_test` is marked
  `requires_cycombine` and **skips visibly**. The Stage 6 gate is recorded as
  **BLOCKED**, not passed.

---

## B-005 — `CytoNorm` cannot be executed, and the backend is deliberately not implemented

**Stage affected:** 6 (Tier B).

Same root cause as B-004, plus `flowCore` and `FlowSOM`.

Beyond availability, the audited API shapes the work: `CytoNorm.train(files, labels,
channels, transformList, ...)` and `CytoNorm.normalize(model, files, ...)` operate on
**FCS files**, so a working adapter needs a tested FCS round trip and a `transformList`
consistent with the declared processing manifest. That is real work, not a wrapper, and
it is not attempted here. The backend registers itself and raises with an explicit
deferred reason; a configuration naming it gets a status, not a fabricated correction.

---

## B-006 — `CytoBatchFlagR` cannot be executed

**Stage affected:** 5, 7 (flag-resolution row).

The checkout is available (`scripts/fetch_upstream.sh` clones commit `ad21e3b1` from
GitHub, which is reachable), and the 1.5 IQR rule and the three IQR metrics were read
directly from `R/CytoBatchFlagR_functions.R`. What is missing is its ~35 CRAN
dependencies (B-001).

The adapter is written and reports `not_available` when either the checkout or the
dependencies are absent. Its eligibility constraints — a control sample in every batch,
and at least four batches — are encoded in Stage 3 and are separately proven:
`tests/scientific/test_stage3_design.py::test_iqr_fence_cannot_flag_with_three_values`
shows arithmetically, over 2000 random draws, that the Tukey fence can never fire for
three observations.

---

## B-007 — OpenProblems metric implementations were not extracted

**Stage affected:** 7 (all OpenProblems rows), 8 (functional-marker component).

The metric code at commit `37151479` lives inside Viash components
(`src/metrics/*/config.vsh.yaml` plus `script.py` / `script.R`) that are built and run
through Nextflow with containerised deployment. Whether those implementations can be
extracted without the Nextflow runtime was listed as unresolved in the plan; **it was
not attempted in this build**, so no metric here claims to be an upstream
compatibility component, and `scripts/audit_upstream.py` fails the build if one ever
does.

What *was* done: the metric inventory was reconciled against the pinned commit. Nine
metric names exist across six components. **The `cms` (CellMixS) metric named in the
plan is not present**: the CHANGELOG records it being added in PR #79 and later
modified in PR #83, but `src/metrics/cms` does not exist at `37151479`. It is recorded
as unavailable in `batchdx.cyto.metrics.registry.UPSTREAM_METRICS_NOT_FOUND` and shown
in the report.

---

## B-008 — no permissibly usable real cytometry example is included

**Stage affected:** 12.

No real dataset with reference or replicate structure was obtainable under the egress
policy above, and none was included from memory. The scientific suite runs entirely on
deterministic synthetic fixtures from `batchdx.cyto.simulate`, which are labelled as
fixtures and are not a scientific claim.

---

## B-009 — PeacoQC-based preprocessing is not implemented

**Stage affected:** 2 (optional preprocessing clause).

PeacoQC is an R/Bioconductor package (B-002). The configuration schema declares the
preprocessing block, refuses an `asinh` transformation without a per-panel cofactor, and
refuses to transform data that the processing manifest says is already transformed.
Gating, PeacoQC and per-removed-event QC tracking are **not implemented**; Tier A uses
the already-QCed bypass path.
