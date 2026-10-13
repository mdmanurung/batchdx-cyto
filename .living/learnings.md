# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### 2026-10-03 — Full package tests need the isolated R runtime and writable caches

- **Category:** test environment
- **What happened:** A sandboxed full pytest run reported backend failures when `Rscript` was absent from `PATH`; Snakemake also tried to write below the read-only home cache. With `renv/bin` on `PATH` and caches in `/tmp`, cyCombine still needed a BiocParallel localhost port that the sandbox denied. The isolated cyCombine test passed with local socket access.
- **Why it matters:** These failures can be mistaken for changes in correction algorithms or workflow behavior.
- **Resolution:** Run integration checks with the repository's `renv/bin`, writable `XDG_CACHE_HOME`, `MPLCONFIGDIR` and `NUMBA_CACHE_DIR`, and local socket access for BiocParallel.
- **Tags:** testing, R, Snakemake, sandbox
**mitigation_type**: convention
- **structural_mitigation_candidate:** A test launcher could set the local runtime and caches when this environment recurs; it cannot grant socket access inside a restricted sandbox.

### 2026-10-03 — Fresh-clone build needs a writable uv cache and fetched build backend

- **Category:** build environment
- **What happened:** The disposable fresh clone built and installed after `UV_CACHE_DIR` was directed to `/tmp` and the declared `hatchling>=1.27` backend was fetched. Offline build with an empty cache could not resolve hatchling.
- **Why it matters:** An offline failure with no cached build backend is an environment limit, not evidence that the source distribution is broken.
- **Resolution:** The fresh clone produced a wheel and sdist, installed `.[dev]` into a new virtual environment, imported version `0.1.0`, fetched all pinned upstream repositories, and passed the Python upstream audit.
- **Tags:** build, release, cache
**mitigation_type**: convention
- **structural_mitigation_candidate:** A release CI job can build from a fresh checkout with the declared build backend and inspect the resulting archives.

### 2026-10-03 — CI PATH shadowed the tested Python after R setup

- **Category:** CI environment
- **What happened:** The first draft-PR jobs exposed a GitHub API HTTP 401 when R packages were installed remotely, an R regex escape during setup verification, and an R environment `python` without pip shadowing `actions/setup-python` after the R `bin` directory was added to `PATH`.
- **Why it matters:** A local green suite does not establish that CI provisions the same interpreter and R packages; setup failures prevented the test step entirely.
- **Resolution:** Fetch pinned upstream checkouts first, install cyCombine and CytoNorm with `R CMD INSTALL`, use an escape-free R regex, and invoke the setup-python interpreter by absolute path in CI. The complete CI job remains pending verification.
- **Tags:** ci, R, python, environment
**mitigation_type**: structural
- **structural_mitigation_candidate:** Validate the GitHub Actions YAML and run the complete draft-PR job before release sign-off.

### 2026-10-03 — The public spectral-control record has complete file-to-control metadata

- **Category:** data access
- **What happened:** Zenodo record 15388817 exposes 72 FCS files and a 72-row sample table with four control labels in each of 18 batches. All filenames match and both metadata attachments match their published MD5 checksums. FlowRepository's HTTPS certificate failed validation from this host.
- **Why it matters:** The Zenodo record can seed a multi-batch technical-control validation, but its `sample_id` is a batch/control label and does not establish donor identity or biological contrasts.
- **Resolution:** Freeze the complete source manifest and the unchanged attachments; keep eligibility provisional until all FCS files and channel headers are inspected.
- **Tags:** data, provenance, spectral-flow, validation
**mitigation_type**: structural
- **structural_mitigation_candidate:** Require a design audit that maps every FCS filename and channel to the frozen source tables before this dataset supports a claim.

### 2026-10-03 — Interrupted-run hashes exposed runtime fields in scientific results

- **Category:** workflow reproducibility
- **What happened:** The complete CI suite passed 375 tests but two resilience tests found different bytes after SIGTERM and SIGKILL recovery. The changed cyCombine correction and native metric payloads included elapsed walltime; the report inherited the correction payload. The local sandbox also blocked BiocParallel ports, causing one additional local-only cyCombine failure.
- **Why it matters:** Runtime belongs in Snakemake benchmark logs; embedding it in scientific result JSON makes a correct resume look nonreproducible and changes the offline report bytes.
- **Resolution:** Remove elapsed fields from the two cyCombine R result payloads and from the correction result resources. The targeted cyCombine integration test passed with local ports enabled. The full CI resilience rerun is pending.
- **Tags:** workflow, ci, reproducibility
**mitigation_type**: structural
- **structural_mitigation_candidate:** Keep the full interrupted-run hash checks in CI and record performance separately in `benchmarks/`.

### 2026-10-03 — Public spectral controls span six assay layouts

- **Category:** data contract
- **What happened:** All 72 FCS files contained 27,022,765 events, but the headers resolved to six distinct layouts. A complete 24-file four-control grid across batches 13–18 has one exact layout.
- **Why it matters:** The apparent 18-batch grid is not a single panel; some fluorophores label different markers, so a common-channel correction across all files would be invalid.
- **Resolution:** Freeze the six-batch homogeneous subset and preserve a per-file SHA-256, layout map, and explicit panel exclusions. Biological donor identity remains unestablished.
- **Tags:** data, panel, validation
**mitigation_type**: structural
- **structural_mitigation_candidate:** Require an explicit assay-marker mapping before using other batches.

### 2026-10-03 — Inline SVG metadata made resumed reports byte-unstable

- **Category:** workflow reproducibility
- **What happened:** After elapsed fields were removed from cyCombine results, the next complete draft-PR CI run passed 377 tests but both interrupted-run tests still found a changed `report/report.html`. Repeated rendering of one Matplotlib figure showed a current `<dc:date>` and random SVG element IDs.
- **Why it matters:** A report rebuilt from unchanged analytical artifacts must have stable bytes for the resume hash check; random SVG metadata obscured that property.
- **Resolution:** Set a stable figure-specific SVG hash salt and omit the SVG date metadata. A repeatability regression test and the Stage 10 report tests passed locally. The complete CI rerun is pending because the push of commit `b313006` was rejected by automatic approval review.
- **Tags:** workflow, ci, reproducibility, report
**mitigation_type**: structural
- **structural_mitigation_candidate:** Keep the SVG repeatability test alongside the full interrupted-run CI checks.

### 2026-10-10 — Prototype defects survive passing library tests

- **Category:** package validation
- **What happened:** A fresh focused unit/scientific run passed 242 tests (four deselected), while isolated probes confirmed four defects: report.html.j2 is not autoescaped by the configured predicate; contract version 999 is accepted; a test-only policy recommends after joint/MAD guardrails are removed; and the public OpenProblems wrapper passes baseline values to both candidate EMD components.
- **Why it matters:** The first usable prototype needs tests at the public pipeline and decision boundaries. Existing library parity and uncalibrated-policy tests do not establish correct candidate dispatch or fail-closed evidence handling.
- **Resolution:** Saved the source-grounded proposed prototype checklist and receipts under review/2026-10-09/. Defects remain unfixed; shipped policies remain uncalibrated, and the dispatch probe mocks numerical component execution.
- **Tags:** prototype, pipeline, validation, report, decisions
**mitigation_type**: structural
- **structural_mitigation_candidate:** Add candidate-sensitive metric-stage tests and reachable-policy missing-evidence tests before closing prototype scoring tasks.

### 2026-10-10 — OpenProblems metric splits and control filtering change prototype data requirements

- **Category:** benchmark input contract
- **What happened:** The clean pinned upstream checkout defines six components with nine scores. The metric API requires unintegrated plus two integrated split inputs; evaluation removes controls. Input metadata describe actual donor/sample/batch/group and benchmark split assignments, distinct from train/tune/validate splits. Functional scores use a two-batch significant baseline, and FlowSOM uses lineage channels and declared grid/cluster parameters.
- **Why it matters:** A control-only dataset with no established donor identities cannot validate the complete upstream suite. Random cell halves or adapted grouping cannot silently establish benchmark-equivalent replicated samples or full component parity.
- **Resolution:** Promoted the full metric suite and two-split execution into prototype requirements and recorded a source-hashed roster receipt. Require suitable real data and component-level parity; no metric calculation was executed in this revision.
- **Tags:** prototype, metrics, openproblems, data-contract, validation
**mitigation_type**: structural
- **structural_mitigation_candidate:** Assert all nine names, input eligibility, split cell IDs, filtering counts, and full-component parity before closing the metrics gate.

### 2026-10-10 — Benchmark visual style and aggregation authority are separate

- **Category:** benchmark visualization
- **What happened:** The public cytometry leaderboard describes control-anchored normalization and an overall scaled mean. The installed and public scIB renderer instead computes a 0.4 batch/0.6 biology total and draws grouped score circles and summary bars. The local pinned task delegates reporting to an unpopulated common submodule.
- **Why it matters:** Reusing scIB's presentation can silently change the ranking rule. Missing scores, degenerate controls, and changing candidate-dependent scaling can also create misleading rankings despite attractive figures.
- **Resolution:** Required explicit anchor/aggregation provenance, complete-scope ranks, visible failure rows, normalization QC, and separate recommendation status in BENCHMARK_VISUALIZATION_SPEC.md. The current preview is synthetic and makes no scientific claim.
- **Tags:** prototype, visualization, report, openproblems, ranking
**mitigation_type**: structural
- **structural_mitigation_candidate:** Test scaling directions, missing scores, ties, fixed anchors and plot/data rank agreement at the public report boundary.

### 2026-10-10 — Validate failure isolation through both stages and the actual DAG

- **Category:** software failure containment
- **What happened:** A non-finite candidate escaped stage_correct before the fix. New regressions cover candidate execution/export errors, stale retries, failure-probe errors, cancellation and configuration boundaries. A real Snakemake execution with a simulated interpreter that omits correction output now completes with a failed candidate, empty metrics, identity output and an inspectable report.
- **Why it matters:** Passing adapter checks alone did not prove the public stage contained malformed output. The broader related suite also exposed an unchanged CytoNorm availability assertion, reproduced independently from a HEAD archive, so successful P01 checks cannot be described as a green full package gate.
- **Resolution:** Saved all pre/post, targeted, related, workflow and baseline reports with source hashes in prototype-evidence/2026-10-10-p01-verification.json. P01 is verified; strict R payload validation and real backend/scientific validation remain separate.
- **Tags:** prototype, pipeline, workflow, validation, failure-handling
**mitigation_type**: structural
- **structural_mitigation_candidate:** Retain the missing-output workflow regression alongside stage-level invalid-output/retry checks; resolve the recorded CytoNorm availability gate under the intended environment without hiding it.

### 2026-10-10 — Process success and list shape require a checked R boundary

- **Category:** cross-language result handling
- **What happened:** Six payload regressions failed before the fix; exit zero accepted a missing result or non-object/status-invalid JSON, while malformed JSON escaped the shared failure type. R auto_unbox converted single string vectors into scalars. Execution timeouts escaped before log persistence.
- **Why it matters:** A detector could report an empty successful assessment or iterate a warning's characters; subprocess failure handlers need valid structured results and partial logs, independently of backend algorithms.
- **Resolution:** Shared writer array preservation and checked envelope/list parsing now cover all six script callers. Forty focused Python tests and 110 R assertions pass, including real R-to-Python list exchange and actual interpreter-fixture timeout/failure paths. Three related CytoNorm failures reproduce on an untouched HEAD snapshot under the same default PATH.
- **Tags:** prototype, r-boundary, validation, failure-handling
**mitigation_type**: structural
- **structural_mitigation_candidate:** Keep envelope/list/timeout regressions and real R exchange tests with caller checks; do not substitute these receipts for real backend or measured-data validation.
