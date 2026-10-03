# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### 2026-10-03 — Full package tests need the isolated R runtime and writable caches

- **Category:** test environment
- **What happened:** A sandboxed full pytest run reported backend failures when `Rscript` was absent from `PATH`; Snakemake also tried to write below the read-only home cache. With `renv/bin` on `PATH` and caches in `/tmp`, cyCombine still needed a BiocParallel localhost port that the sandbox denied. The isolated cyCombine test passed with local socket access.
- **Why it matters:** These failures can be mistaken for changes in correction algorithms or workflow behavior.
- **Resolution:** Run integration checks with the repository's `renv/bin`, writable `XDG_CACHE_HOME`, `MPLCONFIGDIR` and `NUMBA_CACHE_DIR`, and local socket access for BiocParallel.
- **Tags:** testing, R, Snakemake, sandbox
- **mitigation_type:** convention
- **structural_mitigation_candidate:** A test launcher could set the local runtime and caches when this environment recurs; it cannot grant socket access inside a restricted sandbox.

### 2026-10-03 — Fresh-clone build needs a writable uv cache and fetched build backend

- **Category:** build environment
- **What happened:** The disposable fresh clone built and installed after `UV_CACHE_DIR` was directed to `/tmp` and the declared `hatchling>=1.27` backend was fetched. Offline build with an empty cache could not resolve hatchling.
- **Why it matters:** An offline failure with no cached build backend is an environment limit, not evidence that the source distribution is broken.
- **Resolution:** The fresh clone produced a wheel and sdist, installed `.[dev]` into a new virtual environment, imported version `0.1.0`, fetched all pinned upstream repositories, and passed the Python upstream audit.
- **Tags:** build, release, cache
- **mitigation_type:** convention
- **structural_mitigation_candidate:** A release CI job can build from a fresh checkout with the declared build backend and inspect the resulting archives.

### 2026-10-03 — CI PATH shadowed the tested Python after R setup

- **Category:** CI environment
- **What happened:** The first draft-PR jobs exposed a GitHub API HTTP 401 when R packages were installed remotely, an R regex escape during setup verification, and an R environment `python` without pip shadowing `actions/setup-python` after the R `bin` directory was added to `PATH`.
- **Why it matters:** A local green suite does not establish that CI provisions the same interpreter and R packages; setup failures prevented the test step entirely.
- **Resolution:** Fetch pinned upstream checkouts first, install cyCombine and CytoNorm with `R CMD INSTALL`, use an escape-free R regex, and invoke the setup-python interpreter by absolute path in CI. The complete CI job remains pending verification.
- **Tags:** ci, R, python, environment
- **mitigation_type:** structural
- **structural_mitigation_candidate:** Validate the GitHub Actions YAML and run the complete draft-PR job before release sign-off.

### 2026-10-03 — The public spectral-control record has complete file-to-control metadata

- **Category:** data access
- **What happened:** Zenodo record 15388817 exposes 72 FCS files and a 72-row sample table with four control labels in each of 18 batches. All filenames match and both metadata attachments match their published MD5 checksums. FlowRepository's HTTPS certificate failed validation from this host.
- **Why it matters:** The Zenodo record can seed a multi-batch technical-control validation, but its `sample_id` is a batch/control label and does not establish donor identity or biological contrasts.
- **Resolution:** Freeze the complete source manifest and the unchanged attachments; keep eligibility provisional until all FCS files and channel headers are inspected.
- **Tags:** data, provenance, spectral-flow, validation
- **mitigation_type:** structural
- **structural_mitigation_candidate:** Require a design audit that maps every FCS filename and channel to the frozen source tables before this dataset supports a claim.

### 2026-10-03 — Interrupted-run hashes exposed runtime fields in scientific results

- **Category:** workflow reproducibility
- **What happened:** The complete CI suite passed 375 tests but two resilience tests found different bytes after SIGTERM and SIGKILL recovery. The changed cyCombine correction and native metric payloads included elapsed walltime; the report inherited the correction payload. The local sandbox also blocked BiocParallel ports, causing one additional local-only cyCombine failure.
- **Why it matters:** Runtime belongs in Snakemake benchmark logs; embedding it in scientific result JSON makes a correct resume look nonreproducible and changes the offline report bytes.
- **Resolution:** Remove elapsed fields from the two cyCombine R result payloads and from the correction result resources. The targeted cyCombine integration test passed with local ports enabled. The full CI resilience rerun is pending.
- **Tags:** workflow, ci, reproducibility
- **mitigation_type:** structural
- **structural_mitigation_candidate:** Keep the full interrupted-run hash checks in CI and record performance separately in `benchmarks/`.

### 2026-10-03 — Public spectral controls span six assay layouts

- **Category:** data contract
- **What happened:** All 72 FCS files contained 27,022,765 events, but the headers resolved to six distinct layouts. A complete 24-file four-control grid across batches 13–18 has one exact layout.
- **Why it matters:** The apparent 18-batch grid is not a single panel; some fluorophores label different markers, so a common-channel correction across all files would be invalid.
- **Resolution:** Freeze the six-batch homogeneous subset and preserve a per-file SHA-256, layout map, and explicit panel exclusions. Biological donor identity remains unestablished.
- **Tags:** data, panel, validation
- **mitigation_type:** structural
- **structural_mitigation_candidate:** Require an explicit assay-marker mapping before using other batches.
