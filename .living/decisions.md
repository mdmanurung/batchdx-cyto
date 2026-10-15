# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### 2026-10-03 — Add Mycelium without moving application files

- **Context:** This is an existing Python/R package with imports, packaged workflow assets, tests, and substantial uncommitted work.
- **Decision:** Initialize the Mycelium living layer and manifests in place. Keep `src/`, `workflow/`, `docs/`, `config/`, and existing project files at their current paths.
- **Rationale:** Moving those files would break package paths and workflow references without improving the requested initialization.
- **Consequences:** The Mycelium manifests begin as registries for future analyses and data; the existing package layout remains authoritative.
- **Tags:** mycelium, repository-layout

### 2026-10-03 — Use Python 3.12 as the Alpha CI status environment

- **Context:** The workflow environment pins Python 3.12, while the prior CI job selected 3.11 and the generated status came from a local 3.12 run.
- **Decision:** Select Python 3.12 in the primary CI job and run the R exchange suite there with the job's Python interpreter.
- **Rationale:** The status comparison needs an executed environment matching the local report, and cross-language tests must not silently skip.
- **Consequences:** Python 3.10, 3.11 and 3.13 remain declared-only until an executed matrix establishes support; the CI job remains unverified until it runs.
- **Tags:** ci, release, python

### 2026-10-03 — Freeze Zenodo spectral-control metadata before raw-data validation

- **Context:** FlowRepository access failed TLS validation, while the public Zenodo record offered a complete, checksum-bearing fallback with 72 FCS files across 18 batches.
- **Decision:** Commit the complete 74-file source manifest and unchanged small metadata attachments, plus an all-file downloader; keep raw FCS outside Git and defer eligibility until the full data and channel headers are audited.
- **Rationale:** The source metadata establishes file-to-control coverage but not donor identity, marker-channel correspondence, or endpoint validity.
- **Consequences:** This is B1 preparation only. No real-data correction or biological-preservation claim follows from it.
- **Tags:** data, provenance, validation

### 2026-10-03 — Restrict the spectral-control validation to an exact-layout subset

- **Context:** All 72 public FCS files passed source MD5 and local SHA-256 auditing, but their batches use six distinct channel layouts. Batches 13–18 alone contain all four measured controls under one exact 33-channel layout.
- **Decision:** Freeze those 24 files as the first technical-control subset and record the remaining 48 as panel exclusions. Do not pool the full 18-batch collection until marker harmonization is specified and validated.
- **Rationale:** Shared fluorophore names map to different markers in some panels; pooling by channel name would silently mix biological measurements.
- **Consequences:** The subset can enter technical-control workflow review after a marker manifest and raw preprocessing exist. It cannot support donor-disjoint or biological-preservation claims.
- **Tags:** data, panel, validation

### 2026-10-03 — Make inline report SVG stable per figure

- **Context:** Complete CI recovery tests found only `report/report.html` changed after resume. Matplotlib SVG output included a creation date and random element IDs even for identical plot inputs.
- **Decision:** Omit SVG date metadata and use a stable, figure-specific hash salt for element IDs.
- **Rationale:** Repeated reports need identical bytes, while figures embedded together must retain separate element IDs.
- **Consequences:** Report rendering is repeatable under the tested Matplotlib version; the complete interrupted-run CI gate remains pending until the local fix is pushed and rerun.
- **Tags:** report, workflow, reproducibility

### 2026-10-10 — Require the full pinned OpenProblems metric suite for the prototype

- **Context:** The first-prototype checklist originally allowed deferring FlowSOM and the two upstream functional-marker metrics and using custom scores as the core evaluation path.
- **Decision:** The user instructed that OpenProblems cytometry batch integration supply the full measured metric set. Require all nine score outputs across the six components at pinned commit 37151479aaa89edf37ae23068de8fbb41a1b7bb4. Custom common/native metrics remain supplemental diagnostics, not substitutes.
- **Rationale:** A fixed upstream suite defines the coverage and comparability target; six existing local compatibility implementations do not establish full component/workflow parity.
- **Consequences:** P06–P09 must include all nine scores; P13 must execute the two benchmark correction splits; P16–P17 need a genuine full-suite eligible measured dataset. Controls are excluded by upstream evaluation, so the current control-only spectral subset cannot alone close the gate. Plan estimate expands to 12–18 engineer-days plus preparation/access delays. No implementation or release was performed or authorized by this documentation update.
- **Tags:** prototype, metrics, openproblems, validation

### 2026-10-10 — Require OpenProblems/scIB-style benchmark rank visualizations

- **Context:** After requiring all nine OpenProblems cytometry scores, the user also requested the OpenProblems/scIB benchmark ranking visualization style.
- **Decision:** Make the grouped method-by-metric rank table, category/overall bars, normalization QC, and offline/vector/data exports part of P08/P20. Use a declared complete-scope mean of control-scaled scores for overall ranking, with no implicit scIB 40/60 weights. Record scaling controls and keep partial/failing candidates visible and unranked overall.
- **Rationale:** The cytometry website specifies ranking by scaled overall mean; scIB provides the desired grouped circle/bar table layout. The visual layout does not itself establish numerical benchmark parity or calibrated correction recommendations.
- **Consequences:** Added BENCHMARK_VISUALIZATION_SPEC.md and a clearly labelled synthetic layout preview with a reproducing review script. Reference reporting/scaling code must be pinned and verified before parity claims; the local task's common rendering submodule is unpopulated. No package renderer, measured corrections, or recommendation calibration was implemented.
- **Tags:** prototype, report, visualization, openproblems, ranking

### 2026-10-10 — Contain expected candidate failures at the shared correction stage

- **Context:** The user requested “next” after the prototype plan identified P01 as the first bounded implementation ticket. A regression reproduced non-finite output escaping the stage and stopping other candidates.
- **Decision:** Contain backend/value/file/timeout errors around candidate execution and export, retaining failed/not_available records and tracebacks. Invalidate only that candidate's previous result/H5AD/metrics before retry, and validate its directory ID and backend early. Keep shared-input loading, cancellation, unexpected defects and inability to persist evidence fatal.
- **Rationale:** The same stage serves CLI and Snakemake. Expected candidate failures should remain inspectable without producing baseline fallbacks or scores from previous successful attempts. A software probe must not mask the original failure.
- **Consequences:** P01 acceptance checks pass, including a simulated missing-output Snakemake run. A broader related test run has one unchanged CytoNorm availability assertion failure independently reproduced on HEAD; full validation remains unapproved by evidence. P02 payload/list/timeout checks, nine-score integration and rank plots are still pending.
- **Tags:** prototype, pipeline, failure-handling, validation

### 2026-10-10 — Validate R envelopes and bound subprocess execution

- **Context:** P02 regressions reproduced missing/malformed payloads being accepted or escaping the detector failure handler despite exit zero, and the shared R writer could unbox single warnings/node IDs.
- **Decision:** Require the existing JSON envelope and declared string arrays, preserve list cardinality in the shared writer, and use RScriptFailed for malformed results and execution timeouts. Keep partial logs and report the actual reason in all six callers. Default execution to 14,400 seconds, matching the workflow correction budget; bound version probes to 30 seconds. Keep domain-output validation with callers.
- **Rationale:** One shared change reaches correction, detection, preprocessing and native metrics without a new RPC protocol or dependency. Timed-out processes have no observed exit code, so record None rather than a synthetic status.
- **Consequences:** Scoped Python and R exchange checks pass. Broader CytoNorm failures under the default PATH reproduce on HEAD and remain open; real backend/scientific validation and the full Python gate are separate. P03 freshness, full OpenProblems scoring and rank visuals remain pending.
- **Tags:** prototype, r-boundary, failure-handling, validation
