# Remaining work to a first usable batchdx-cyto prototype

Reviewed 2026-10-10. Branch `alpha-ci-20261003`, HEAD
`b313006ea175c53c7ce709f7b894d1768b7db0e2`, with pre-existing uncommitted
documentation, data-metadata, gate, and review changes. This is a proposed task
list, not implementation authorization or a release sign-off.

**Scope update, user instruction 2026-10-10:** use the full OpenProblems cytometry
batch-integration metric suite as the primary measured evaluation set. All nine
scores at the existing upstream pin are prototype requirements. FlowSOM and
functional-marker metrics are no longer deferred extensions. Custom common/native
scores are supplemental diagnostics, not replacements for this suite.

**Start with P01–P05. The first milestone is safe execution and trustworthy
artifacts. Adding correction methods is not the immediate bottleneck.**

## 1. Proposed prototype boundary

The first usable prototype should let a researcher supply one already-QCed,
explicitly transformed, single-panel dataset; validate its sample/marker/processing
manifests; diagnose batch structure; compare identity with cyCombine; save a
corrected H5AD; inspect a self-contained report; and rerun or recover without
silently using stale data. It should run locally in one documented environment,
with all nine OpenProblems scores implemented and measured on eligible inputs.

The prototype offers exploratory research evidence. It does not automatically
approve corrections or claim biological preservation from technical controls.
Unsupported inputs, absent measurements, partial correction, failed metrics, and
ineligible designs must remain visible.

The fastest route uses **preprocessed H5AD**. A documented, reproducible preparation
script may create that baseline from the already-audited spectral-control subset;
this does not require turning the whole raw-FCS preprocessing library into a
product workflow first. Record every preparation operation and retain source
event identity. Do not call supplied FCS “already transformed” merely because it
has been gated. Do not silently change the frozen dataset selection or marker
exclusions.

The spectral-control subset remains useful for supplemental technical diagnostics,
but it cannot be the sole full-suite validation dataset: the upstream evaluation
removes controls and requires actual donor/sample/biological-group and benchmark
split metadata. Add a suitable benchmark/cohort dataset for the full measured suite.

CytoNorm and automated raw-FCS ingestion are extensions with their own acceptance
checks below. If either is included in the prototype's advertised supported path,
its extension tasks become prototype requirements. Diagnostic-only local use
could be delivered earlier, but it would not satisfy the correction prototype.

| Milestone | Tasks | Completion outcome | Planning estimate |
|---|---|---|---|
| M1: safe execution | P01–P05 | Bad candidates do not invalidate other results; reports and outputs are safe and versioned | 2–3 engineer-days |
| M2: full OpenProblems scoring | P06–P10 | All nine upstream scores are implemented with correct split/input semantics and parity evidence | 6–9 engineer-days |
| M3: correction and diagnostic semantics | P11–P15 | The configured biology, marker masks, fit scope, coverage, and diagnostic limits are honored | 2–3 engineer-days |
| M4: measured-data usability | P16–P20 | A full-suite eligible measured-data demonstration and an installable, documented prototype | 3–5 engineer-days |

Allow **13–20 engineer-days** for the full metric suite plus requested rank visuals, excluding data preparation/access
delays and CI infrastructure repair. These are planning estimates, not measured
delivery times. Tasks overlap; do not sum every old reviewer estimate. A usable
preprocessed dataset available now substantially reduces M4 uncertainty.

### Required metric roster and benchmark contract

Authority: `config/upstream_lock.yaml`, OpenProblems commit
`37151479aaa89edf37ae23068de8fbb41a1b7bb4`. The clean local checkout was verified
against that pin on 2026-10-10. Six component configurations declare nine scores:

| Component | Required score | Upstream score direction |
|---|---|---|
| emd | `emd_mean_ct_horiz` | lower |
| emd | `emd_mean_ct_vert` | lower |
| average_batch_r2 | `average_batch_r2_ct` | lower |
| lisi | `iLisi` | higher |
| lisi | `cLisi` | higher for the upstream normalized score |
| ratio_consistent_peaks | `ratio_consistent_peaks` | higher |
| flowsom_mapping_similarity | `flowsom_mean_mapping_similarity` | higher |
| functional_marker_preservation | `functional_marker_preservation_wilcoxon` | higher |
| functional_marker_preservation | `functional_marker_preservation_cohens_d` | lower |

The upstream metric API reads unintegrated data plus **two integrated benchmark
splits**. Its metadata include cell type, sample, donor, biological group, batch,
control status, split, marker type and correction mask; FlowSOM additionally uses
declared grid dimensions and cluster count. Benchmark splits are not the package's
train/tune/validate splits. Honor the upstream sample/donor/batch pairing and
control/unlabelled-cell filtering; merely cutting each candidate's cells in half
is not proof of upstream-equivalent splitting or independent technical replication.

Functional-marker preservation uses baseline-significant marker/cell-type pairs
in the two-batch upstream design, consistent biological-group orientation, and
both integrated splits. Implement its exact retention/effect-drift definitions,
including its declared raw-p-value rule; any donor-bootstrap/FDR extension must
be separately named supplemental evidence. Do not silently select the first two
batches of a multi-batch cohort. Freeze an eligible two-batch scope or explicitly
record a separately validated adaptation.

Every configured candidate, including identity, must have all nine named records
in the report/JSON: either a valid score with its measured scope, or an explicit
ineligible/undefined/failed reason. Missing metadata never becomes a zero score
or an omitted row. **Implementation completeness requires all nine to execute on
eligible validation inputs; a dataset where everything is ineligible does not
close the prototype gate.** The measured full-suite demonstration needs suitable
real inputs, including a defined functional baseline. `cms` is absent at this pin
and is not added to this roster without an explicit upstream-version change.

Verified roster receipt:
[`prototype-evidence/2026-10-10-openproblems-metric-contract.json`](prototype-evidence/2026-10-10-openproblems-metric-contract.json).

**User visualization requirement (2026-10-10):** include an OpenProblems/scIB-style
ranked benchmark table for all nine scores, grouped metric glyphs and summary bars,
identity and failed-candidate rows, normalization QC, and offline/exportable values.
Use documented control-scaled overall means; do not silently import scIB's 40/60
weights or turn descriptive rank into calibrated approval. The complete contract,
layout preview, ranking rules, and checks are in
[BENCHMARK_VISUALIZATION_SPEC.md](BENCHMARK_VISUALIZATION_SPEC.md).

## 2. What already exists

Reuse the current implementation rather than starting again:

1. Python CLI and a Snakemake DAG with diagnose, benchmark, selected correction,
   and report modes; strict configuration; sample, marker, and processing records;
   H5AD/FCS I/O; baseline hashing; separate corrected H5AD output.
2. Identity, cyCombine, and CytoNorm adapters; isolated R subprocesses; pinned
   upstream source manifests; event-order and excluded-marker checks. Availability
   on a new machine still needs an executed environment check.
3. Frozen labels and evaluation sampling, split machinery, common metrics,
   decision states, and an offline HTML report with metric tables and plots.
4. Library implementations for gating, asinh, PeacoQC, OpenProblems compatibility
   metrics, and sample/donor-level contrast statistics. **These are not all called
   by the normal workflow. Library tests do not establish that a user run uses them.**
5. Extensive synthetic tests and a checksummed spectral-data inventory. The saved
   `build/report.json` reports 355 passed, 1 skipped, and 23 deselected; it is
   historical non-slow evidence, not a complete current validation run. The saved
   dataset audit is preparation evidence, not a correction result.

## 3. Prototype tasks, in execution order

P01 and P02 are implemented and verified as of 2026-10-10; remaining items are outstanding
unless implementation and acceptance evidence say otherwise. The parenthesized L# identifiers refer to the existing review lenses;
they identify provenance, not independent confirmation of every reviewer claim.

### M1 — Safe execution and artifacts

1. **P01 — Contain candidate failures and retain their evidence.**
   **Status: implemented; scoped acceptance checks passed.**
   [P01 implementation and limits](P01_FAILURE_CONTAINMENT.md) and
   [receipt](prototype-evidence/2026-10-10-p01-verification.json).
   The reviewed stage let output-validation errors, missing files, malformed JSON,
   timeouts and output-writing errors escape. It now normalizes expected candidate
   execution/export errors into explicit records and logs, validates backend names
   and directory IDs early, and invalidates old outputs/results/metrics before a
   retry. Cancellation, shared-input failures and unexpected defects still raise.
   The broader suite has an existing CytoNorm availability assertion failure,
   reproduced on untouched HEAD; do not call full validation green.
   **Done when:** one deliberately broken candidate leaves a `failed` record,
   identity and another candidate finish, the report names the failure, and no
   invalid or stale corrected output is scored.
   **Reach:** adapters/base/validate, pipeline, candidate output handling, stage-level
   integration tests. Evidence: `pipeline.py:256–347`, `adapters/validate.py:19`.
   (L1-03, L10-08.)

2. **P02 — Make R results a checked boundary.**
   **Status: implemented; scoped Python and R exchange checks passed.**
   [P02 implementation and limits](P02_R_BOUNDARY.md) and
   [receipt](prototype-evidence/2026-10-10-p02-verification.json).
   Require a valid result payload after a successful R exit. Preserve empty and
   one-element lists as lists, including warnings and SOM node IDs; the existing
   `auto_unbox` path can turn strings into character-by-character entries.
   Reuse the existing conversion/validation mechanisms rather than inventing a
   new RPC layer. Resolve a bounded timeout and record timeout/log outcomes.
   **Done when:** zero, one, and multiple warnings/node IDs round-trip correctly;
   missing or invalid `result.json` and a timeout become the appropriate failure.
   **Reach:** `r_runner.py`, affected R writers and Python consumers, boundary tests.
   Evidence: `r_runner.py:99–146`, `cycombine.py:215–218`.
   (L1-04/06/14, L2-14, L10-08.)

3. **P03 — Declare the files the workflow actually consumes and produces.**
   The validate rule depends only on YAML, while stages reread data and manifests.
   Correction H5AD is neither a declared correction output nor a metric input.
   Represent each candidate's output bundle in the DAG so both a successful H5AD
   and an explicit unsuccessful outcome are legitimate; do not require fake H5AD
   files for failed candidates. Declare labels, thresholds, sample/marker/processing
   files, input data, and relevant execution identity. Move exploratory override
   records to per-candidate artifacts so parallel corrections cannot lose entries;
   clear/tie them to the current analysis when eligibility changes.
   **Done when:** changing an input invalidates its descendants, deleting a successful
   corrected H5AD repairs it, a no-op does no analysis, and two overrides both survive.
   **Reach:** Snakefile, RunLayout, override reader/writer, workflow resilience tests.
   Evidence: `workflow/Snakefile:67–85,149–184`, `pipeline.py:131–161,350–368`.
   (L1-06, L7-09/10/11/12.)

4. **P04 — Escape report metadata and repair sharing/redaction behavior.**
   `select_autoescape(["html"])` does not match `report.html.j2`; arbitrary metadata
   currently enters HTML unescaped. Turn escaping on while retaining deliberately
   trusted SVG and correctly encoded JSON. Keep the reversible identifier map
   outside the shareable report bundle. Make overlapping-ID replacement
   deterministic and prevent private filenames/free text from leaking identifiers.
   Use a persisted per-run salt if stronger cross-study separation is implemented;
   generating a new salt for each render would break repeatability. Describe the
   result as pseudonymisation, not anonymisation.
   **Done when:** HTML metacharacters render as text, the intended SVG/JSON still
   works offline, S1/S10 redact consistently, repeated reports match, and a report
   export excludes its reversal map.
   **Reach:** reporting/html/redact/template, report integration tests, sharing docs.
   Evidence: `reporting/html.py:532–542`, `redact.py:24–84`.
   (L10-01/02/03.)

5. **P05 — Refuse unsafe initialization and incompatible configuration/results.**
   `init` currently writes into an existing destination without an overwrite
   guard. Contract records stamp but do not enforce `schema_version`. Analytical
   commands should validate configuration consistently before launching a workflow,
   reject unknown backend names and invalid resource values, and give actionable
   errors. Reject inactive preprocessing settings on the prototype route rather
   than accepting them as if they affect a run. Avoid a mandatory new migration
   framework; bump versions only where a changed contract requires it.
   **Done when:** non-empty initialization leaves existing files intact; an
   unsupported schema version is refused; a typo fails before any correction;
   unsupported raw-preprocessing settings name the unsupported operation.
   **Reach:** CLI/examples/config/contracts, input validation, focused tests.
   Evidence: `examples.py:180–201`, `contracts.py:31–49`, `config.py:346–348`.
   (L1-05/07, L6-02, L10-07.)

### M2 — Full OpenProblems metrics and honest decisions

1. **P06 — Fix and wire the six existing compatibility scores.**
   Correct the EMD public-wrapper baseline/candidate dispatch. Reconcile EMD,
   average batch R2, LISI, and peak consistency with the actual upstream two-split
   input, filtering, pairing, aggregation, and normalization contract. Existing
   helper parity is not proof of full component equivalence: the current local
   grouping adaptations must be removed or explicitly validated and labeled.
   **Done when:** both identity and actual correction run these six named scores
   through the normal workflow, respond to changed candidate values, and match
   pinned component outputs on the same eligible inputs within recorded tolerances.
   **Reach:** openproblems metrics, config/input adapter, pipeline, public-entry
   integration/parity tests. Evidence: `openproblems.py:337–368`,
   upstream `src/api/comp_metric.yaml` and component scripts. (L4-01/03/10, L8-08.)

2. **P07 — Implement FlowSOM mapping similarity as a required score.**
   Use the pinned bidirectional donor-paired split mapping definition, lineage
   clustering channels, seed, grid dimensions, and cluster count. Reuse available
   FlowSOM/R boundary infrastructure; verify the needed environment rather than
   inheriting the stale package-absence deferral. Include donor/direction detail
   and the exact parameter provenance.
   **Done when:** `flowsom_mean_mapping_similarity` executes in the user workflow,
   matches the pinned full component on shared inputs, and refuses unsupported
   pairings or missing lineage/parameter data with a named result record.
   **Reach:** metric/R adapter, config, typed results/report, parity tests.
   Evidence: upstream `src/metrics/flowsom_mapping_similarity/{script.R,helper.R}`.
   (L4-09/11; previously deferred V01.)

3. **P08 — Deliver the nine-score contract and OpenProblems/scIB rank visuals.**
   Declare the roster from the pinned component metadata; retain upstream score
   directions and grouping/scale semantics, including normalized cLISI. Every
   candidate has all nine records plus identity comparison where meaningful;
   errors/unsupported/undefined states never disappear from the table. Record
   input/split identities, exclusion counts, sample/donor coverage, implementation
   version, and parity provenance. Custom MAD/joint/frequency/native metrics may
   remain supplemental diagnostics; do not make additional custom scoring families
   a prerequisite or substitute them for missing OpenProblems scores.
   Add the ranked table, metric glyphs/category/overall bars, per-metric
   normalization QC, offline details/sorting, and SVG/PDF/PNG exports described in
   [BENCHMARK_VISUALIZATION_SPEC.md](BENCHMARK_VISUALIZATION_SPEC.md). Execute and
   record benchmark control methods for metric-specific scaling; verify and pin
   the reference renderer's scaling rules. Keep artificial benchmark controls
   distinct from physical control samples and user-selectable corrections.
   Preserve raw values and fixed anchors; only fully evaluated candidates receive
   an overall rank on a declared common scope. Partial/failing rows stay visible.
   **Done when:** report, CSV, and JSON agree on all nine names/statuses/scopes;
   an incomplete score cannot masquerade as suite coverage; optional diagnostics
   remain clearly separate; rank visuals and normalization QC reproduce the same
   data with correct directions, ties, complete coverage and scaling provenance.
   **Reach:** registry/contracts/report/plots/template/workflow/status/config,
   complete-roster, ranking/scaling, and export integration checks.
   (L1-01, L4-03/04/10/15; user visualization requirement.)

4. **P09 — Implement both upstream functional-marker preservation scores.**
   Reuse verified per-sample arithmetic where appropriate, but implement the pinned
   two-batch baseline significance intersection, consistent group orientation,
   and evaluation across both integrated splits. Wilcoxon is retained baseline
   pair fraction; Cohen's d is mean absolute effect-size drift from baseline.
   The package's isolated donor-level contrast statistics are not replacements.
   Preserve upstream significance and aggregation rules; report unsupported design,
   no functional markers, or an empty significant baseline as non-numeric statuses.
   **Done when:** both named scores execute end to end and match pinned component
   outputs on eligible inputs, with destructive and retained-contrast controls.
   **Reach:** metric implementation/input adapter/pipeline/report, component parity
   and edge-case tests. Evidence: upstream functional-marker `script.py/helper.py`.
   (L4-03/12/13, L5-04; previously deferred V01.)

5. **P10 — Assess candidates using declared OpenProblems evidence and scope.**
   Update the policy integration: it currently reads only custom EMD/MAD/joint
   metrics. Expose the complete upstream evidence and the declared descriptive
   benchmark ranking without inventing universal approval thresholds. Required score failures, unsupported
   comparisons, missing split metadata, and partial correction must constrain any
   favorable claim. Keep explicit eligibility/identity/failure checks and provisional
   review outcomes; calibrating approval remains separate. If supplemental guards
   are used, missing guard values must still fail closed.
   **Done when:** the decision/report consumes the canonical nine-score results;
   dropping required evidence cannot improve an assessment; an exploratory override
   cannot be approved; no automatic recommendation rests on uncalibrated thresholds.
   **Reach:** policy/pipeline/config/results/report and reachable-policy tests.
   Evidence: `policy.py:231–323`, upstream component metadata. (L5-03/05/10/12/14.)

### M3 — Correction and diagnostic semantics

1. **P11 — Make biological protection explicit in cyCombine.**
   `design.protect` is used by the audit, but cyCombine only uses the separately
   supplied `covar`. Define the one supported covariate choice when a protected
   factor is declared; refuse or clearly warn about unsupported multiple factors
   instead of promising protection. Record when upstream drops covariate protection
   inside SOM nodes. Do not treat auditing a protected factor as protecting it in
   the correction model.
   **Done when:** a supported protected factor reaches R, a missing factor refuses,
   and a synthetic biological contrast remains visible after an actual cyCombine
   correction; the shipped candidate configuration matches the documented behavior.
   **Reach:** config/examples, cyCombine adapter/R script, actual-backend tests.
   Evidence: `cycombine.py:76,133–145`, `examples.py:120`,
   `cycombine_correct.R:65–76`. (L2-05/08/09.)

2. **P12 — Honor the correction/clustering masks and only working parameters.**
   Backends cluster on correction markers despite the separate clustering mask.
   Pass the intended clustering subset separately and retain the normalization
   mask; verify excluded markers stay unchanged. Reject cyCombine `anchor` until
   its required input is actually supplied, or implement and test that input.
   Record effective defaults in configuration/specs even for ineligible outcomes.
   **Done when:** changing only the clustering mask changes the intended upstream
   call; invalid parameters fail early; no accepted option silently does nothing.
   **Reach:** adapters/R scripts/config/example manifests, mask and parameter tests.
   Evidence: `cycombine.py:121,158–178`, `cytonorm.py:311,379–397`.
   (L1-10, L2-04/06.)

3. **P13 — Implement benchmark split execution and distinguish fit scope.**
   Freeze the upstream two-split donor/sample/batch layout and correction runs,
   including documented control inclusion for fitting and exclusion for evaluation.
   Store and score two integrated outputs per candidate rather than merely labeling
   one joint correction as two independent runs. The benchmark split contract is
   separate from train/tune/validate; neither automatic donor-disjoint splitting nor
   random cell halves reproduces its replicate design by itself. Retain honest
   joint-fit versus fit/apply labels for cyCombine/CytoNorm and describe any validated
   grouping adaptation explicitly.
   **Done when:** candidate input/output cell IDs and scope match both frozen
   benchmark splits, donor/batch pairings are valid, applicable fit/apply evaluation
   does not leak training units, and identity follows the same score path.
   **Reach:** sampling/config/input adapter/candidate output layout, Snakefile,
   metrics/policy/report. Evidence: upstream `src/api/file_unintegrated.yaml`,
   `src/data_processors/process_dataset/script.py`; `cycombine.py:62–63`.
   (L2-07, L5-01/06/12.)

4. **P14 — Measure correction applicability and disclose uncorrected regions.**
   A successful process can leave batches, samples, or SOM nodes untouched.
   Record supported/unsupported cells and samples and covariate-dropped nodes;
   distinguish applicability from numerical change, since a valid correction may
   legitimately leave a value unchanged. Use those records in assessment and
   report coverage instead of relying only on warnings or a success status.
   **Done when:** an intentionally unsupported region cannot disappear from a
   denominator or support a whole-dataset conclusion; a correct identity-like
   result is not classified as failed merely because it changes nothing.
   **Reach:** adapter payload/result, policy, report, partial-correction fixtures.
   Evidence: `cytonorm.py:482–545`, `cycombine_correct.R:65–76`.
   (L2-02/09/12, L5-05.)

5. **P15 — Make diagnostic flags interpretable on small real designs.**
   Composition permits two samples per batch; absence of a significant result
   must not imply adequate power. Add a support rule or an appropriate finite-sample
   test, retain descriptive frequencies, and disclose multiplicity for the tested
   family. Label CytoBatchFlagR fences as uncalibrated exploratory flags with the
   number of assessments; do not convert them to probabilities or evidence of
   absence. Distinguish actual acquisition time from an event-order proxy; enable
   a step flag only with a declared threshold. Show available density/composition
   plots and diagnostic status, without making new embeddings a prerequisite.
   **Done when:** a thin design yields an insufficient/descriptive result, flags
   show their thresholds/origins and family size, and event-index drift is labeled
   correctly. **Reach:** detectors/context/config/report, negative and thin fixtures.
   Evidence: `composition.py:20,83–117`, `io/dataset.py:422–428`,
   `diagnostics/drift.py`, `diagnostics/cytobatchflagr.py`.
   (L3-01/02/07/08/09/10/13/14/15.)

### M4 — Measured-data use, tests, installation, and documentation

1. **P16 — Freeze a real input eligible for the full OpenProblems suite.**
   Select an upstream benchmark dataset or equivalent measured cohort with genuine
   donor/sample identities, cross-batch replicate layout, stable cell-type labels,
   exactly supported biological groups, lineage/functional marker roles, correction
   mask, and benchmark split assignments. Verify access/rights, checksums, transform,
   processing provenance, and FlowSOM parameters. The existing control-only spectral
   subset can support supplemental technical analysis but cannot establish full-suite
   coverage: upstream removes controls and these files have no established donors.
   Do not fabricate donor or split metadata to make the suite run.
   **Done when:** frozen H5AD/manifests support all nine scores with a defined
   functional baseline and a reproducible preparation recipe.
   **Reach:** data metadata/preparation/config, upstream-format adapter and input
   eligibility tests. Evidence: upstream metric API/file schema, spectral DATASET.md.
   (L6-04, L8-04; full-suite coverage promoted from V04.)

2. **P17 — Measure the complete suite on actual corrections.**
   Before looking at candidate scores, freeze the real-data benchmark scope,
   endpoints, pairings, filters, and limitations. Run identity and actual cyCombine
   on both benchmark splits and measure all nine scores; retain per-score detail,
   canonical JSON/CSV, both integrated outputs, report, hashes, logs, failures, and
   resource receipts. Negative results remain valid software/scientific observations
   if faithfully reported. A genuinely undefined baseline must not be replaced with
   a fabricated value; identify another eligible input before closing full coverage.
   **Done when:** another researcher reproduces all nine eligible measured scores,
   with upstream parity evidence and appropriately bounded interpretation. This does
   not establish universal policy thresholds or validation of every modality.
   **Reach:** full-suite run configuration, immutable validation artifacts and docs.
   Evidence: this measured run has not yet been executed. (L8-03/04/11.)

3. **P18 — Test through the supported CLI and pipeline, using actual methods.**
   Add full nine-score component/public-pipeline parity and integration checks for
   candidate-sensitive metrics, failure containment,
   preservation reachability, covariate/mask behavior, missing-evidence refusal,
   and status semantics. Keep idealized and destructive fixtures as controls, but
   supplement them with actual cyCombine; a test-only oracle is not its validation.
   Exercise a clean run, no-op, artifact deletion, input change, report regeneration,
   and interruption/resume. Run the changed stochastic claims at several declared
   seeds; an entire five-seed scientific suite is a later stronger gate.
   **Done when:** public entry-point regressions fail for the confirmed bugs, all
   scoped integration checks pass, and each omitted check is disclosed.
   **Reach:** current unit/scientific/integration suites, seed fixtures, gate registry.
   (L8-01/02/03/06/07/08/09/10/13/14.)

4. **P19 — Verify one installation path outside the checkout.**
   Support and document Python 3.12 plus one executed R dependency environment
   for the prototype. Keep the pinned backend-test runtime distinct from standalone
   R analysis in the user's R4_51 environment; any runtime change needs parity
   verification. Record actual Python/R/package versions rather than assuming the
   environment YAMLs are what runs. Build/install a wheel in isolation and exercise
   init/validate/diagnose/benchmark/report from another directory. Either ship the
   audit scripts as package assets or explicitly label `audit` checkout-only.
   **Done when:** a clean install runs the documented workflow and finds all assets;
   the supported path has no hidden imports or dependence on the template's pixi
   setup. **Reach:** packaging/paths/CLI/environment docs, install smoke test.
   Evidence: `pyproject.toml:73–89`, `cli.py:253–255`,
   `.github/workflows/batchdx-cyto-tests.yaml`. (L7-07/08/20, L10-06/07/10.)

5. **P20 — Finish status, documentation, and the complete validation record.**
   In diagnose/selected-correction mode, unrequested candidates must be skipped,
   not failed; stale outputs from another mode must not establish execution. Make
   documented logs contain useful content. Repair getting-started paths, base-dir
   instructions, availability requirements, contradictory execution claims,
   scope/fit labels, and enabled-tool notices. Explain corrected H5AD location and
   expression layer. Generate status from a complete current report, including
   the slow workflow checks and R exchange suite, with matching lint/audits and CI
   evidence. Do not relabel library-only or calibrated-on-synthetic capabilities
   as workflow/scientific validation. Tag/PyPI publication is not required for
   a local usable prototype and needs its own explicit authorization.
   Document the benchmark figure's controls, metric grouping, equal-metric overall
   mean, complete-scope ranking, and exports. A real measured run must demonstrate
   the requested visuals; a synthetic layout preview alone cannot close P08/P20.
   **Done when:** a new user follows the instructions successfully; a report-only
   rerender causes no correction; every scoped gate has a current receipt; known
   deferred functionality stays labeled. **Reach:** status/CLI/docs/examples,
   notices/citations, gates/status generator/CI. (L7-13/14, L8-05, L9-01–22,
   L10-04/07/11.)

## 4. Extensions after the scoped prototype

These are genuine remaining package tasks. Their priority changes if the prototype
is advertised to support the corresponding input, backend, claim, or deployment.

### C — CytoNorm as a supported prototype backend (about 2–4 additional days)

1. **C01 — Require scientifically valid training-control identity/coverage.**
   `_training_controls` checks role and batch coverage, not shared `reference_id`.
   Define valid common control material for the fitted batches; disconnected
   controls must not silently become interchangeable. Reuse the reference graph
   where appropriate, but graph connectivity alone does not prove every material
   can be pooled into the same quantile target. Test unrelated control materials.
   (L2-01; `cytonorm.py:236–288`, `design/audit.py:190–192`.)
2. **C02 — Classify design refusals as ineligible and agree with the design audit.**
   No training control, inadequate training-batch coverage, unbalanced splits,
   and leakage should not be reported as missing software. Test the actual frozen
   training layout, not only full-dataset control counts. (L2-03/10, L5-02.)
3. **C03 — Supply a genuinely runnable reference example.**
   Replace the insufficient-control demonstration with enough distinct reference
   units for the declared split, without fabricating extra independent donors.
   Require successful train-only fit, validation application, and an explicit
   coverage table. Keep the thin-design refusal example separately. (L9-10/11.)
4. **C04 — Honor clustering/normalization channels and supported targets.**
   Pass the intended `colsToUse`; test any exposed goal/range parameters. A reusable
   frozen-model/new-batch workflow is a separate capability: implement persistence,
   model/input compatibility, provenance, and apply tests before advertising it.
   (L2-06/13 and fit/apply gap.)
5. **C05 — Make float32 FCS accounting consistent.**
   Use measured round-trip tolerance for unchanged-sample/marker accounting;
   distinguish serialization drift from correction. Validate identity-function
   clusters and unfitted samples as partial coverage, with end-to-end report checks.
   (L2-02/11/12, L9-12; `cytonorm.py:466–545`, `io/fcs.py:37`.)

### R — Automated raw-FCS-to-report workflow (about 4–7 additional days)

1. **R01 — Add preprocessing before validation.**
   Conditional rule/stage for declared gates, transform, and PeacoQC; write the
   immutable baseline plus processing/removal/QC artifacts. Keep the preprocessed
   bypass intact. Define supported compensation/unmixing prerequisites. (L6-02.)
2. **R02 — Repair removal provenance and accounting.**
   The loader currently demands removed events be present among surviving baseline
   events. Validate them against original source counts/checksums/identity instead;
   require kept + removed to match source events and reject overlap/duplicates.
   Accept a valid empty removal log. Show removals and per-sample QC failures.
   (L6-05/14; `io/dataset.py:310–356`.)
3. **R03 — Validate FCS semantics and actual acquisition time.**
   Surface and test Time/TIMESTEP, amplification/gain, channel metadata, and stored
   compensation/unmixing state. Refuse unsupported representations rather than
   silently treating raw stored values as analysis-scale values. Confirm PeacoQC's
   time/scale behavior against its installed implementation. (L6-01/03/07.)
4. **R04 — Make channel and transform mappings explicit.**
   Check both channel names and antigen labels per file; resolve the inert
   `fcs_glob` contract; validate per-channel/per-panel cofactors, forward/inverse
   convention, gating coordinate scale, and the correction/clustering/evaluation
   masks. Do not harmonize the six spectral layouts just by matching channel names.
   (L6-04/06/08/11.)
5. **R05 — Add realistic raw-route integration fixtures and a measured run.**
   Include negatives, scatter/time, metadata mismatches, supported transform
   variation, zero removals, and QC failure cases. Prove event identity, no double
   processing, immutable baseline across candidates, and input-change/rerun behavior.
   FCS export stays unsupported until the separate inverse/metadata/export path
   passes; a callable export helper does not establish a supported CLI feature.
   (L6-12/13/14.)

### V — Stronger scientific evaluation (scope in separate bounded plans)

1. **V01 — Validate supplemental metrics only when they remain exposed.**
   All nine OpenProblems scores are now required by P06–P09, not deferred here.
   Native cyCombine identity-only scoring and any two-run flag metric need correct
   candidate inputs and tests if advertised as supplemental diagnostics. Keep CMS
   absent at the existing pin; expanding to a different upstream version needs an
   explicit roster/contract review. (L4-02/09/10/11.)
2. **V02 — Add supplemental independent-unit uncertainty and distortion summaries.**
   The custom Stage 8 donor-contrast workflow remains optional supplemental
   evidence, not the upstream functional preservation implementation. Add
   sample/donor-balanced technical metrics, independent-unit intervals, replicate
   agreement, finite-sample/null floors, scale-aware thresholds, rare populations,
   per-marker worsening, and multiplicity across declared contrasts. Validate
   incomplete paired/unpaired cases, bootstrap strata, successful-resample
   denominators, and sample-versus-donor labels. Avoid a large metric catalog until
   each new metric answers a prespecified scientific question. (L4-05–13.)
3. **V03 — Strengthen confounding and comparison-scope diagnosis.**
   Audit partial imbalance, design rank, disconnected controls, nesting, and
   metadata completeness; connect these to contrast-level eligibility. Preserve
   unsupported batches in whole-run coverage. Require an executed test for any
   claimed stronger identifiability criterion. (L3-03, L5-05/08/09/12.)
4. **V04 — Validate distinct data roles and modalities separately.**
   Obtain suitable technical replicates and independent biological-contrast
   cohorts with rights/provenance checks, measured donor identities, stable labels,
   and supported designs. Spectral controls alone cannot validate conventional
   flow, CyTOF, or biological protection. Accessions and literature candidates
   suggested by unfinished external reviews require verification before selection.
   (L8-03/04/11; capability matrix modality validation rows.)
5. **V05 — Calibrate without reusing evaluation data.**
   Freeze endpoints; estimate thresholds/flag error behavior on tune data; test
   once on untouched validation units with uncertainty and failure denominators.
   Run the complete relevant stochastic suite at five declared seeds. A configurable
   threshold and a “calibrated” string are not a calibration artifact. Freeze
   policy calibration provenance before allowing recommendations. (L5-07/14,
   L8-01/07/10/11; old handoff C1.)

### D — Reproducible distribution and HPC (not needed for local prototype)

1. **D01 — Lock the executed environment and test support claims.**
   Choose one source of truth for Python/R dependency versions; distinguish the
   micromamba prefix from R `renv`; test claimed Python/dependency bounds or narrow
   claims. Check dirty upstream trees as well as SHAs, tool-download integrity,
   CI permissions, and action pins. Quarantine/template-scope unrelated pixi and
   CI jobs rather than redesigning the repository as a prototype prerequisite.
   (L7-07/08/19/22, L8-12, L10-05/06/09.)
2. **D02 — Build and verify the actual container dependency closure.**
   Align it with the executed R recipe, make missing-package installs fail the
   build, resolve needed system libraries through an actual build, pin versions,
   run adapter smokes, and record image/SIF digest/checksum. A digest-pinned base
   is already in source; stale prose that calls it only a tag must be corrected.
   (L7-02/03/04/05/06; `containers/Dockerfile.r-backends`.)
3. **D03 — Make a site profile usable and execute it.**
   Replace placeholders, verify current site/account settings, activate dependencies
   on compute nodes, align CPU/thread declarations, and reject unresolved template
   values. Document direct Snakemake invocation; add CLI profile support only if
   part of the supported user path. Capture a real smoke job's ID/logs. (L7-15/16.)
4. **D04 — Verify resources and recovery on the target execution platform.**
   Prefetch for offline nodes, pass thread limits to real computation, prevent
   nested BLAS oversubscription, isolate run state where needed, write meaningful
   logs, measure peak memory/walltime, and compare an interrupted/resumed benchmark
   with a clean run. Local synthetic recovery is not a cluster receipt.
   (L7-13/17/18; old handoff C3.)
5. **D05 — Prepare a release from final executed evidence.**
   Fresh wheel/sdist checks, notices/CFF validation, compatible versions and assets,
   current generated gates, consistent docs, complete CI, and retained validation
   artifacts. Do not require every full-package extension for a honestly scoped
   Alpha. Commit/push/tag/publish remain separate authorized actions; this review
   performs none of them. (L10-04/06/12, RELEASE.md, old handoff A3/C4.)

## 5. Capability choices beyond defect repair

These gaps matter for a broadly complete package, but they are not all necessary
for the first prototype. Confirm demand and verify maintained APIs/primary sources
before adding a dependency or algorithm. The incomplete Biomni summaries are not
adjudicated evidence and do not mandate implementation.

1. **More correction methods:** one simple reference alignment method, a
   reference-free landmark method, or a replicate-based method could add design
   coverage after the existing methods work correctly. Select a specific user
   design first; do not implement the entire external benchmark method roster.
2. **Modality-specific preprocessing:** CyTOF bead normalization/debarcoding and
   specialized cleanup; flow compensation consistency; spectral unmixing/control
   and autofluorescence provenance. Either integrate a verified existing tool or
   define upstream preprocessing as an explicit supported prerequisite.
3. **Richer diagnostics:** sample-level PCA/MDS and batch-versus-biology summaries,
   per-sample outliers, cross-file run-order drift, reagent/lot effects, and native
   marker-distance visualization. Add modality/platform/run metadata only when
   used; separate a biological composition difference from technical intensity drift.
4. **More panels and downstream exports:** shared-marker mapping across layouts,
   explicit excluded-marker/channel roles, supported FCS inverse/export semantics,
   reusable models and incremental/new-batch workflows. Each needs its own contract
   and measured-data acceptance checks.
5. **Embedding or graph integration:** introduce a separate output type and
   evaluation path if needed. Do not force an embedding into the current
   same-shape, same-scale corrected-expression contract. Label transfer,
   neighborhoods, and rare-population checks are optional new evidence families,
   not reasons to rebuild the existing registry before a prototype.

No prototype work is required solely to modularize the metric dispatcher, add
`py.typed`, make version metadata dynamic, implement static checkpoints, remove
every template artifact, or fix reviewer citation-line nits. Retain the existing
version-consistency tests. L9-19 records verified existing behavior, not a defect.
The existing lens files retain individual nits and inferred claims for later
adjudication; this list consolidates overlapping findings into outcome-based tasks.

## 6. Acceptance and first bounded action

The prototype is usable when all five outcomes hold:

1. A clean installation runs a documented single-panel dataset through the CLI,
   producing an inspectable corrected H5AD and offline report.
2. Identity and actual cyCombine use the intended masks/covariates and both
   benchmark splits; all nine upstream scores execute with parity evidence, and
   failed/partial/ineligible outcomes stay visible.
3. Metrics report their scope, independent units, coverage, missing evidence,
   and exploratory limits; no uncalibrated automatic approval is possible.
4. A real full-suite eligible benchmark/cohort run has reproducible measured
   outputs for all nine scores and appropriate scientific boundaries; control-only
   or fabricated-donor demonstrations do not close this gate.
5. Complete scoped tests, R exchange, lint/audits, install checks, and current CI
   receipts agree with the generated status and documentation; no stale results
   survive input changes or missing outputs.

**Next implementation ticket:** P03, declare consumed files and candidate output
bundles in the DAG. Start with input-change/deleted-H5AD regressions, including
failed candidates that correctly have no H5AD. P01/P02 scoped acceptance checks
pass; broader CytoNorm availability/guard failures remain recorded for validation
recovery. Earlier time estimates are planning history, not measured remaining time.

## 7. Evidence and limits of this review

The central workflow, public metric wrapper, common metrics, decision policy,
adapter boundaries, CytoNorm training-control selection, configuration/records,
report escaping/redaction, initialization/packaging, existing review lenses,
data metadata, and local validation reports were inspected.

Fresh executable checks and their exact scope are recorded in
[`prototype-evidence/2026-10-10-verification.json`](prototype-evidence/2026-10-10-verification.json).
The probes confirm: escaping is false for the actual template suffix; schema
version 999 is accepted; missing guardrails can yield a recommendation under a
test-only policy; and the public EMD wrapper dispatches the baseline panel to
both EMD functions. The EMD probe isolates dispatch by mocking component execution;
it is not an upstream numerical parity test. The test-only policy is not exposed
through shipped YAML, whose policies remain uncalibrated.

At the original review checkpoint, before subsequent P01 implementation, no
source fixes, dependencies, full R/backend run, measured-data correction,
container build, scheduler submission, external-review continuation, GitHub
mutation, or release publication were performed. Current remote CI was not queried;
saved CI/handoff states are historical. Detailed external-method claims remain
unverified. The prior full-package roadmap and unfinished Biomni review are not
silently declared complete by this prototype-focused audit.

The 2026-10-10 user scope update expands the metrics gate. The earlier 8–13-day
scoped estimate is superseded by 13–20 days plus data-access/preparation delays.
The upstream checkout was inspected locally; online retrieval of the pinned
GitHub sources failed. No metric implementation or full-suite execution was
performed during this documentation revision.
