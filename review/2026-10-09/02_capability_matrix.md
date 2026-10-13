# 02 Capability matrix: target vs implemented vs real-data vs gap (2026-10-09)

STATUS: COMPLETE (sections A-E written 2026-10-09).

Repo HEAD `b313006` + uncommitted doc edits. Read-only review. Sources for "implemented": sibling lenses
`lenses/L1..L6` (cited as L#-NN, each of which cites path:line), `docs/validation_results.md`, `RELEASE.md`,
`docs/implementation_status.md`, `data/DATA_MANIFEST.md`, plus targeted greps done for this file. Nothing was executed.

## How to read this

- **Row = one concrete capability.** A `Mod` column says which of M (mass cytometry / CyTOF), F (conventional flow),
  S (spectral flow) the target applies to. Capabilities that need modality-specific physics get their own row; a row
  marked M/F/S is the same code path for all three (the repo has no `modality` field at all, so nothing in the code is
  conditioned on modality: L3 coverage matrix; `grep -i "modality\|CyTOF"` in `schemas/contracts.py` returns nothing).
- **implemented** uses three words and never merges them:
  `exists` = code in `src/` or `workflow/scripts/`; `wired` = reachable from `workflow/Snakefile` -> `pipeline.py` in a
  user run; `tested` = a pytest/testthat test touches it (always on simulated or fixture data, see next bullet).
  `yes` below means exists+wired+tested; anything less is spelled out.
- **tested on real data**: the repo has NO real-data result for any capability. `data/raw` and `data/processed` are empty;
  the only real data is a downloaded, checksummed, audited public spectral technical-control set (Zenodo 15388817; 72 FCS,
  27,022,765 events; 24-file/6-batch homogeneous subset) with **no workflow run on it**
  (`data/metadata/cytobatchflagr_spectral/DATASET.md`; `docs/validation_results.md:24-26,109-110,132-137`;
  `RELEASE.md` header "real-cohort calibration ... not verified"). So every cell in the column is "no" unless
  stated; there is no mass-cytometry or conventional-flow real data anywhere in the repo.
- **Priority**: P0 = scientific validity (can produce a wrong/unsupported conclusion), P1 = usability / completeness a user
  would expect, P2 = nice to have.
- Citations in `target` are keyed `[Rn]` and listed in section D with title/DOI/URL. "unverified" = I could not confirm
  the claim from a primary source during this review. All literature searches used generic method names only.

## Headline

1. Wired in a user run today: identity, cyCombine and CytoNorm correction; within-file drift, composition (Kruskal),
   CytoBatchFlagR and cyCombine-EMD detectors; 6 "common" metrics; decision policy; HTML/JSON report. All on simulated data only.
2. Code exists but is **not wired**: raw-FCS gating / asinh / PeacoQC (`preprocess/`, L6-02), the 6 OpenProblems-compatible
   metrics incl. iLISI/cLISI (`metrics/openproblems.py`), and all Stage 8 donor-level preservation statistics
   (`metrics/preservation.py`) (L4-03, L1-01). Native cyCombine EMD/MAD run only identity-vs-identity (L4-02).
3. Absent: any bead/drift (mass) or compensation/unmixing (flow/spectral) handling, any modality field, reference-free flow
   landmark methods, replicate-based RUV, embedding methods, label-transfer / kBET / silhouette over-correction checks,
   sample-level uncertainty on decision metrics, real-data calibration of any threshold.

## A. Cross-modality capabilities (Mod = M/F/S unless stated)

### A1. Pre-processing and QC (feeds diagnostics; "area" = diagnostics)

| ID | Mod | target (what a complete package should have, with citation) | implemented (yes/partial/no; exists vs wired; path:line) | tested on real data | gap + priority |
|---|---|---|---|---|---|
| A1.1 | M/F/S | Event-level signal-stability QC: remove unstable-flow / margin / outlier events per file, report what was removed. PeacoQC (mass + flow + spectral) [R12]; flowAI [R13] | **partial**. Exists: `preprocess/peacoqc.py:393` `run_peacoqc` + `workflow/scripts/peacoqc_qc.R`, removal log `preprocess/removals.py`. **Not wired**: no Snakefile rule or pipeline caller (L6-02; `preprocess/__init__.py:3-6` "deliberately outside the run"). flowAI: absent. | no | P1: wire raw-FCS stage (handoff B2, 3-5 d). Until then the package cannot start from raw FCS. |
| A1.2 | M/F/S | Gating of debris/doublets/dead cells with an auditable gate table (FlowJo/openCyto-style) [R14 unverified]. | **partial**: exists `preprocess/gating.py:405` `apply_gates` (polygon/rectangle gate tables), not wired (L6-02). No automated gating. | no | P1 (same fix as A1.1). Automated gating is P2. |
| A1.3 | M/F/S | Transform to a variance-stabilised scale with per-channel cofactor (asinh; cofactor 5 mass, ~150-ish flow; CATALYST docs [R7]); cofactor choice is a documented, checked parameter. | **partial**: exists `preprocess/transform.py:154` `asinh_transform` with convention checks, not wired; the correction pipeline refuses to transform (`io/dataset.py:143-152`). No cofactor estimation. | no | P1: wire with B2; add per-modality cofactor guidance/estimation (spectral negatives and wide range untested: L2 gap 6). |
| A1.4 | M/F/S | Acquisition **time** read from FCS (`Time`, `$TIMESTEP`) so drift is vs real time, and run order across files. | **partial/incorrect label**: `acquisition_order = np.arange(n)` event index; Time channel never read (L6-01; `io/dataset.py:422-428`; L3-07). | no (spectral subset has a time channel: DATASET.md) | P1 (0.5-1 d). Becomes P0 for any drift claim. |
| A1.5 | M/F/S | Ingest with channel/marker roles (lineage / functional / excluded / clustering) and a **modality** field so checks can be conditioned. | **partial**: roles exist; `use_for_clustering` mask ignored by both backends (L2-06); no modality field (L3-11, grep above). | no | P1: add `modality` + `platform` to manifest (prerequisite of every modality row below). |

### A2. Diagnostics (pre-correction: is there a batch effect, what kind, is it identifiable?)

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| A2.1 | M/F/S | Sample-level visualisation of batch structure: per-sample marker-median PCA/MDS/UMAP coloured by batch and biology; MDS/UMAP and density plots (cyCombine `detect_batch_effect` provides MDS, UMAP, EMD heatmap) [R1]; CytofRUV diagnostic plots [R3]. | **no (in the report)**. cyCombine's own PNGs are generated via R but the report prints their paths only (L3-04; `reporting/templates/report.html.j2:278-280`). `plots.drift_lines` defined, never called. | no | P1 (3-5 d, Python-native so it works without R). |
| A2.2 | M/F/S | Marker-wise batch-shift quantification: EMD/KS per marker x batch pair with uncertainty and a technical-replicate null [R1 uses EMD; R2 uses EMD]. | **partial**: EMD only via R `cyCombine::compute_emd`, table only, no flag, no threshold, pooled over cells, doubled rows (L3-05, L3-14; `diagnostics/cycombine_detect.py:134-148`). exists+wired; unavailable without R. | no | P1/P0-leaning: unthresholded numbers invite over-reading; Python-native `marker_shift` with sample-level bootstrap (3-4 d). |
| A2.3 | M/F/S | Technical-reference-based flagging (control sample across batches; outlier batch detection): CytoBatchFlagR [R15 unverified title/DOI; repo pins upstream commit]. | **yes (wired)** `diagnostics/cytobatchflagr.py`, `workflow/scripts/cytobatchflagr_detect.R`. But 1.5-IQR fence has measured FPR 0.207/test at 4 batches, not shown in report (L3-01). | **no run on the real spectral set** (set was audited for this purpose; no workflow run) | **P0**: attach expected-false-flag rate and n_tests to every assessment (1-2 d); calibrate on the 18-batch technical-control null (3-5 d). |
| A2.4 | M/F/S | Composition / frequency shift across batches with compositional statistics and multiplicity control (CLR / count models, e.g. crumblr-type; unverified). | **partial**: Kruskal per population, no multiplicity, power-gate bug at 2 samples/batch (L3-02; `diagnostics/composition.py:20,83-117`); frozen-label tier may be batch-driven (L3-13). wired. | no | **P0** (flag can be silently "ok, no flag" when test cannot fire): raise floor / BH / CLR (1-2 d). |
| A2.5 | M/F/S | Graded batch-vs-biology confounding diagnosis (design-matrix rank, Cramer's V, variance partition of sample medians) [R1 states covariate/biology must not be confounded; unverified exact wording]. | **partial**: binary declared-metadata perfect-confounding check only; `design_full_rank` computed but not used for identifiability (L3-03; `design/audit.py:150-170,412-418`). wired. | no | **P0** (partial confounding is the common real case): L3-03 (2-3 d). |
| A2.6 | M/F/S | Per-sample outlier/QC score so one bad sample is not read as a batch effect. | **no** (L3 coverage matrix row "per sample QC score: none"). | no | P1 (2 d). |
| A2.7 | M/F/S | Within-file drift (signal vs time) with a step/trend flag. | **partial**: exists+wired `diagnostics/drift.py`; event-index not time (A1.4); step statistic unflagged; thresholds not configurable (L3-08, L3-09). | no | P1. |
| A2.8 | M/F/S | Cross-sample drift over run order/date using reference-sample or QC-channel medians (Levey-Jennings-style). | **no**: `SampleRecord` has no run-order field; no QC role (L3-06). | no | **P0 for mass** (see B-M rows), P1 otherwise. |
| A2.9 | M/F/S | Panel/antibody-lot/reagent aliasing check against batch. | **no**: second `panel_id` aborts; lots only as free `technical_factors`, no detector (L3-12). | no | P1 (1 d). |
| A2.10 | M/F/S | Eligibility gating: which corrections the design can support (reference layout, replicates, balance). | **partial**: exists+wired `design/audit.py:43-115`, hand-synchronised table, CytoNorm shows eligible on a layout the adapter refuses (L2-03, L5-02). | no | P1. |

### A3. Corrections (Mod = M/F/S; method availability)

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| A3.1 | M/F/S | A no-correction (identity) control that passes through the identical validate/score path, so every metric has a null. (Good-practice; scib uses controls, OpenProblems lists 5 control methods [R8].) | **yes** exists+wired+tested `adapters/identity.py:63-71` (L2 "what is solid"). | no | none. |
| A3.2 | M/F/S | cyCombine: reference-free, SOM-cluster-wise ComBat on asinh data, optional covariate, no technical replicates needed [R1]. | **partial**. exists+wired `adapters/cycombine.py`, `workflow/scripts/cycombine_correct.R`; routed correctly (L2 "solid"). Defects: protected factors not passed as `covar` (L2-05); `anchor` declared but non-functional (L2-04); single-batch / confounded SOM nodes returned uncorrected yet status SUCCESS and the policy never reads it (L2-02, L2-09); transductive (fits on validation cells) but scored as held-out (L2-07); smoke-level tests only (L2-08). | no (simulated fixtures) | **P0**: L2-02 (completeness in the decision), L2-07 (`fit_scope`), L2-05 (1-2 d each). |
| A3.3 | M/F/S | CytoNorm: reference/control-based, per-FlowSOM-cluster quantile splines fitted on shared controls, applied to all samples [R2]. Needs the *same* control material in every batch. | **partial**. exists+wired `adapters/cytonorm.py`, `cytonorm_correct.R` with strong output checks (L2 "solid"). **Does not require controls to share a `reference_id`** (L2-01, `cytonorm.py:236-269`); partial correction reported as SUCCESS (L2-02); clustering ignores `use_for_clustering` (L2-06); `goal`/`limit` not exposed (L2-13); design limits typed `not_available` not `ineligible`, and the shipped `paired_reference` example cannot fit it (L2-10, `docs/blocked.md:185-186`). | no | **P0**: L2-01 (control identity), L2-02. P1: L2-06, L2-10, L2-13. |
| A3.4 | M/F/S | Fit-once / apply-to-new-batch with a frozen, serialisable normalisation model (CytoNorm model object; CytoNormPy/2.0 port, details unverified [R2b]). | **no** as a user workflow: CytoNorm fit-on-train/apply-on-validate happens inside a run only; cyCombine returns a table, not a model (L2 section "Fit/apply ... gap"). | no | P1 (needed for longitudinal/clinical use; 2-3 d). |
| A3.5 | F (also M/S without controls) | Reference-free per-channel landmark/peak registration: gaussNorm and warping (Hahne et al. 2010; flowStats `gaussNorm`, `warpSet`; `fdaNorm` itself no longer shipped, per secondary source: unverified) [R6]. In OpenProblems benchmark (GaussNorm) [R8]. | **no** (L2 gap 1). | no | P1 (3-5 d). The only reference-free option today is cyCombine; flow studies without a control have nothing else. |
| A3.6 | M | BatchAdjust: per-protein per-batch scaling from anchor/reference replicates in every barcode batch, applied to all samples [R5]. In OpenProblems (no/one/all controls) [R8]. | **no** (L2 gap 4). | no | P1 (2-3 d); closes the documented "one control per batch" limitation. |
| A3.7 | M | CytofRUV: replicate-based RUV (technical replicates across batches) with diagnostic plots [R3]. | **no** (L2 gap 2). Design layer already models replicate/anchor roles, so the data model fits. | no | P1 (4-6 d). Exact replicate requirements: unverified. |
| A3.8 | M/F/S | iMUBAC: Harmony applied to healthy-control cells only, cluster, then batch-specific classifiers transfer labels back to uncorrected patient data (patients never corrected) [R4]. | **no**. | no | P2 (design pattern more than a backend; relevant to over-correction concerns). |
| A3.9 | M/F/S | Embedding-level integration (Harmony [R9], rPCA, scVI, limma `removeBatchEffect`); in OpenProblems [R8]. | **no**; contract = same-shape corrected matrix, cannot carry embeddings (L2 gap 7). | no | P2 (4-6 d contract change + 1-2 d/method). |
| A3.10 | M/F/S | Safe, reproducible execution: pinned upstream commits with licence boundary, parameter accept/reject lists, typed failure states, output validation. | **yes** exists+wired+tested (L2 "solid"; `config/upstream_lock.yaml`). Gaps: non-`BackendError` exceptions crash the rule and leave no `failed` record (L1-03); `event_identity_preserved` tautological (L2-12); no R timeout (L2-14). | no | P1: L1-03 (0.5-1 d). |
| A3.11 | M/F/S | Method plug-in cost (new backend = one module). | **partial**: 4-6 manual edit sites; eligibility not derived from capabilities (L1-07, L2-03). | n/a | P2. |

### A4. Evaluation (post-correction: did batch effect fall, was biology kept, what do we recommend?)

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| A4.1 | M/F/S | Batch-alignment metric: EMD between batches per marker (and per population) before vs after, as in CytoNorm's EMD evaluation on held-out controls [R2] and cyCombine `evaluate_emd` [R1]. | **partial**. exists+wired `metrics/common.py:106-140`, `metrics/deltas.py:56-101` (`emd_relative_reduction`, identical comparison sets: solid). But pooled over cells, no sample/donor uncertainty (L4-05), finite-sample bias floor ignored (L4-08), and computed on **all** frozen cells including training samples and CytoNorm control aliquots, not held-out (L5-01, `policy.py:170-198`). | no | **P0**: L5-01 + L4-05 (held-out scoring; sample-level EMD with donor-block bootstrap, 3 d). |
| A4.2 | M/F/S | Upstream-native cyCombine EMD / MAD scores for each candidate (the community-standard numbers [R1]). | **no (wired wrongly)**. exists `metrics/native.py`; `pipeline.py:399-419` runs it only when candidate is identity, on identical corrected/uncorrected frames, so every real candidate gets none (L4-02, L1-02); identity gets degenerate/NaN. | no | **P0** (L4 "blocker"; 0.5 d fix). |
| A4.3 | M/F/S | OpenProblems cyto-integration metric set (EMD per cell type horizontal/vertical, peak/R2-type, LISI-type etc.) so results are comparable with the community benchmark [R8]. | **no (not wired)** `metrics/openproblems.py:337` `compute_openproblems_metrics`: no caller in src/workflow/scripts (L4-03). Even if called, both EMD rows score the uncorrected baseline (`openproblems.py:360-368`; L4-01). Adapted grouping, `upstream_verified` is a hand-set flag (L4-10). | no | **P0**: L4-01 (0.5 d) before any wiring; wiring itself P1 (3-4 d). |
| A4.4 | M/F/S | Neighbourhood-mixing metrics: iLISI/cLISI [R9], kBET [R10], graph connectivity and silhouette-based batch/label ASW (scib set [R11]). | **partial**: iLISI/cLISI exist `metrics/openproblems.py:907-1011` but not wired (L4-03). kBET, silhouette/ASW, graph connectivity, PCR: **absent** (L4-11). | no | P1 (kBET/ASW 2-3 d; note scib metrics were designed for scRNA-seq: applicability to marker-limited cytometry is an assumption, unverified). |
| A4.5 | M/F/S | Biology-preservation: label-transfer / classifier accuracy on frozen labels, cluster/population-frequency stability (scib bio-conservation family [R11]; iMUBAC uses classifier transfer [R4]). | **no**. `population_frequency_shift` never reads the candidate matrix (identical for all candidates: L4-04; `metrics/common.py:179-206`). | no | **P0**: a package that scores only batch removal can recommend over-correction (6-10 d with A4.6). |
| A4.6 | M/F/S | Donor-level preservation of biological contrasts (effect size / test before vs after, restricted permutation, donor-block bootstrap). | **partial / not wired**: strong design in `metrics/preservation.py:252-428,664-720` (Stage 8, 18/18 tests) but no caller, no contrast source in config (L4-03, L1-01); no baseline-delta, no FDR, bootstrap anti-conservative at n=3-6 (L4-12). | no | **P0**: wire + baseline-delta + FDR (3-4 d + 2 d). |
| A4.7 | M/F/S | Over-correction guardrails (variance collapse, correlation/joint structure, rare-population loss). | **partial**: wired MAD-compression and joint-structure proxies (`common.py:142-251`; `policy.py:43-58`); pooled, magnitude-blind, uncalibrated (L4-06, L4-07); **fail open when the metric is missing** and the regression test cannot fail (L5-03, `policy.py:231-233,265-274`). Rare-population checks absent. | no | **P0**: L5-03 (fail closed, 0.5 d), L4-06 (1.5 d). |
| A4.8 | M/F/S | Replicate-based evaluation: distance between technical replicates / blinded control across batches before vs after (CytoNorm's validation design [R2]; CytoBatchFlagR technical-control logic [R15]). | **no** (L4-11: no replicate-CV / paired technical-replicate distance). | no | **P0 for reference designs**; it is the one evaluation that does not depend on frozen labels. 2-3 d. |
| A4.9 | M/F/S | Uncertainty on every decision metric (donor/sample bootstrap) with intervals reported. | **no** for common metrics (`uncertainty_*` null; L4-05); exists only in unwired Stage 8. Stale text says "not implemented" (`config.py:92-100`). | no | **P0** (covered by A4.1, 3 d). |
| A4.10 | M/F/S | Recommendation / abstention logic with calibrated thresholds. | **partial**: exists+wired `decisions/policy.py`; honest refusal vocabulary and deterministic; built-in policies cannot reach `recommended_under_policy` by construction (L5 verdict 1); all thresholds "uncalibrated engineering default"; refusal is design-wide not contrast-wide (L5-05); in-split control shortage mislabelled (L5-02). | no | **P0** (no calibration against any real or null data; the package currently only abstains). |
| A4.11 | M/F/S | Report: before/after embeddings and density, per-metric intervals, run provenance, redaction for safe export. | **partial**: exists+wired `reporting/html.py`, `redact.py`, JSON exports; before/after density for 3 markers only, no embeddings, detector PNGs not embedded (L3-04); metrics intervals absent. | no | P1. |
| A4.12 | M/F/S | Reproducibility: seeds, hashed inputs, versioned schemas, status record generated from real test reports. | **yes** (L1 verdict 1; `RELEASE.md` section 4; `scripts/generate_status.py:88-90`). Gaps: `schema_version` not enforced on read (L1-05). | n/a | P2. |

## B. Mass cytometry (CyTOF / IMC-adjacent suspension data), modality-specific capabilities

Repo-wide fact for all rows in B, C and S: the schema has no `modality` field, `compensated`/`unmixed` are declared booleans only
(`schemas/contracts.py:124-130`), and nothing in `src/` mentions beads (grep `bead` in `src` returns no hits outside docstrings).

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| B.1 | M | Bead-standard normalisation of signal drift within and across runs (Finck et al. 2013; CATALYST `normCytof`) [R7]; this reduced median-signal fluctuation 4.9x to 1.3x in the paper's one-month series. | **no**. Cannot be a backend under the current contract (needs raw, untransformed data and bead channels; L2 gap 3; `adapters/cytonorm.py:11-13`). No `BEAD`/`QC` marker role (L3-06). | no | **P1** (P0 if the package claims to diagnose mass drift): 3-5 d, requires keeping a raw layer in the dataset contract. At minimum record "bead-normalised yes/no" as required metadata. |
| B.2 | M | Run-order drift diagnostic from bead or reference-sample channels (Levey-Jennings style). | **no** (L3-06; no run-order field). | no | P1 (3-4 d). |
| B.3 | M | Channel crosstalk / isotopic spillover compensation (Chevrier et al. 2018, bead-based, CATALYST `compCytof`; spillR as a newer alternative) [R16]. | **no**; only the declared flag `compensated` and an export note "do not compensate again" (`io/export.py:104-105`). | no | P1 (needed to audit that inputs are compensated: 1 d; implementing: 3-4 d). |
| B.4 | M | Barcoding / debarcoding and carry-over QC (batch = barcode set; ~20 samples per barcode run per BatchAdjust paper [R5]). | **no**. Barcode set as batch is expressible in the design layer; debarcode quality is not checked. Primary-source for debarcoding: unverified. | no | P2. |
| B.5 | M | Anchor / reference replicate designs in every barcode batch, supported with identity of material tracked (BatchAdjust, CytofRUV, CytoNorm all depend on it [R2][R3][R5]). | **partial**: design layer models `reference_role` and `shared_reference_ids` (L2 gap 2; `design/audit.py:105-114,190-192`) but CytoNorm ignores the shared id (L2-01). exists+wired. | no | **P0** (= L2-01). |
| B.6 | M | asinh transform, cofactor ~5 for mass (CATALYST convention [R7]); applied before correction. | **partial**: `preprocess/transform.py:154` exists, not wired (A1.3). No cofactor estimation. | no | P1. |
| B.7 | M | Real CyTOF validation set(s) with technical replicates across batches (e.g. data released with CytofRUV [R3] / BatchAdjust [R5] / iMUBAC [R4]; availability and licence of each: unverified). | **no**. No mass-cytometry data of any kind in the repo or in the audit (`data/DATA_MANIFEST.md`). | **no** | **P0** (single biggest validation gap for the mass modality; 5-10 d including a reference-design analysis). |

## C. Conventional flow cytometry, modality-specific capabilities

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| C.1 | F | Compensation (spillover matrix) applied and verified before batch correction; consistency of spillover across batches/instrument settings. Primary citation: unverified (not searched). | **no**; declared boolean + "do not compensate again" export policy (`io/fcs.py:5,28`). | no | P1 (compensation-consistency diagnostic 2-3 d). |
| C.2 | F | Transform: logicle / biexponential / asinh with per-channel cofactor (OpenProblems accepts arcsinh or logicle [R8]; Logicle primary paper: unverified). | **partial**: asinh only; `logicle` deliberately removed from config because no tested inverse (`config.py:209-214`). | no | P1 (declared-not-faked behaviour is correct; add logicle with an inverse test, 2 d). |
| C.2b | F | Reference-free per-channel landmark normalisation for studies with no control sample. | **no** = A3.5. | no | P1 (as A3.5). |
| C.3 | F | Instrument standardisation tracking: reference-bead MFI / PMT-voltage over time; target values and Levey-Jennings charts (standardisation initiatives, e.g. EuroFlow; primary citation unverified). | **no** (L3-11, L3-06). | no | P1 (shares machinery with B.2). |
| C.4 | F | Flow-rate / margin / dynamic-range QC per file (flowAI [R13], PeacoQC [R12]). | **partial**: PeacoQC wrapper exists, not wired (A1.1); flowAI absent. | no | P1 (wiring B2). |
| C.5 | F | Real conventional-flow validation data with a shared control across batches (e.g. CytoNorm's flow datasets [R2]; availability: unverified). | **no**. | **no** | **P0**. |

## S. Spectral flow cytometry, modality-specific capabilities

Literature caveat: my searches found no peer-reviewed benchmark of batch correction specific to spectral cytometry
(the OpenProblems task and the CytoNormPy preprint [R2b] are modality-agnostic). Spectral-specific "targets" below rest on
a workflow paper [R18] and on vendor/preprint guidance [R17], and are therefore labelled **unverified** as community consensus.

| ID | Mod | target | implemented | real data | gap + priority |
|---|---|---|---|---|---|
| S.1 | S | Unmixing quality diagnostics: unmixing residuals / spread, autofluorescence signature per batch, reference-control quality and matching to stained samples, tandem-dye lot effects (AutoSpectral preprint [R17], workflow paper [R18]; both only partly verified). | **no**; `unmixed` is a declared boolean (`schemas/contracts.py:124-130`; L3-11). | no | P1; **P0** if spectral correction is advertised without it. |
| S.2 | S | Cross-batch panel handling: shared-marker/channel-intersection mode when layouts differ between batches (OpenProblems requires all markers returned [R8]). | **no**; more than one `panel_id` aborts load, "shared-marker mode ... not implemented" (`io/dataset.py:124-130`). Public set: six channel layouts over 18 batches, only batches 13-18 homogeneous, 48 of 72 files excluded (`DATASET.md`). | partial: layouts audited only | P1 (3-5 d); without it 2/3 of the only real data is unusable. |
| S.3 | S | Spectral-aware transform (wide dynamic range, negative values after unmixing, per-channel cofactor) and a check that CytoNorm spline range / cyCombine behave on it. | **no spectral-specific test** (L2 gap 6); transform not wired (A1.3). | no | P1 (1-2 d audit on the 24-file subset). |
| S.4 | S | Acquisition-time-aware drift on spectral files (Time channel present in the public set). | **no**: Time channel never read (L6-01, A1.4). | no (channel exists in DATASET.md audit) | P1 (0.5-1 d). |
| S.5 | S | A first end-to-end real run: raw FCS -> QC -> diagnose (CytoBatchFlagR on technical controls) -> correct -> evaluate, on the 24-file, 6-batch subset. | **no**: data audited and subset frozen (`data/metadata/cytobatchflagr_spectral/technical_control_subset.tsv`); blocked on raw-FCS preprocessing wiring and a marker manifest; "no correction result has been validated" (DATASET.md; `docs/validation_results.md:24-26`). | **no (data audit only)** | **P0** (first real-data evidence for any capability). |
| S.6 | S | Establish whether the 4 controls per batch (`tcu1/tcu2/tcs3/tcs4`) are the same material across batches, i.e. whether CytoNorm's reference assumption holds. | **unknown**: DATASET.md says table `sample_id` is a batch/control label and "does not establish a biological donor ID". Combined with L2-01 (CytoNorm never checks reference identity) the first real CytoNorm run could silently fit unlike material. | no | **P0** (verify from the Zenodo record before running; 0.5 d). |

## D. References (targets). All titles/DOIs from web searches done 2026-10-09 unless flagged

- [R1] Pedersen CB et al. "cyCombine allows for robust integration of single-cell cytometry datasets within and across technologies." Nat Commun 2022;13:1698. doi:10.1038/s41467-022-29383-5. https://pmc.ncbi.nlm.nih.gov/articles/PMC8971492 (verified by search; MDS/UMAP/EMD detail from the repo's own use of upstream, not re-read).
- [R2] Van Gassen S et al. "CytoNorm: A Normalization Algorithm for Cytometry Data." Cytometry A 2020;97(3). doi:10.1002/cyto.a.23904. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7078957/ (verified).
- [R2b] "CytoNormPy enables a fast and scalable removal of batch effects in cytometry datasets." bioRxiv doi:10.1101/2024.07.19.604225 (preprint; says it implements "CytoNorm 2.0"; what 2.0 changes, including any control-free mode, is **unverified**).
- [R3] Trussart M et al. "Removing unwanted variation with CytofRUV to integrate multiple CyTOF datasets." eLife 2020. doi:10.7554/eLife.59630. https://elifesciences.org/articles/59630 (verified title/venue; replicate-based mechanism **unverified** from the search snippets).
- [R4] Ogishi M et al. "Multi-batch cytometry data integration for optimal immunophenotyping" (iMUBAC). J Immunol 2021 (bioRxiv 2020.07.14.202432; https://pmc.ncbi.nlm.nih.gov/articles/PMC7855665; code github.com/casanova-lab/iMUBAC). Mechanism per search summary: Harmony on healthy-control cells only, label transfer to uncorrected data (verified at summary level). Journal volume/DOI unverified.
- [R5] Schuyler RP et al. "Minimizing Batch Effects in Mass Cytometry Data" (BatchAdjust). Front Immunol 2019. doi:10.3389/fimmu.2019.02367. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6803429/ (verified).
- [R6] Hahne F et al. "Per-channel basis normalization methods for flow cytometry data." Cytometry A 2010;77(2):121-131, PMID 19899135; implemented as flowStats `gaussNorm` / `warpSet` (rdrr.io/github/RGLab/flowStats/man/gaussNorm.html). `fdaNorm` availability: **unverified**.
- [R7] Finck R et al. "Normalization of mass cytometry data with bead standards." Cytometry A 2013;83:483-494. https://pmc.ncbi.nlm.nih.gov/articles/PMC3688049 ; CATALYST `normCytof` vignette https://bioconductor.org/packages/release/bioc/vignettes/CATALYST/inst/doc/preprocessing.html (verified).
- [R8] OpenProblems "Cyto batch integration" task (in development; 21 methods incl. CytoNorm, cyCombine, Harmony, Seurat RPCA, ComBat, GaussNorm in no/one/all-control variants; 8 metrics). https://openproblems.bio/benchmarks/cyto_batch_integration ; github.com/openproblems-bio/task_cyto_batch_integration. No finished leaderboard seen.
- [R9] Korsunsky I et al. "Fast, sensitive and accurate integration of single-cell data with Harmony." Nat Methods 2019;16:1289-1296. doi:10.1038/s41592-019-0619-0 (verified; defines LISI, cited from my knowledge of the paper: the search snippet does not mention LISI: **LISI attribution unverified**).
- [R10] Buttner M et al. "A test metric for assessing single-cell RNA-seq batch correction" (kBET). Nat Methods 2019;16:43-49. doi:10.1038/s41592-018-0254-1 (verified).
- [R11] Luecken MD et al. "Benchmarking atlas-level data integration in single-cell genomics" (scIB). Nat Methods 2022;19:41-50. doi:10.1038/s41592-021-01336-8 (verified; 14 metrics, scRNA-seq/ATAC scope).
- [R12] Emmaneel A et al. "PeacoQC: Peak-based selection of high quality cytometry data." Cytometry A 2022;101:325-338. doi:10.1002/cyto.a.24501 (verified; flow and mass).
- [R13] Monaco G et al. "flowAI: automatic and interactive anomaly discerning tools for flow cytometry data." Bioinformatics 2016;32:2473-2480. doi:10.1093/bioinformatics/btw191 (verified).
- [R15] cytoFlagR / CytoBatchFlagR: "A comprehensive framework to objectively assess high-parameter cytometry data for batch effects." bioRxiv doi:10.1101/2025.05.27.656370; a Cytometry A listing (2026) was seen via a secondary source: **journal citation unverified**.
- [R16] Chevrier S et al. "Compensation of signal spillover in suspension and imaging mass cytometry." Cell Syst 2018;6(5):612-620. doi:10.1016/j.cels.2018.02.010 (verified); spillR (Bioconductor) newer alternative.
- [R17] "AutoSpectral improves spectral flow cytometry accuracy through optimised spectral unmixing and autofluorescence-matching at the cellular level." bioRxiv doi:10.1101/2025.10.27.684855 (preprint, not peer reviewed; plus vendor blogs: treat as practice guidance, **unverified as consensus**).
- [R18] "How to prepare spectral flow cytometry datasets for high dimensional data analysis: a practical workflow." Front Immunol 2021 (PMC8640183; DOI believed 10.3389/fimmu.2021.768113: **DOI unverified**); covers QC, cleaning, transformation, batch correction.
- Not searched / left unverified: gating standards (openCyto), Logicle display method, EuroFlow instrument standardisation, debarcoding (Zunder), MMD/EMD-based sample-distance evaluation, crumblr (CLR count models; named from the local skill list, not a literature search).

## E. Top 10 gaps, ranked

Ranking = scientific-validity risk first, then how much each one blocks the others. "Wired" status is from L1-L6.

1. **No real-data run of anything; no mass or conventional-flow data at all** (B.7, C.5, S.5). Every threshold and every "PASSED" gate is on simulation. First step: execute the 24-file spectral subset end to end (blocked on 3 below). **P0**.
2. **Decision metrics are scored on all frozen cells, not held-out, with no sample/donor uncertainty** (A4.1, A4.9; L5-01, L4-05). The headline held-out guarantee is weaker than stated. **P0**.
3. **Raw-FCS preprocessing (gating, asinh, PeacoQC) not wired**, so no real file can enter the pipeline (A1.1-A1.3; L6-02). Gate for gap 1. **P0 enabler / P1 on its face**.
4. **Native cyCombine metrics run only identity-vs-identity, and OpenProblems EMD rows score the baseline instead of the candidate** (A4.2, A4.3; L4-01, L4-02). Mechanical fixes (0.5 d each) that currently make metric families meaningless or invariant. **P0**.
5. **Biology-preservation is not evaluated in any run** (A4.5, A4.6; Stage 8 unwired, `population_frequency_shift` candidate-invariant, guardrails fail open: L4-03, L4-04, L5-03, L5-04). A corrector that merges populations can still pass. **P0**.
6. **CytoNorm does not check that controls are the same material across batches, and the public spectral control identity is itself unestablished** (A3.3, B.5, S.6; L2-01, DATASET.md). **P0**.
7. **Partial / transductive correction is reported as SUCCESS and scored like a full held-out fit** (A3.2, A3.3; L2-02, L2-07). Needs a completeness / `fit_scope` field read by the policy. **P0**.
8. **Diagnostics can mislead: CytoBatchFlagR flags with an unreported ~0.2-0.3 false-flag rate, composition test that cannot fire at 2 samples/batch, binary-only confounding** (A2.3-A2.5; L3-01, L3-02, L3-03). **P0**.
9. **No modality awareness and no modality-specific handling**: no bead normalisation or drift tracking (mass), no compensation/spillover or instrument-standardisation checks (flow), no unmixing/autofluorescence/panel-layout handling (spectral) (B.1-B.3, C.1, C.3, S.1-S.3; L3-11, L3-06, L2 gap 3/6). **P1** (P0 for any claim of spectral/mass support).
10. **Method coverage is 2 of ~9 benchmark families**: no reference-free flow landmark method (gaussNorm/warp), BatchAdjust, CytofRUV, iMUBAC-style or embedding integration; and missing mixing/over-correction metrics (kBET, ASW, label-transfer, replicate-CV) (A3.4-A3.9, A4.4, A4.8; L2 gaps 1,2,4,7; L4-11). **P1**.

Not ranked but cheap: run-order/Time-channel drift (A1.4, L6-01), config-tunable detector thresholds (L3-09), `ineligible` vs `not_available` (L2-10).

## Limits of this matrix

- Static: built from sibling lenses and targeted greps; nothing executed. Where a lens said "high confidence, not reproduced", that carries over.
- Literature column: verified at title/venue/DOI level from web search results (standard mode), not by reading full texts; mechanism statements are marked unverified where the search snippet did not support them. No spectral-specific correction benchmark was found; that absence is itself reported as unverified-absence.
- L7-L10 lenses were not present when this was written (only L1-L6 exist); rows may need revisiting after them.
