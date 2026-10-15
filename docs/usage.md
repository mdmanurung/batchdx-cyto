# Using batchdx-cyto

## 1. Installation

```bash
pip install -e '.[dev]'     # from a checkout
batchdx-cyto info           # prints the backends, detectors and metrics this build ships
```

`batchdx-cyto info --json` is the authoritative answer to "can this build do X". It
reports, per backend, whether fit-and-apply is supported *and whether that support has
been verified*, which parameters are accepted, and which are refused and why.

The R backends need an R environment on top. See
`containers/Dockerfile.r-backends` — **not built or tested in this build** — and
`docs/blocked.md`.

## 2. Preparing input

The package reads **preprocessed** data and refuses to guess what has been done to it.

### The expression object

Either an H5AD with the values in a **named layer** (declared as
`input.expression_layer`; `X` is never used implicitly), or a directory of FCS files
with one file per sample named in the sample manifest.

### `samples.tsv`

One row per measured sample. Only `sample_id` is required; everything else enables
something specific.

| Column | Enables |
|---|---|
| `sample_id` | required |
| `batch` | every batch contrast |
| `donor_id` | donor-matched metrics, and `n_donors` on every metric record. It does **not** buy donor-level uncertainty intervals; none are computed in this build (Stage 8 is not implemented) |
| `biospecimen_id`, `aliquot_id` | recorded on `SampleRecord` and kept distinct in every emitted manifest, so the levels cannot be silently merged. Nothing in this build reads them |
| `group`, `stimulation`, `timepoint` | protected factors and stratified contrasts |
| `reference_id` + `reference_role` | reference-supported evidence, CytoNorm, CytoBatchFlagR |

Any other column is kept as a technical factor and can be named as
`design.correction_batch` or `design.nested_within`.

`reference_role` is one of `none`, `technical_replicate`, `anchor`,
`biological_pseudo_control`. **Declaring a role without a `reference_id` is an error**,
because reference identity is never inferred from file names. Only
`technical_replicate` and `anchor` count as measured technical evidence;
`biological_pseudo_control` never does, whatever the configuration says.

The same donor at two timepoints, or under two stimulations, is **not** a technical
replicate. The package will not treat it as one.

### `markers.tsv`

| Column | Meaning |
|---|---|
| `feature_id` | required, unique |
| `channel` | required for FCS input; must map one-to-one to `feature_id` |
| `marker`, `panel_id` | labels; more than one `panel_id` is rejected — run panels separately |
| `role` | `lineage`, `functional`, `technical`, `viability`, `barcode`, `other` |
| `use_for_clustering`, `use_for_correction`, `use_for_evaluation` | three **independent** masks |

The three masks are independent on purpose. Clustering a correction backend does
internally is a different thing from the evaluation markers, and neither is
automatically the correction scope.

### `processing.json`

```json
{
  "compensated": true, "unmixed": false,
  "compensation_source": "acquisition-time spillover matrix",
  "qc_history": ["PeacoQC v1.2 upstream of this workflow"],
  "transformation": "asinh", "transformation_parameters": {"cofactor": 5.0},
  "scale": "asinh_cofactor5", "units": "arcsinh units", "corrected": false
}
```

`scale` is mandatory. `corrected: true` is refused: the input is the uncorrected
baseline, and chained correction is not supported by default. There is **no universal
cofactor** — declare it per panel, and do not apply spectral-flow defaults to CyTOF.
The community conventions from the cyCombine paper are cofactor **5 for CyTOF**,
**150 for conventional flow** and **6000 for spectral flow**; they are conventions,
not defaults, so declare the value that matches the instrument and panel. A `scale` of
the form `asinh_cofactor<N>` must agree with `transformation_parameters.cofactor`; a
disagreement is rejected at load time.

## 3. Configuring a run

`config/README.md` lists three shipped configurations; `batchdx-cyto init` materialises
any of them with its data. Validation is strict and unknown keys are rejected:

```bash
batchdx-cyto validate --config config/config.yaml
```

Validation refuses three different things, and the difference matters. An **unknown** key
is a typo. A **removed** key is one this build used to accept and then never read, so
setting it produced a run that looked configured and was not; `preprocessing.time_channel`
is the only one in that category now, and `validate` names it and says what to use
instead. A **refused value** is a key that stays but cannot take every value:
`evaluation.sampling.unit` is `sample_id` and nothing else,
`evaluation.representation.scaling` is `none`, and `report.offline` is `true`, because
this build has exactly one sampler, one coordinate system and one offline renderer. A key
is never accepted at a value the package does not honour.

`preprocessing.gating` and `preprocessing.peacoqc` were in the removed category until
2026-09-30 and are now live keys, read as **paths**: to a gate table and to a PeacoQC
settings file, both in the shape the pinned reference implementation's own configuration
files declare. They are read by `batchdx.cyto.preprocess`, which is **not** yet called by
a run — preprocessing operates on raw events before a baseline exists, and wiring it into
a pipeline stage, a rule and the report is still to be done
([`docs/blocked.md`](blocked.md) B-009). `preprocessing.transformation` no longer offers
`"logicle"`: this build implements no logicle, has no tested inverse for one, and the
reference implementation contains none either.

### Reference requirements, by capability

| Capability | Needs |
|---|---|
| descriptive diagnosis | nothing beyond a batch column |
| drift diagnostic | an acquisition axis: the FCS Time channel (ingested as `acquisition_time`) when present, otherwise the event index |
| composition diagnostic | frozen labels and ≥ 2 samples per batch |
| `reference_supported` evidence | a measured reference id shared by both batches of the pair |
| cyCombine cluster detector | 3–30 batches (upstream restriction, quoted in the eligibility manifest) |
| CytoBatchFlagR | a control sample in **every** batch, and ≥ 4 batches |
| CytoNorm (strongest tier) | a measured reference in every batch |
| CytoNorm (control-free tier) | nothing beyond ≥ 2 batches; weaker evidence (CytoNorm 2.0 `QuantileNorm.train goal="mean"` or a manual goal distribution) |

### Confounding

If every batch contains exactly one level of a protected factor, the batch term and the
factor term span the same subspace and neither is estimable. The audit marks the
comparison `non_identifiable`, marks the `cycombine` backend ineligible with
`reason_code: batch_not_identifiable`, and marks that entry **non-overridable**. An
exploratory override exists (`design.allow_exploratory_override`) for eligibility
failures that are a matter of judgement; it does not apply to a contrast that is
arithmetically unavailable.

Adding the confounded factor as a covariate does not help. The report says so.

### The design floor

`design.minimum_batches` (default 3) and `design.require_shared_reference` (default
true). Both are configuration values with a declared origin, echoed into
`design_audit.json` as `design_floor_origin` — not hidden constants. Below the floor
every candidate is `below_design_floor` and the system produces descriptive output only.

## 4. Running

```bash
batchdx-cyto diagnose  --config config/config.yaml --outdir results/run1   # no correction
batchdx-cyto benchmark --config config/config.yaml --outdir results/run1   # all candidates
batchdx-cyto correct   --config config/config.yaml --outdir results/run1 --candidate cycombine_default
batchdx-cyto report    --config config/config.yaml --outdir results/run1   # re-render only
```

`correct` requires `--candidate`. There is no default and no automatic choice.

Every command delegates to Snakemake; `README.md` shows the equivalent native
invocation. Pass `--cores` to control parallelism. A Snakemake profile
(`workflow/profiles/local/config.yaml`, or the **untested**
`workflow/profiles/slurm/config.yaml`) applies to the native invocation: the command-line
wrapper has no `--profile` flag, because a profile rewrites scheduler and resource
behaviour and a wrapper that hid that would be a wrapper you could not see through.

## 5. Reading the output

```
results/run1/
  config.resolved.json        resolved configuration + analytical fingerprint + versions
  manifests/                  dataset manifest, sampling manifest, exported JSON Schemas
  design/design_audit.json    contingency, rank, comparisons, eligibility
  diagnostics/<detector>/     detector_result.json, plots, native tables
  candidates/<id>/            correction_result.json, corrected.h5ad
  metrics/<id>/               metric_result.json
  decisions/                  decision_result.json
  status/run_status.json      one entry per planned unit, with its reason
  report/report.html          the offline report
```

### Interpreting metric records

Every record carries its namespace, direction, scale, status, coverage, aggregation
weight, the label tier and confidence, and the five representation settings (space,
marker mask, scaling, distance, k). Three things follow:

* **Namespace matters.** `native` numbers come from upstream running its own
  partitioning (including cyCombine's moving uncorrected baseline). `common` numbers
  come from the one frozen partition. They are not interchangeable and a lower `common`
  EMD is not a better `native` score.
* **Status matters more than value.** `insufficient_data`, `undefined` and `ineligible`
  are distinct, all carry a reason, and none is a zero.
* **`label_tier: clustering` is a caveat, not a footnote.** Those labels came from the
  uncorrected data whose batch structure is under test.

### Interpreting the decision

`not_needed`, `eligible_for_review`, `recommended_under_policy`, `rejected_distortion`,
`insufficient_evidence`, `non_identifiable`, `below_design_floor`, `execution_failed`.

The shipped policies use **uncalibrated** thresholds, and the report says so. A
candidate that clears every gate under an uncalibrated policy reaches
`eligible_for_review`, never `recommended_under_policy`. That is deliberate: these are
engineering defaults, not validated universal cutoffs, and they must not be cited as
such. To get automatic approval you must supply a policy whose `threshold_origins` state
a real calibration — which means doing the calibration.

Retaining the uncorrected data is always an available outcome.

### Calibration, if you want it

There is no calibration in this build. Calibrating means: choose a design like yours,
establish ground truth (measured technical replicates are the only strong form), sweep
`min_emd_reduction`, `max_mad_compression_fraction` and `max_joint_structure_change`,
and record what error rate each choice buys you. Then declare that in
`threshold_origins`. Until then, review.

## 6. Failure recovery

* A scientific skip is **not** a workflow failure. An ineligible or unavailable method
  produces a result file with a status and a reason, the run succeeds, and the report
  shows it under "Failures and limitations".
* A genuine runtime failure keeps failure status; nothing is suppressed with `|| true`
  or a blanket `except`.
* Snakemake reruns what is missing. Deleting `metrics/<id>/metric_result.json` and
  `decisions/decision_result.json` forces both to be recomputed, because `decide`
  declares the metric files as real inputs rather than depending on a sentinel.
* Changing a scientifically material setting changes the analytical fingerprint.
  Changing the report title does not.

## 7. Memory and scale

`correct` holds all relevant cells jointly for cyCombine: ComBat is fitted per SOM node
across the whole dataset, so **correcting files independently is not equivalent** and
chunking does not make it out-of-core. Size `resources.mem_mb` for the full matrix.
Metrics and diagnostics are memory-bounded per group and parallelise across candidates.

No production-scale performance claim is made. The largest thing exercised in this build
is a synthetic fixture of a few tens of thousands of cells.

## 8. Privacy

`report.redact_identifiers` is on by default: sample, donor and reference identifiers in
the report become stable pseudonyms. The pseudonym mapping goes to
`report/identifier_map.json`, deliberately a separate file, so a report can be circulated
without it.

Identifiers are redacted. **Paths are redacted too, on a stated rule.** A path inside the
project becomes `<run>/…`; a path outside it becomes `<redacted>/<file name>`, so the
reader learns which file was read and not where it lives. This covers the resolved
configuration the page prints, and the free-form payload each detector records for
itself — its parameters, requirements, table and plot names, and warnings. A detector's
parameters can carry one entry per evaluation cell (the pinned `cycombine_express` API
groups by a column literally named `sample`, so a per-cell sample vector crosses the
boundary), and a list longer than eight entries is replaced by its length and its first
entries, with the page saying so.

Two limits, both stated on the page rather than here. A path *inside* a captured log line
is not rewritten, because rewriting a log would falsify it — read such a warning as the
machine's own text. And the dataset manifest's provenance is not rendered at all: only
`dataset_id` reaches the page.

## 9. Extending

* **A correction backend** implements the `CorrectionBackend` protocol
  (`adapters/base.py`): declare `BackendCapabilities` honestly — especially
  `supports_fit_apply` *and* `fit_apply_verified` — then `register_backend`.
  Unsupported parameters must go in `rejected_parameters` with a reason.
* **A detector** implements the `Detector` protocol (`diagnostics/base.py`) and returns
  a `DetectorResult` for a scientific skip rather than raising.
* **A metric** registers a `MetricSpec` (`metrics/registry.py`). Use the `common`
  namespace unless you have verified your implementation against a pinned upstream one;
  `scripts/audit_upstream.py` fails the build if anything claims `compatibility` without
  `upstream_verified`.
* Adding an R call means adding it to `adapter_calls` in `config/upstream_lock.yaml`,
  or the Stage 0 test will fail — by design, so no adapter can quietly grow a call the
  audit never checks. That gate is `test_r_audit_runs_and_finds_every_adapter_referenced_function`
  in `tests/integration/test_stage0_audit.py`, and it **skips** when no R interpreter is
  on `PATH`; on a machine without R, adding an adapter call is not caught. Run the R
  environment before relying on it.

## 10. Known limits

Read `docs/implementation_status.md` (generated) and `docs/blocked.md` before relying on
anything. In short: the R backends are not executable in the environment this was built
in; biological preservation (Stage 8) is not implemented, and no candidate here has been
shown to preserve biology; the OpenProblems metrics are deferred; there is no real-data
example; and no clinical or unseen-batch generalisation claim is made.
