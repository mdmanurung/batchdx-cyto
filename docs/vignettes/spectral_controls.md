---
title: "Vignette: benchmarking corrections on a real spectral-cytometry control panel"
---

# Benchmarking corrections on a real spectral-cytometry control panel

The [quickstart](quickstart.md) showed the machinery on synthetic data. This vignette
is the **real-world workflow**: taking a public spectral-flow cytometry panel through
the same pipeline — from raw FCS files to a per-candidate decision — with your own
manifests, your own panel declaration and the memory discipline the R backends need.

!!! warning "What is validated here"

    The workflow, the manifests and the configuration shown below were built and
    schema-validated against the real dataset described, and the design audit facts
    quoted are from that audited sample table. The *decision thresholds* remain
    uncalibrated engineering defaults, and the package ships no real-data validation —
    treat the run's output as a structured diagnosis, not a validated recommendation.

## 1. The dataset

We use the spectral cytometry control panel published with CytoBatchFlagR
([Zenodo record 15388817](https://zenodo.org/records/15388817), DOI
`10.5281/zenodo.15388817`, CC-BY-4.0): unstained and singly-stained control
acquisitions spread over six instrument batches. The subset used here:

* **24 FCS files**: batches 13–18 × control samples `{tcs3, tcs4, tcu1, tcu2}`
* **~6.15 M cells** total; layout of 33 channels: 7 scatter/autofluorescence
  channels, 25 antigen channels, and Time
* **25 declared markers**, including CD45 and the functional markers
  GZMB, GNLY, PFN, PD1, CD39, FOXP3 and TIGIT

Download the record's FCS files into `raw/` in your project directory.

## 2. The sample manifest

One row per file. The control identity is declared explicitly — the package never
infers a reference from a file name:

```tsv
sample_id	batch	donor_id	reference_id	reference_role
B13_tcs3	13	control	tcs3	anchor
B13_tcs4	13	control	tcs4	anchor
B13_tcu1	13	control	tcu1	anchor
B13_tcu2	13	control	tcu2	anchor
B14_tcs3	14	control	tcs3	anchor
...
```

Design properties this encodes, and why each matters:

* **`batch`** — the correction factor. Six levels.
* **`donor_id = control`** for every file. These are technical controls, not
  independent biological donors; declaring them as one donor identity stops the
  metrics from treating tubes as replicates of anything biological.
* **`reference_id` + `reference_role: anchor`** — the same control measured in every
  batch is the measured technical evidence that makes a batch contrast
  `reference_supported`. With this column, all **15 pairwise batch comparisons** of
  the six batches are reference-supported, the design matrix is full rank, and no
  protected factor is confounded with batch — the audit verifies all three.

## 3. The marker manifest

```tsv
feature_id	channel	marker	panel_id	role	use_for_clustering	use_for_correction	use_for_evaluation
...
CD3	Antigen_CD3	CD3	panelA	lineage	1	1	1
GZMB	Antigen_GZMB	GZMB	panelA	functional	0	1	1
Time	Time	Time	panelA	technical	0	0	0
```

The three masks are independent on purpose: clustering scope, correction scope and
evaluation scope are different questions. The Time channel is declared `technical`
and excluded from all three masks. All 25 markers carry the same `panel_id` — the
package runs one panel per project, and (since the panel-id refusal) a marker that
loses its `panel_id` under a declared transformation is an error, not a silent
raw-scale leak into your baseline.

## 4. Processing declaration

`processing.json` declares what has already been done upstream, and the package
refuses to guess:

```json
{
  "compensated": true,
  "unmixed": true,
  "compensation_source": "acquisition-time spectral unmixing",
  "qc_history": [],
  "transformation": "asinh",
  "transformation_parameters": {"cofactor": 150},
  "scale": "asinh_cofactor150",
  "units": "arcsinh units",
  "corrected": false
}
```

**Cofactor 150** is the community convention for conventional/spectral flow from the
cyCombine paper — a convention, not a default; declare the value that matches your
instrument. The `scale` string must agree with the declared cofactor or the file is
rejected at load time. `corrected: true` is refused: this file is the uncorrected
baseline, and chained correction is not supported.

Alternatively, leave the FCS files raw and declare the preprocessing in the config —
the pipeline stage applies per-panel asinh and a PeacoQC pass (declared via
`preprocessing.peacoqc`), refuses a sample whose QC pass removes more than half its
events, and freezes the survivors into the baseline with per-event provenance.

## 5. The run configuration

The material keys of `config/config.yaml`:

```yaml
project:
  seed: 473
input:
  path: "data"
design:
  correction_batch: "batch"
  reference_identity: "reference_id"
  reference_role: "reference_role"
  minimum_batches: 3
evaluation:
  clustering_k: 8
  max_cells_per_sample: 2000
  min_cells_per_population: 50
candidates:
  - uncorrected
  - cycombine_default
  - cytonorm_reference
```

`max_cells_per_sample: 2000` caps the per-sample subsample the frozen partition and
the metrics see — with 24 samples that is a ~48,000-cell evaluation matrix, while the
correction backends see the full ~6 M-cell matrix.

## 6. Running it

```bash
batchdx-cyto benchmark --config config/config.yaml \
    --outdir results/run1 --cores 4
```

or, natively, with a profile:

```bash
snakemake --snakefile /path/to/batchdx-cyto/workflow/Snakefile \
    --profile /path/to/batchdx-cyto/workflow/profiles/local \
    --cores 4 \
    --config config=config/config.yaml outdir=results/run1 \
             base_dir=. mode=benchmark
```

### Memory: size for the full matrix

`correct` holds all relevant cells jointly for cyCombine: ComBat is fitted per SOM
node **across the whole dataset**, so correcting files independently is not
equivalent and chunking does not make it out-of-core. For this 24-file, 6 M-cell
panel, budget **64 GB** of RAM for the correction rule; the metrics and diagnostics
stages are memory-bounded per group and parallelise across candidates. On a cluster,
use the `slurm` profile and give the correction rule the memory; on a laptop, this
panel does not fit comfortably.

## 7. Reading the results

The output layout is the same as in the quickstart. The comparison that matters is
per candidate, across the same frozen partition:

* **Technical improvement** — Σ within-population EMD on the frozen common partition,
  relative to the identity candidate, plus the native cyCombine numbers (namespace
  `native` — not interchangeable with `common`).
* **Preservation** — `mad_compression_fraction` (how much within-group spread
  survived), `joint_structure_distance` (did the neighbourhood geometry survive), and
  `biological_contrast_change` for every prespecified contrast against its predeclared
  tolerance. On a controls-only panel there may be no biological contrast to
  declare — which is itself informative: the decision layer then has no biology gate
  to clear, and the report says so.
* **Composition** — per-batch population frequencies on the frozen labels, with the
  `label_tier: clustering` caveat: those labels came from the uncorrected data whose
  batch structure is under test.
* **CytoBatchFlagR flags** — the acquisition-drift and composition detectors run per
  batch; with a control in every batch and ≥ 4 batches, the eligibility requirements
  are met. The flags' false-positive rate is uncalibrated; read them as a structured
  screen, not a verdict.

The decision layer then assigns each candidate a state with its reason. On this
panel, with all 15 comparisons reference-supported, a candidate that clears the
guardrails reaches `eligible_for_review` — and stops there, because the shipped
thresholds are uncalibrated. The honest reading of the report is: *these are the
measured improvements and the measured costs, this is what the design can support,
and the recommendation is yours to make after review.*

## 8. Scaling the pattern to your study

Everything in this vignette transfers directly:

1. Declare your sample sheet with batch, donor identity, and reference roles — never
   let the tool infer them.
2. Declare your panel, your masks, and the processing that already happened, with the
   cofactor your instrument actually implies.
3. Run `validate` first; read its refusals as design review, not bureaucracy.
4. Benchmark, then read the decision reasons — including for the candidates that
   failed, because a failed method with a logged reason is evidence too.
5. If you need an automatic recommendation, calibrate: choose a design like yours,
   establish ground truth from measured technical replicates, sweep
   `min_emd_reduction`, `max_mad_compression_fraction` and
   `max_joint_structure_change`, record the error rates, and declare them in
   `threshold_origins`. Until then, review.

## Next

* [A first diagnosis, end to end](quickstart.md) — the same pipeline on synthetic
  data, with every number reproducible.
* [The user guide](../usage.md) — every manifest column, config key and refusal.
