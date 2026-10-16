---
title: "Vignette: a first diagnosis, end to end"
---

# A first diagnosis, end to end

This vignette runs the package once on its shipped synthetic example and reads every
artifact it produces. **Every number below is from a real executed run** (version
0.1.0, seed 473, 1,500 cells per sample); if you rerun it, you should reproduce them.

The scenario: three batches (B1–B3), a measured reference aliquot in every batch,
a balanced case/control biology sitting on CD69, and a technical batch offset on CD4
and HLA-DR. The biology is real and worth keeping; the batch offset is technical and
worth removing. The question the run answers is: *can any correction method remove the
offset without destroying the biology — and can the package tell?*

## 1. Materialise the example project

```bash
pip install -e '.[workflow]'          # core + Snakemake
batchdx-cyto init my-project --example paired_reference
```

The project contains the input data (an H5AD with the values in a named layer), the
three manifests the package reads — `samples.tsv`, `markers.tsv`, `processing.json` —
and a resolved `config/config.yaml`. `batchdx-cyto info` prints what this build can
do; `batchdx-cyto validate --config config/config.yaml` refuses a misconfigured run
before any compute happens.

## 2. Run the benchmark

```bash
batchdx-cyto benchmark --config config/config.yaml \
    --outdir results/run1 --cores 4
```

`benchmark` runs every configured candidate — here `uncorrected` (the identity),
`cycombine_default` and `cytonorm_reference` — plus the design audit, the diagnostics,
the metrics and the decision layer, as one Snakemake DAG. The run takes a few minutes
on four cores. (`diagnose` runs the same pipeline with the correction candidates
switched off.)

## 3. What the design audit said

`results/run1/design/design_audit.json`:

* 3 batches, 3 population-free batch comparisons.
* **All 3 comparisons are `reference_supported`** — both batches of every pair share a
  measured reference aliquot. This is the strongest evidence class the package
  recognises: a technical interpretation of the batch contrast is attributable.

## 4. What the metrics said

Total within-population EMD over the frozen common partition (lower is better):

| candidate | Σ EMD |
|---|---|
| uncorrected | 27.55 |
| cytonorm_reference | 5.55 |

CytoNorm removed roughly **80%** of the batch-associated distributional difference. On
a leaderboard, that is the winner. The package does not run a leaderboard.

## 5. What the decision said

`results/run1/decisions/decision_result.json` — `recommended_candidate_id` is **null**:

| candidate | state | reason |
|---|---|---|
| `uncorrected` | `eligible_for_review` | the identity candidate is always retained as a valid outcome |
| `cytonorm_reference` | `rejected_distortion` | `mad_compression_fraction 0.517 exceeds 0.25` |
| `cycombine_default` | `execution_failed` | cyCombine returned non-finite values; they are not silently repaired |

Read these one at a time:

* **CytoNorm's 80% EMD reduction is real, and it is not enough.** The
  `mad_compression_fraction` guardrail asks how much of each marker's within-group
  spread survived the correction. CytoNorm compressed 51.7% of it — more than the
  declared 25% ceiling — which means the correction did not just move the batches
  together, it flattened the biology's spread. The candidate is rejected *by name and
  by number*, and the report shows both the improvement and the rejection.
* **cyCombine's failure is a refusal, not a crash.** The R backend produced non-finite
  values; the adapter refused to hand them downstream and recorded the R log tail in
  `correction_result.json`. The run itself still succeeds — a failed method is a
  *result with a status*, never a silent zero and never a suppressed error.
* **The identity candidate stays on the table.** "Do nothing" is always a valid
  outcome, and here it is the only one left standing.

This is the package's core behaviour in one run: **it will not trade a technical
improvement for an unmeasured biological cost, and it will not recommend under
uncalibrated thresholds.** Even a candidate that had cleared every guardrail would
have reached `eligible_for_review`, not `recommended_under_policy` — the shipped
thresholds are engineering defaults, and the report says so next to the recommendation.

## 6. The report

```text
results/run1/report/report.html
```

One offline HTML file: the resolved configuration (paths redacted on a stated rule),
the design audit, every detector and metric record with its status and coverage, the
per-candidate decisions with reasons, and a "Failures and limitations" section where
the cyCombine failure appears with its log. Sample and donor identifiers are
pseudonymised by default; the mapping goes to a separate file so the report can
circulate without it.

## 7. Where the numbers live

```text
results/run1/
  config.resolved.json        resolved configuration + analytical fingerprint
  manifests/                  dataset manifest, sampling manifest, JSON Schemas
  design/design_audit.json    contingency, rank, comparisons, eligibility
  diagnostics/<detector>/     detector_result.json, plots, native tables
  candidates/<id>/            correction_result.json, corrected.h5ad
  metrics/<id>/               metric_result.json
  decisions/                  decision_result.json
  status/run_status.json      one entry per planned unit, with its reason
  report/report.html          the offline report
```

Every metric record carries its namespace (`native` vs `common` — not interchangeable),
its status (`ok`, `insufficient_data`, `undefined`, `ineligible` — a missing metric is
never a zero), its coverage, and the five representation settings that make it
reproducible.

## Next

* [When the honest answer is "don't correct"](confounded.md) — what the package does
  when the design itself makes correction impossible.
* [The real-world workflow](spectral_controls.md) — the same pipeline on a
  spectral-cytometry control panel, with your own manifests.
* [The user guide](../usage.md) — every manifest column, every config key, and the
  three different things validation refuses.
