---
title: Home
---

# batchdx-cyto

Batch-effect **diagnosis**, **correction** and **evaluation** for mass cytometry,
conventional flow cytometry and spectral flow cytometry, orchestrated by Snakemake.

**Scope: synthetic fixtures only.** This package ships no real cytometry dataset or
real-data example. Its scientific checks use simulated data; its thresholds are not
calibrated for a real cohort. See [validation limits](validation_results.md).

Built for cross-cohort immune studies, where biological differences and technical
batches are often correlated and where the honest answer is frequently *"this contrast
cannot be corrected"*.

> **Status: 0.1.0 — research preview.** The Python core, the design audit, the
> diagnostics, the metrics (including seven verified against pinned upstream
> components), the decision layer, the report and the Snakemake workflow run end to
> end and are tested, including real Snakemake executions and the R-backed cyCombine
> and CytoNorm paths against their pinned upstream checkouts. Three things you must
> know before using it:
>
> 1. **Every decision threshold is an uncalibrated engineering default.** The report
>    says so next to every recommendation. Do not cite the thresholds as validated
>    cutoffs.
> 2. **No real-data validation yet.** All evidence comes from synthetic fixtures of
>    250–2,000 cells per sample. Production CyTOF acquisitions are orders of
>    magnitude larger; expect the memory and runtime profile to be untested.
> 3. **CytoNorm and PeacoQC execute for real** at their pinned versions, and CI
>    builds the full R environment (cyCombine, CytoNorm, PeacoQC, flowVS) from the
>    pinned local checkouts, so the R-backed tests run there too. They remain the
>    least exercised backends, and their GPL licences mean they are called, never
>    vendored.

---

## What it does

1. **Diagnoses** batch-associated differences and asks whether a *technical*
   interpretation is actually supported by the design.
2. **Runs eligible correction methods independently**, always including an identity
   candidate, because "no correction" is a valid result.
3. **Evaluates** technical improvement and, separately, what the correction did to
   the data's structure — including whether the biological contrast of interest
   survived.
4. **Reports**: one offline HTML file, structured JSON/CSV records, and corrected
   datasets with an explicit scale.

Three things it deliberately refuses to do:

* **It will not turn association into attribution.** A batch contrast is classified as
  `reference_supported`, `biologically_stratified_conditional`, `descriptive_only` or
  `non_identifiable`, and perfect confounding blocks automatic recommendations for the
  affected contrast. Adding the confounded factor as a covariate does not fix it, and
  the audit says so.
* **It will not invent evidence.** An unavailable dependency, an invalid comparison, an
  absent control or a failed method gets an explicit status and a reason. A missing
  metric is never a zero, and a near-zero denominator is `undefined`, not an epsilon.
* **It will not recommend under uncalibrated thresholds.** A candidate that clears
  every gate under the shipped policies reaches `eligible_for_review`, never
  `recommended_under_policy`. Automatic approval requires a policy whose thresholds
  carry a declared calibration — which means doing the calibration.

## Quick start

```bash
pip install -e '.[workflow]'
batchdx-cyto init my-project --example paired_reference
batchdx-cyto benchmark --config my-project/config/config.yaml \
    --outdir my-project/results/run1 --cores 4
```

Then open `my-project/results/run1/report/report.html` — one offline file, no server.
The [quickstart vignette](vignettes/quickstart.md) walks through a full run and reads
the decision with you.

## Where to go next

| I want to… | Go to |
|---|---|
| run my first diagnosis | [Quickstart vignette](vignettes/quickstart.md) |
| understand what the package refuses to do | [When the honest answer is "don't correct"](vignettes/confounded.md) |
| prepare my own FCS files, manifests and config | [User guide](usage.md) |
| benchmark corrections on a real spectral panel | [Real-world vignette](vignettes/spectral_controls.md) |
| see what is validated and what is not | [Validation results](validation_results.md) |
| know why the code is shaped this way | [Design decisions](decisions/0001-internal-data-contract.md) |
| find a function or contract | [API reference](api/index.md) |

## Citation

`CITATION.cff` in the repository root carries the citation metadata.
