# Validation results

What the scientific suite establishes, and — more importantly — what it does not.

`docs/implementation_status.md` is generated from an executed test report and is the
authority on which gates passed. This file explains what those tests *mean*.

## The fixtures

All of it runs on deterministic synthetic data from `batchdx.cyto.simulate`. Five
populations with declared marker means, a shared per-cell latent factor within each
population (so markers genuinely co-vary, which is what makes a joint-structure metric
measurable at all), and separately controllable perturbations:

| Knob | Injects |
|---|---|
| `batch_shift` | a technical offset per batch and marker |
| `batch_spread` | a technical dispersion change |
| `group_shift` | a genuine biological difference |
| `drift` | a within-file trend across acquisition order |
| per-sample `frequencies` | a composition difference |
| `within_population_correlation` | co-expression structure to preserve or destroy |

**These are fixtures, not a scientific claim.** No real cytometry dataset was available
under this environment's egress policy (`docs/blocked.md` B-008), so nothing here
demonstrates behaviour on real data.

## Validation-only backends

`tests/fixtures/backends.py` registers three controls that are **not** user correction
methods and are deliberately not registered by the package:

* `fixture_offset_removal` — the *correct* answer for the paired-reference fixture:
  subtract each batch's offset as measured on its reference aliquot.
* `fixture_collapse` — replaces every value with the population's median across all
  batches. Alignment becomes perfect; all spread is destroyed.
* `fixture_joint_scramble` — permutes each marker independently within each
  (population, batch) group. Every marginal is preserved **exactly**; co-expression is
  destroyed.
* `fixture_failure` — always raises, to exercise the failure path.

## What the tests establish

### Metrics respond to what they claim to measure

* Identical distributions give near-zero EMD; a known translation moves it.
* A pure scale change shows up in `mad_signed_change` as a **negative** value. The
  cyCombine MAD score takes an absolute difference and cannot tell restoration from
  compression, which is why the signed metric exists under a different name.
* A constant marker gives `undefined`, never zero.
* A rare population gives `insufficient_data` with the required cell count in the
  reason.
* A near-zero denominator gives `undefined`; no epsilon is applied.
* `median_abs_deviation(scale="normal")` reproduces the 1.4826 consistency constant, so
  the dispersion metric is on the same footing as R's `stats::mad`.

### Composition is separable from intensity

Dropping 80 % of one population's cells from one batch's evaluation set moves the
composition metric by more than 0.05 while every within-population EMD for the other
populations moves by **less than 1e-9**. A frequency shift cannot masquerade as an
intensity artefact.

### The joint-structure metric earns its place

`fixture_joint_scramble` leaves every EMD and every MAD unchanged to within 1e-9, and
the joint-structure metric still detects it. Without that row, a corruption that
preserves all marginals while destroying co-expression would be invisible — which is
precisely the Stage 8 requirement.

### Guardrails reject destructive corrections

`fixture_collapse` achieves better alignment than the baseline *and* is rejected
(`rejected_distortion`), because the compression guardrail fires. A system that rewarded
alignment alone would recommend it.

### The design audit refuses unsupportable claims

Perfect confounding yields `non_identifiable`, marks the cyCombine backend ineligible
with `reason_code: batch_not_identifiable`, and marks that entry **non-overridable**. A
purely biological shift present in every batch does not inflate batch-pair EMD.

### Upstream constraints are derived, not asserted

The four-batch minimum for CytoBatchFlagR is proven arithmetically: over 2000 random
draws of three values, the 1.5 IQR Tukey fence is never exceeded, and a four-value
counterexample is shown. The 3–30 batch range for cyCombine's cluster detector is read
from the pinned source and quoted in the eligibility reason.

### Abstention works

Uncalibrated thresholds produce `eligible_for_review`, never
`recommended_under_policy`. A policy whose `threshold_origins` state a real calibration
*can* recommend — tested — which shows the gate is a gate and not a permanent refusal.
Removing a candidate's own failing guardrail metric does not promote it.

## What these tests do NOT establish

Stated plainly, because the absence of a test is not evidence of correctness:

1. **No biological preservation has been demonstrated.** Stage 8 is not implemented.
   There are no prespecified biological contrasts, no paired or repeated-measure models,
   no donor-block bootstrap, and no uncertainty quantification. The joint-structure
   metric tests within-population co-expression, which is not a biological contrast.
2. **No upstream correction method has been executed.** cyCombine, CytoNorm and
   CytoBatchFlagR could not be installed (`docs/blocked.md` B-001 to B-006). Every
   correction tested here came from the identity backend or a validation fixture. The
   "outputs match a direct upstream call within a documented tolerance" clauses of
   Stages 5 and 6 are **not** satisfied.
3. **No metric has been checked against an upstream implementation.** Nothing in this
   build is registered in the `compatibility` namespace, and the audit script fails the
   build if anything ever claims to be without verification.
4. **No real data.** Behaviour on heavy tails, multimodality, marker dropout, instrument
   artefacts and genuine spillover is untested.
5. **Single seed.** The scientific suite runs at one seed per scenario. Multi-seed
   stochastic checks are not implemented.
6. **No cluster execution.** The SLURM profile has never run on a cluster and says so.
   The container has never been built.
7. **No calibration.** Every recommendation threshold is an uncalibrated engineering
   default. They are not validated cutoffs and must not be cited as such.
8. **No performance claim.** The largest thing exercised is a few tens of thousands of
   synthetic cells.

## Tolerances used

| Check | Tolerance | Why |
|---|---|---|
| marginal preservation under scrambling | 1e-9 | the permutation is exact; anything larger would be a bug |
| excluded-marker invariance | 1e-9 (`INVARIANCE_ATOL`) | backends run in float64; larger is a real modification |
| identical-distribution EMD | < 1e-3 | sampling noise on ~1e3 cells per group |
| offset removal on affected markers | after < 0.5 × before | the offset is 0.6–0.8 on a spread of 0.35 |
| collapse compression fraction | > 0.9 | collapse destroys spread in essentially every group |
| joint-structure change under scrambling | > 0.1 | `rho² = 0.36` co-expression is removed |

These are predeclared in the tests. They were not loosened to obtain a favourable
result; where a test initially failed, the cause was a fixture that did not contain the
structure being tested for (no co-expression to destroy; a collapse that did not
actually align), and the **fixture** was corrected, not the tolerance.
