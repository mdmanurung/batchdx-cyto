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
demonstrates behaviour on real data. **This package ships no real cytometry data and no
real-data example, as a stated scope decision** — see item 4 below and
`examples/README.md`, where the same statement sits on the face of the documentation
rather than only in a blocked-items file.

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

1. **Biological preservation is demonstrated on synthetic fixtures only.** Stage 8 is
   now implemented: prespecified contrasts, paired and unpaired sample-level models,
   donor-block bootstrap, restricted permutation, and an uncertainty interval on every
   emitted record. **But** the fixtures are deterministic and synthetic, so this
   establishes that the machinery detects, refuses and quantifies correctly — it is not a
   calibration for any real cohort, and no threshold in it should be cited as one. The
   single most important limit on this package is that it has never seen real cytometry
   data (item 4).
2. **Both correction backends and the CytoBatchFlagR detector execute.**
   cyCombine 0.3.0 and CytoNorm 2.0.12 run at their pinned commits; CytoBatchFlagR's
   detector executes against a four-batch fixture and distinguishes its injected shift
   from the corresponding unshifted design. Its 1.5-IQR flag remains uncalibrated for
   multiple tests, so a flag is not a validated false-positive-controlled decision.
   CytoNorm only normalises batches with a fitted control; the result reports the
   samples it left unchanged. See `docs/blocked.md` B-005 and B-006.
3. **Six upstream metric components are verified; three metric names are not.** Six of
   the OpenProblems components at commit `37151479` are registered in the
   `compatibility` namespace with `upstream_verified=True`, each proven by a test that
   imports the **pinned helper file** and compares outputs on a shared fixture at the
   tolerances tabulated below: `emd_mean_ct_horiz`, `emd_mean_ct_vert`,
   `ratio_consistent_peaks`, `average_batch_r2_ct`, `iLisi` and `cLisi`. Extraction did
   not need the Nextflow runtime. What is **not** verified is the case-assembly loop in
   the components' `script.py` files, which read `.h5ad` resources the pinned checkout
   does not carry; those loops are documented transcriptions with the classifier tested
   directly. The remaining three names are in `DEFERRED_METRICS` with a currently-true
   obstacle, and `openproblems.cms` is recorded as genuinely absent at the pinned
   commit. The audit script fails the build if anything claims `compatibility` without
   verification.
4. **No real data, and this is a scope decision rather than a gap.** Behaviour on heavy
   tails, multimodality, marker dropout, instrument artefacts and genuine spillover is
   untested, and every scientific result here rests on `batchdx.cyto.simulate`. No
   permissibly licensed real dataset was obtainable under this environment's egress
   policy, and none was written from memory (`docs/blocked.md` B-008). **This package
   ships no real cytometry data and no real-data example, deliberately** — no gate depends
   on real data, so the honest move is to state the boundary where a reader will meet it
   rather than leave it as a not-done item that never converges. Read any number in this
   file as a property of the code on synthetic input, not as a calibration for any cohort.
   The same statement is on the face of `examples/README.md`,
   `docs/compatibility_matrix.md` and `README.md`.
5. **Multi-seed, but only over the claims that can be checked across seeds.**
   `tests/scientific/test_multiseed.py` runs **25 tests over five seeds**
   (`473, 7, 101, 20250930, 65537`), and it holds at every one of them: a correct offset
   removal reduces alignment error on both shifted markers; a metric recomputed on the
   same data is bit-identical; the same seed reproduces both the data and the frozen
   evaluation draw; different seeds produce different data *and* different draws; a
   confounded design is `non_identifiable` and recommends nothing; and a calibrated
   policy *can* still recommend on the reference-supported design, so the refusal above is
   the design and not a permanently closed gate. The seeds also caught a real defect: the
   shipped examples wrote a configuration declaring `seed: 473` regardless of the seed
   their data was generated with, so `materialise_example(..., seed=7)` shipped data and a
   run record that disagreed.

   **What this is not.** It is not a five-seed re-run of the scientific suite. It cannot
   be one: the existing scientific tests take their data from the `frozen` fixture in
   `tests/conftest.py`, which calls `materialise_example` without a seed, so their data is
   473 whatever a parametrisation says. The cost is not the obstacle either —
   `tests/scientific` takes 1m47s, so five passes would be under ten minutes. Making it a
   genuine suite-wide sweep needs a seed-aware conftest, which is a change to a shared file
   and therefore a central change, not something the sweep can do for itself. Until that
   lands, **Stage 7 metrics, Stage 8 preservation and Stage 9 decisions remain
   single-seed**, and this file says so rather than letting a passing 25-test file imply
   otherwise.
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
| pinned `emd/helper.py:bin_array` | exact (`np.array_equal`) | pure `np.histogram`; any difference is a porting error, not sampling noise |
| pinned `emd/helper.py:compute_emd`, per marker | < 1e-9 relative | a weighted `wasserstein_distance` over a fixed 2000-bin grid; only float summation order can differ |
| pinned `emd/helper.py` horizontal and vertical means | < 1e-9 relative | the same distances, pooled over more markers in a different order |
| pinned `ratio_consistent_peaks/helper.py` density and peak count | < 1e-12 absolute, exact counts | pure `gaussian_kde` and `find_peaks`; deterministic, so the bar is the port being correct |
| pinned `average_batch_r2/helper.py`, per marker and aggregate | < 1e-9 relative | measured agreement 2e-16. Only correct once one dummy is fitted **per batch**; a single 0/1 contrast is a different model and a test asserts the two differ by more than 0.05 |
| pinned `lisi/helper.py` per-cell iLisi and cLisi | exact (`rtol=0, atol=0`) | `pynndescent` is deterministic at the pinned seed, so per-cell parity is achievable rather than merely statistical; the normalisation is compared at 1e-12 relative |
| CytoNorm FCS round trip | to the float32 storage bound | FCS stores float32, so a round trip is not bit-exact and the bound is measured (`FCS_STORAGE_DTYPE`, `round_trip_tolerance`) rather than assumed |
| CytoNorm correction magnitude, injected marker vs others | affected marker > 0.4 mean absolute; others ≤ 0.05 | only CD4 carries an injected offset on the fixture, so a correction that moved every marker would be a red flag rather than a success |

These are predeclared in the tests. They were not loosened to obtain a favourable
result; where a test initially failed, the cause was a fixture that did not contain the
structure being tested for (no co-expression to destroy; a collapse that did not
actually align), and the **fixture** was corrected, not the tolerance.
