---
title: "Vignette: when the honest answer is 'don't correct'"
---

# When the honest answer is "don't correct"

Most batch-correction tools answer one question: *how do I remove the batch effect?*
This vignette runs the case where the correct answer is **you cannot, and here is
why** — and shows what a tool that refuses to say so would have done to your study.

**Every state and reason below is from a real executed run** (version 0.1.0, seed
473). The scenario is the shipped `confounded` example: two batches, no measured
reference, and — the fatal property — the biological group of interest is perfectly
confounded with the batch. Every case sample was measured in batch A; every control in
batch B.

## 1. Run it

```bash
batchdx-cyto init my-project --example confounded
batchdx-cyto benchmark --config my-project/config/config.yaml \
    --outdir my-project/results/run1 --cores 4
```

## 2. What the decision said

`recommended_candidate_id` is **null**:

| candidate | state | reason |
|---|---|---|
| `uncorrected` | `below_design_floor` | 2 batch(es) present; the declared minimum is 3; no measured reference sample is shared across any pair of the batches under comparison |
| `cycombine_default` | `insufficient_evidence` | at least one protected factor is perfectly confounded with the batch factor |

## 3. Why this is the correct answer

With group and batch perfectly confounded, the group term and the batch term span the
same subspace. Any "batch effect" the diagnostics find is *also* the biological effect
of interest, and no algorithm — however clever — can tell them apart, because the data
contains no information that distinguishes them. A correction that "succeeds" here has
removed the biology, or kept it and called it batch; there is no third option.

Three design properties work together to make the refusal automatic and legible:

1. **The design audit classifies the contrast before any correction runs.** The
   contingency of batch × group is computed from the sample manifest, the comparison
   is marked `non_identifiable`, and the cyCombine backend is marked ineligible with
   `reason_code: batch_not_identifiable`. This entry is **non-overridable**: the
   exploratory override (`design.allow_exploratory_override`) exists for eligibility
   failures that are a matter of judgement, and arithmetic unavailability is not one.
2. **The design floor stops the run from pretending.** Two batches is below the
   declared `design.minimum_batches` (default 3), so every candidate is
   `below_design_floor` and the system produces descriptive output only. The floor is
   a configuration value with a declared origin, echoed into `design_audit.json` —
   not a hidden constant.
3. **Adding the factor as a covariate does not help, and the report says so.** A
   common folk remedy for confounding is to throw the confounded factor into the
   model as a covariate. With *perfect* confounding there is nothing to estimate: the
   design matrix is rank-deficient. The audit states this explicitly instead of
   letting the user discover it in a reviewer comment.

## 4. What a leaderboard would have done

A tool that scores corrections and picks the best number would have "corrected" this
data — the batch effect and the biology are the same signal, so removing one is
trivially easy and looks spectacular on every alignment metric. The published study
would then report that the group effect vanished after batch correction, which is
exactly backwards: the group effect *was* the batch effect, and no one can say which.

The package's answer is the least dramatic possible output — a null recommendation,
two named states, two reasons — and it is the only output that is scientifically
honest. **Abstention is a first-class outcome, not a failure path.**

## 5. What you should do with such a design

The package will not tell you how to fix your study, but the audit's structure points
at what is missing: batches that contain *both* biological groups (so group and batch
are estimable separately), and a measured technical reference shared across batches
(so the technical component has ground truth). The
[paired-reference example](quickstart.md) is the minimal design that has both.

## Next

* [A first diagnosis, end to end](quickstart.md) — the same pipeline on a design that
  *can* support correction.
* [The user guide: confounding](../usage.md) — the reference requirements table and
  the design floor.
