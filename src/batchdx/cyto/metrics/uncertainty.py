"""Donor-block bootstrap confidence intervals and restricted permutation tests.

Cells are not independent replicates: they come from donors, and donors are the units
that can be resampled or exchanged.  Everything here is therefore stratified by donor
block, and the permutation test never exchanges batch labels across strata that carry
different protected-factor levels, because such an exchange would manufacture exactly
the confounding the design audit warns about.

Both procedures are deterministic given the project seed: the RNG is seeded once per
call and consumed in a fixed order.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.schemas.contracts import MetricRecord, SamplingManifest
from batchdx.cyto.schemas.enums import MetricStatus

#: Number of donor-block resamples per confidence interval.  A package default,
#: declared here rather than hidden in a call site; the interval is a percentile
#: interval and its width is not a calibrated coverage statement.
BOOTSTRAP_RESAMPLES = 200
#: Declared origin of the bootstrap interval, written into every record reason.
BOOTSTRAP_ORIGIN = (
    f"package default: percentile bootstrap over donor blocks, {BOOTSTRAP_RESAMPLES} "
    "resamples; uncalibrated as a coverage statement"
)
#: Number of restricted permutations per test.
PERMUTATION_RESAMPLES = 1000
#: Donor blocks smaller than this cannot be resampled: with fewer than two donors a
#: bootstrap resample would always draw the same donor and the interval would be a
#: zero-width lie about precision.
MIN_DONORS_PER_GROUP = 2
#: The uncertainty kind written on every record that carries a bootstrap interval.
BOOTSTRAP_KIND = "donor_block_bootstrap_ci"


def donor_block_bootstrap_emd(
    x: np.ndarray,
    y: np.ndarray,
    donors_x: np.ndarray,
    donors_y: np.ndarray,
    *,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int,
    level: float = 0.95,
) -> tuple[float, float] | None:
    """Percentile CI for the two-sample EMD, resampling donors within each group.

    ``x`` and ``y`` are the two batches' cell values on one marker; ``donors_x`` and
    ``donors_y`` give each cell's donor.  Donors are resampled with replacement within
    their own batch, never across batches.  Returns ``None`` when either group has
    fewer than :data:`MIN_DONORS_PER_GROUP` donors.
    """
    ux = np.unique(donors_x)
    uy = np.unique(donors_y)
    if len(ux) < MIN_DONORS_PER_GROUP or len(uy) < MIN_DONORS_PER_GROUP:
        return None
    rng = np.random.default_rng(seed)
    idx_x = {d: np.flatnonzero(donors_x == d) for d in ux}
    idx_y = {d: np.flatnonzero(donors_y == d) for d in uy}
    drawn = np.empty(n_resamples, dtype=float)
    for b in range(n_resamples):
        pick_x = rng.choice(len(ux), size=len(ux), replace=True)
        pick_y = rng.choice(len(uy), size=len(uy), replace=True)
        cells_x = np.concatenate([idx_x[ux[i]] for i in pick_x])
        cells_y = np.concatenate([idx_y[uy[i]] for i in pick_y])
        drawn[b] = stats.wasserstein_distance(x[cells_x], y[cells_y])
    alpha = (1.0 - level) / 2.0
    lo, hi = np.quantile(drawn, [alpha, 1.0 - alpha])
    return float(lo), float(hi)


@dataclass(frozen=True)
class PermutationResult:
    """Outcome of a restricted permutation test, including why it may not run."""

    observed: float
    p_value: float | None
    n_permutations: int
    n_exchangeable_strata: int
    reason: str | None = None


def restricted_permutation_test(
    values: np.ndarray,
    batch_labels: np.ndarray,
    strata: np.ndarray,
    *,
    n_permutations: int = PERMUTATION_RESAMPLES,
    seed: int,
) -> PermutationResult:
    """Two-sided permutation test for a batch difference in means, within strata.

    Batch labels are permuted only inside a stratum (a donor block at fixed
    protected-factor levels).  A stratum contributes exchanges only when it contains
    both batch labels; every other stratum is held fixed.  The test therefore never
    swaps labels across donor blocks and never across different protected levels —
    when no stratum holds both labels there is no valid exchange, and the p-value is
    reported as unavailable rather than as 1.0.
    """
    values = np.asarray(values, dtype=float)
    labels = np.asarray(batch_labels, dtype=object)
    strata = np.asarray(strata, dtype=object)
    levels = sorted(set(labels.tolist()))
    if len(levels) != 2:
        return PermutationResult(
            observed=float("nan"), p_value=None, n_permutations=0,
            n_exchangeable_strata=0,
            reason=f"permutation test needs exactly two batch levels, got {levels}",
        )
    a, b = levels
    mask_a = labels == a
    if not mask_a.any() or not (~mask_a).any():
        return PermutationResult(
            observed=float("nan"), p_value=None, n_permutations=0,
            n_exchangeable_strata=0,
            reason="one batch level has no values; nothing is exchangeable",
        )
    observed = float(values[mask_a].mean() - values[~mask_a].mean())

    idx_by_key = {k: np.flatnonzero(strata == k) for k in np.unique(strata)}
    exchangeable = [
        idx for idx in idx_by_key.values()
        if set(labels[idx].tolist()) == {a, b}
    ]
    if not exchangeable:
        return PermutationResult(
            observed=observed, p_value=None, n_permutations=0, n_exchangeable_strata=0,
            reason=(
                "no stratum contains both batch labels, so no exchange is possible "
                "within donor blocks at fixed protected-factor levels; a permutation "
                "test would either do nothing or permute across different protected "
                "levels, and neither is valid"
            ),
        )

    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(n_permutations):
        perm = labels.copy()
        for idx in exchangeable:
            perm[idx] = rng.permutation(labels[idx])
        stat = values[perm == a].mean() - values[perm != a].mean()
        if abs(stat) >= abs(observed):
            extreme += 1
    p_value = (1.0 + extreme) / (1.0 + n_permutations)
    return PermutationResult(
        observed=observed, p_value=float(p_value), n_permutations=n_permutations,
        n_exchangeable_strata=len(exchangeable),
    )


def attach_emd_uncertainty(
    records: list[MetricRecord],
    dataset: Dataset,
    sampling: SamplingManifest,
    candidate_matrix: np.ndarray,
    *,
    seed: int,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
) -> None:
    """Attach donor-block bootstrap CIs to ``emd_within_population`` records, in place.

    A record whose comparison has at least two donors per batch gains
    ``uncertainty_low``/``uncertainty_high`` and ``uncertainty_kind =
    'donor_block_bootstrap_ci'``.  A record whose comparison does not keeps its value
    but states in the reason that no interval could be computed and why.  Uncertainty
    is never silently omitted and never replaced by a zero-width interval.
    """
    obs = dataset.obs
    keep = obs[CELL_ID].isin(sampling.labels_by_cell.keys())
    obs = obs.loc[keep].copy()
    obs["_label"] = obs[CELL_ID].map(sampling.labels_by_cell)
    batch_of = dataset.sample_batches()
    obs["_batch"] = obs[SAMPLE_ID].map(batch_of).astype(str)
    positions = {cid: i for i, cid in enumerate(dataset.obs[CELL_ID])}
    obs["_row"] = obs[CELL_ID].map(positions)
    donor_of = {sid: rec.donor_id for sid, rec in dataset.samples.items()}
    obs["_donor"] = obs[SAMPLE_ID].map(donor_of).astype(str)

    markers = sampling.evaluation_markers
    if not markers:
        return
    feature_idx = dataset.feature_index(markers)
    col_of = {m: i for i, m in enumerate(markers)}

    for record in records:
        if record.metric_name != "emd_within_population":
            continue
        if record.status is not MetricStatus.OK:
            continue
        if record.population is None or record.marker is None:
            continue
        if record.comparison_id is None or "__vs__" not in record.comparison_id:
            continue
        batch_a, batch_b = record.comparison_id.split("__vs__", 1)
        rows = obs.loc[obs["_label"] == record.population]
        ra = rows.loc[rows["_batch"] == batch_a]
        rb = rows.loc[rows["_batch"] == batch_b]
        n_donors_a = ra["_donor"].nunique()
        n_donors_b = rb["_donor"].nunique()
        if n_donors_a < MIN_DONORS_PER_GROUP or n_donors_b < MIN_DONORS_PER_GROUP:
            record.reason = (
                (record.reason + "; " if record.reason else "")
                + f"no donor-block bootstrap CI: this comparison has {n_donors_a} "
                f"donor(s) in {batch_a} and {n_donors_b} in {batch_b}, and at least "
                f"{MIN_DONORS_PER_GROUP} donors per batch are required to resample"
            )
            continue
        col = feature_idx[col_of[record.marker]]
        ci = donor_block_bootstrap_emd(
            candidate_matrix[ra["_row"].to_numpy(), col],
            candidate_matrix[rb["_row"].to_numpy(), col],
            ra["_donor"].to_numpy(),
            rb["_donor"].to_numpy(),
            n_resamples=n_resamples, seed=seed,
        )
        if ci is None:  # pragma: no cover - guarded by the donor count above
            continue
        record.uncertainty_low = ci[0]
        record.uncertainty_high = ci[1]
        record.uncertainty_kind = BOOTSTRAP_KIND
        record.reason = (
            (record.reason + "; " if record.reason else "")
            + f"donor-block bootstrap CI [{ci[0]:.4g}, {ci[1]:.4g}] "
            f"({BOOTSTRAP_ORIGIN})"
        )
