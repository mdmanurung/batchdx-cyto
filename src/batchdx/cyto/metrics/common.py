"""Common-partition metrics.

Every metric here is computed on the *same* frozen evaluation cells, labels, marker
subset and comparison pairs for every candidate, so that two candidates' numbers are
comparable.  Nothing in this module re-derives a partition from a candidate's own
output.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from batchdx.cyto.diagnostics.composition import frequency_table
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.metrics.records import metric_record as _record
from batchdx.cyto.metrics.registry import (
    EMD_WITHIN_POPULATION,
    JOINT_STRUCTURE,
    MAD_SIGNED_CHANGE,
    POPULATION_FREQUENCY_SHIFT,
)
from batchdx.cyto.metrics.uncertainty import attach_emd_uncertainty
from batchdx.cyto.sampling.partitions import UNGATED
from batchdx.cyto.schemas.contracts import (
    DesignAudit,
    MetricRecord,
    MetricResult,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import ComparisonSupport, MetricStatus

#: Minimum cells in a (population, batch) group before a distance is reported.
MIN_GROUP_CELLS = 30
#: Minimum independent samples per (population, batch) group.  Cells are not
#: independent replicates; a group supplied by a single sample has no replication
#: basis, whatever its cell count.
MIN_SAMPLES_PER_GROUP = 2
#: A MAD change below this is treated as no change rather than as compression.
MAD_TOLERANCE = 1e-9


def _grouping(dataset: Dataset, sampling: SamplingManifest):
    """Index the frozen evaluation cells by (batch, population) and by sample."""

    obs = dataset.obs
    keep = obs[CELL_ID].isin(sampling.labels_by_cell.keys())
    obs = obs.loc[keep].copy()
    obs["_label"] = obs[CELL_ID].map(sampling.labels_by_cell)
    batch_of = dataset.sample_batches()
    obs["_batch"] = obs[SAMPLE_ID].map(batch_of).astype(str)
    positions = {cid: i for i, cid in enumerate(dataset.obs[CELL_ID])}
    obs["_row"] = obs[CELL_ID].map(positions)
    return obs


def _comparisons(design: DesignAudit) -> dict[tuple[str, str], str]:
    out: dict[tuple[str, str], str] = {}
    for c in design.comparisons:
        if c.population is not None:
            continue
        if c.support is ComparisonSupport.NON_IDENTIFIABLE:
            continue
        out[(c.batch_a, c.batch_b)] = c.comparison_id
    return out


def _donors(sample_ids, donors_of: dict[str, str | None]) -> int | None:
    """Distinct donors behind ``sample_ids``, or ``None`` when there are none.

    Every family in this module reports it, so the report's ``n_donors`` column means
    the same thing on a dispersion row as on a distance row.  A run in which no sample
    declares a donor reports ``None`` and the page shows an empty cell; it never
    reports zero, because "no donor was declared" and "zero donors were available" are
    different facts and only the first is true.
    """

    found = {donors_of.get(str(s)) for s in sample_ids}
    found.discard(None)
    return len(found) or None


def compute_common_metrics(
    dataset: Dataset,
    candidate_matrix: np.ndarray,
    sampling: SamplingManifest,
    design: DesignAudit,
    *,
    run_id: str,
    candidate_id: str,
) -> MetricResult:
    """All common-namespace metrics for one candidate, on the frozen partition."""

    obs = _grouping(dataset, sampling)
    markers = sampling.evaluation_markers
    if not markers:
        return MetricResult(
            run_id=run_id, candidate_id=candidate_id, metric_group="common",
            status="ineligible",
            reason="no marker is flagged use_for_evaluation in the marker manifest",
        )

    feature_idx = dataset.feature_index(markers)
    baseline = dataset.expression()
    records: list[MetricRecord] = []
    # The residual 'ungated' pool is a heterogeneous leftover, not a biological
    # population: it is reported descriptively (counts, composition) but is excluded
    # from the quantitative population metrics, where it could otherwise dominate the
    # worst-marker summaries that drive the no-effect gate.
    populations = sorted(
        pop for pop in set(sampling.labels_by_cell.values()) if pop != UNGATED
    )
    pairs = _comparisons(design)

    donors_of = {sid: rec.donor_id for sid, rec in dataset.samples.items()}

    # -- distribution alignment within frozen populations ----------------------------
    for pop in populations:
        pop_rows = obs.loc[obs["_label"] == pop]
        for (batch_a, batch_b), comparison_id in sorted(pairs.items()):
            rows_a = pop_rows.loc[pop_rows["_batch"] == batch_a]
            rows_b = pop_rows.loc[pop_rows["_batch"] == batch_b]
            n_a, n_b = len(rows_a), len(rows_b)
            samples = sorted(set(rows_a[SAMPLE_ID]) | set(rows_b[SAMPLE_ID]))
            n_samples_a = rows_a[SAMPLE_ID].nunique()
            n_samples_b = rows_b[SAMPLE_ID].nunique()
            donors = sorted({donors_of.get(s) for s in samples if donors_of.get(s)})
            coverage = float(min(n_a, n_b)) / MIN_GROUP_CELLS if MIN_GROUP_CELLS else 1.0

            for marker in markers:
                common = dict(
                    run_id=run_id, candidate_id=candidate_id, marker=marker, population=pop,
                    comparison_id=comparison_id, n_samples=len(samples),
                    n_donors=len(donors) or None, n_cells=n_a + n_b,
                    coverage=min(coverage, 1.0), aggregation_weight=1.0,
                )
                if n_a < MIN_GROUP_CELLS or n_b < MIN_GROUP_CELLS:
                    records.append(_record(
                        EMD_WITHIN_POPULATION, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=(
                            f"population {pop!r} has {n_a} and {n_b} evaluation cells in "
                            f"{batch_a} and {batch_b}; {MIN_GROUP_CELLS} are required in each"
                        ),
                        **common,
                    ))
                    continue
                if n_samples_a < MIN_SAMPLES_PER_GROUP or n_samples_b < MIN_SAMPLES_PER_GROUP:
                    records.append(_record(
                        EMD_WITHIN_POPULATION, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=(
                            f"population {pop!r} has {n_samples_a} and {n_samples_b} "
                            f"independent sample(s) in {batch_a} and {batch_b}; cells are "
                            f"not independent replicates and {MIN_SAMPLES_PER_GROUP} "
                            "samples per batch are required"
                        ),
                        **common,
                    ))
                    continue
                col = markers.index(marker)
                x = candidate_matrix[rows_a["_row"].to_numpy(), feature_idx[col]]
                y = candidate_matrix[rows_b["_row"].to_numpy(), feature_idx[col]]
                records.append(_record(
                    EMD_WITHIN_POPULATION, sampling, value=float(
                        stats.wasserstein_distance(x, y)
                    ), status=MetricStatus.OK, **common,
                ))

    # -- uncertainty on the EMD records (Stage 8) -------------------------------------
    # Donor-block bootstrap intervals, or an explicit statement of why none could
    # be computed.  Uncertainty is never silently omitted.
    attach_emd_uncertainty(records, dataset, sampling, candidate_matrix,
                           seed=sampling.seed)

    # -- dispersion -------------------------------------------------------------------
    n_compressed = n_total = 0
    for pop in populations:
        pop_rows = obs.loc[obs["_label"] == pop]
        for batch in sorted(set(obs["_batch"])):
            rows = pop_rows.loc[pop_rows["_batch"] == batch]
            idx = rows["_row"].to_numpy()
            for col, marker in enumerate(markers):
                group_samples = sorted(set(rows[SAMPLE_ID]))
                common = dict(
                    run_id=run_id, candidate_id=candidate_id, marker=marker, population=pop,
                    comparison_id=f"batch:{batch}", n_cells=len(idx),
                    n_samples=len(group_samples),
                    n_donors=_donors(group_samples, donors_of),
                    aggregation_weight=1.0,
                )
                if len(idx) < MIN_GROUP_CELLS:
                    records.append(_record(
                        MAD_SIGNED_CHANGE, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=f"{len(idx)} cells, {MIN_GROUP_CELLS} required", **common,
                    ))
                    continue
                if rows[SAMPLE_ID].nunique() < MIN_SAMPLES_PER_GROUP:
                    records.append(_record(
                        MAD_SIGNED_CHANGE, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=(
                            f"{rows[SAMPLE_ID].nunique()} independent sample(s); cells are "
                            f"not independent replicates and {MIN_SAMPLES_PER_GROUP} are "
                            "required"
                        ),
                        **common,
                    ))
                    continue
                f = feature_idx[col]
                mad_base = float(stats.median_abs_deviation(baseline[idx, f], scale="normal"))
                mad_cand = float(
                    stats.median_abs_deviation(candidate_matrix[idx, f], scale="normal")
                )
                delta = mad_cand - mad_base
                n_total += 1
                if delta < -MAD_TOLERANCE:
                    n_compressed += 1
                records.append(_record(
                    MAD_SIGNED_CHANGE, sampling, value=delta, status=MetricStatus.OK,
                    reason=f"MAD baseline {mad_base:.4g}, candidate {mad_cand:.4g}",
                    **common,
                ))

    # -- composition ------------------------------------------------------------------
    freq = frequency_table(dataset, sampling.labels_by_cell)
    batch_of = dataset.sample_batches()
    for pop in populations:
        for (batch_a, batch_b), comparison_id in sorted(pairs.items()):
            samples_a = sorted(s for s in freq if batch_of.get(s) == batch_a)
            samples_b = sorted(s for s in freq if batch_of.get(s) == batch_b)
            fa = [freq[s][pop] for s in samples_a]
            fb = [freq[s][pop] for s in samples_b]
            common = dict(
                run_id=run_id, candidate_id=candidate_id, population=pop,
                comparison_id=comparison_id, n_samples=len(fa) + len(fb),
                n_donors=_donors(samples_a + samples_b, donors_of),
                aggregation_weight=1.0,
            )
            if not fa or not fb:
                records.append(_record(
                    POPULATION_FREQUENCY_SHIFT, sampling,
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=f"no samples in one of {batch_a}, {batch_b}", **common,
                ))
                continue
            if len(fa) < MIN_SAMPLES_PER_GROUP or len(fb) < MIN_SAMPLES_PER_GROUP:
                records.append(_record(
                    POPULATION_FREQUENCY_SHIFT, sampling,
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=(
                        f"{len(fa)} and {len(fb)} independent sample(s) in {batch_a} and "
                        f"{batch_b}; {MIN_SAMPLES_PER_GROUP} per batch are required, "
                        "because the sample is the independent unit"
                    ),
                    **common,
                ))
                continue
            records.append(_record(
                POPULATION_FREQUENCY_SHIFT, sampling,
                value=float(abs(np.mean(fa) - np.mean(fb))), status=MetricStatus.OK,
                reason="samples are weighted equally, so a high-yield sample cannot dominate",
                **common,
            ))

    # -- joint structure --------------------------------------------------------------
    if len(markers) >= 2:
        for pop in populations:
            pop_rows = obs.loc[obs["_label"] == pop]
            for batch in sorted(set(obs["_batch"])):
                rows = pop_rows.loc[pop_rows["_batch"] == batch]
                idx = rows["_row"].to_numpy()
                group_samples = sorted(set(rows[SAMPLE_ID]))
                common = dict(
                    run_id=run_id, candidate_id=candidate_id, population=pop,
                    comparison_id=f"batch:{batch}", n_cells=len(idx),
                    n_samples=len(group_samples),
                    n_donors=_donors(group_samples, donors_of),
                    aggregation_weight=1.0,
                )
                if len(idx) < MIN_GROUP_CELLS:
                    records.append(_record(
                        JOINT_STRUCTURE, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=f"{len(idx)} cells, {MIN_GROUP_CELLS} required", **common,
                    ))
                    continue
                if rows[SAMPLE_ID].nunique() < MIN_SAMPLES_PER_GROUP:
                    records.append(_record(
                        JOINT_STRUCTURE, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                        reason=(
                            f"{rows[SAMPLE_ID].nunique()} independent sample(s); cells are "
                            f"not independent replicates and {MIN_SAMPLES_PER_GROUP} are "
                            "required"
                        ),
                        **common,
                    ))
                    continue
                c_base = _correlation(baseline[np.ix_(idx, feature_idx)])
                c_cand = _correlation(candidate_matrix[np.ix_(idx, feature_idx)])
                if c_base is None or c_cand is None:
                    records.append(_record(
                        JOINT_STRUCTURE, sampling, status=MetricStatus.UNDEFINED,
                        reason="at least one evaluation marker is constant in this group",
                        **common,
                    ))
                    continue
                tri = np.triu_indices(len(markers), k=1)
                records.append(_record(
                    JOINT_STRUCTURE, sampling,
                    value=float(np.mean(np.abs(c_cand[tri] - c_base[tri]))),
                    status=MetricStatus.OK,
                    reason="mean absolute change in within-population marker correlations",
                    **common,
                ))
    else:
        records.append(_record(
            JOINT_STRUCTURE, sampling, run_id=run_id, candidate_id=candidate_id,
            status=MetricStatus.INELIGIBLE,
            reason="joint structure needs at least two evaluation markers",
        ))

    warnings: list[str] = []
    from batchdx.cyto.metrics.registry import MAD_COMPRESSION

    if n_total:
        records.append(_record(
            MAD_COMPRESSION, sampling, run_id=run_id, candidate_id=candidate_id,
            value=float(n_compressed) / n_total, status=MetricStatus.OK,
            n_cells=None, aggregation_weight=1.0,
            reason=f"{n_compressed} of {n_total} (population, batch, marker) groups lost spread",
        ))
    else:
        # A missing metric is never silently absent: with no qualifying group the
        # compression guardrail is unavailable, and that is reported as such.
        records.append(_record(
            MAD_COMPRESSION, sampling, run_id=run_id, candidate_id=candidate_id,
            value=None, status=MetricStatus.INSUFFICIENT_DATA,
            n_cells=None, aggregation_weight=1.0,
            reason="no (population, batch, marker) group met the cell and sample minima, "
                   "so a compression fraction cannot be computed",
        ))

    n_failed = sum(1 for r in records if r.status is MetricStatus.FAILED)
    status = "ok" if not n_failed else "partial"
    return MetricResult(
        run_id=run_id, candidate_id=candidate_id, metric_group="common",
        records=records, status=status, warnings=warnings,
    )


def _correlation(block: np.ndarray) -> np.ndarray | None:
    sd = block.std(axis=0)
    if np.any(sd <= 0):
        return None
    return np.corrcoef(block, rowvar=False)
