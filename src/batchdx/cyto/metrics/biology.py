"""Prespecified biological contrasts, before and after correction (Stage 8).

A biological contrast is declared in the configuration: two levels of a sample-level
factor, one marker, and a preservation tolerance with a declared origin.  For every
contrast this module computes the group-to-group difference on the frozen partition
*before* and *after* a candidate's correction, and reports whether the effect moved
by more than the tolerance.

The outcomes are kept distinct and honest:

* ``ok`` — the effect was computable before and after; the record carries the signed
  change and a verdict against the predeclared tolerance.
* ``insufficient_data`` — too few cells or too few independent units; never a zero.
* ``ineligible`` — the contrast is unassessable (confounded with the correction
  factor, or the groups are measured in disjoint batch sets).  Unassessable is never
  reported as preserved.

Where the design declares donor structure, the effect uses donors as units and, when
donors span both groups, a within-donor paired difference rather than a pooled
comparison.  The model used is stated in the record reason.
"""

from __future__ import annotations

import numpy as np

from batchdx.cyto.config import BiologicalContrastConfig
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.metrics.registry import BIOLOGICAL_CONTRAST
from batchdx.cyto.schemas.contracts import DesignAudit, MetricRecord, SamplingManifest
from batchdx.cyto.schemas.enums import MetricStatus

#: Minimum cells per (population, group) before a contrast effect is computed.  The
#: same floor as the other common-partition metrics, declared here rather than hidden.
MIN_GROUP_CELLS = 30
#: Minimum independent units per group.  Cells are not independent replicates; the
#: unit is the donor when a donor factor is declared, otherwise the sample.
MIN_UNITS_PER_GROUP = 2


def _group_level(sample, factor: str) -> str | None:
    """Resolve a sample-level factor from a manifest field or technical factors."""
    if hasattr(sample, factor):
        value = getattr(sample, factor)
        if value is None:
            return None
        return str(value.value) if hasattr(value, "value") else str(value)
    value = sample.technical_factors.get(factor)
    return None if value is None else str(value)


def _record(sampling: SamplingManifest, **kwargs) -> MetricRecord:
    return MetricRecord(
        metric_namespace=BIOLOGICAL_CONTRAST.namespace,
        metric_name=BIOLOGICAL_CONTRAST.name,
        metric_version=BIOLOGICAL_CONTRAST.version,
        stage=BIOLOGICAL_CONTRAST.stage,
        direction=BIOLOGICAL_CONTRAST.direction,
        scale=BIOLOGICAL_CONTRAST.scale,
        label_tier=sampling.label_source,
        label_confidence=sampling.label_confidence,
        representation=sampling.representation,
        **kwargs,
    )


def compute_biological_contrasts(
    dataset: Dataset,
    baseline: np.ndarray,
    candidate_matrix: np.ndarray,
    sampling: SamplingManifest,
    design: DesignAudit,
    contrasts: list[BiologicalContrastConfig],
    *,
    run_id: str,
    candidate_id: str,
    donor_factor: str | None = None,
) -> list[MetricRecord]:
    """Assess every prespecified contrast on the frozen partition, before and after.

    ``baseline`` and ``candidate_matrix`` must be aligned with the dataset's cell
    order, as every common-partition metric requires.  One record is produced per
    contrast per frozen population.
    """
    if not contrasts:
        return []

    obs = dataset.obs
    keep = obs[CELL_ID].isin(sampling.labels_by_cell.keys())
    obs = obs.loc[keep].copy()
    obs["_label"] = obs[CELL_ID].map(sampling.labels_by_cell)
    batch_of = dataset.sample_batches()
    obs["_batch"] = obs[SAMPLE_ID].map(batch_of).astype(str)
    positions = {cid: i for i, cid in enumerate(dataset.obs[CELL_ID])}
    obs["_row"] = obs[CELL_ID].map(positions)
    if donor_factor:
        donor_of = {sid: rec.donor_id for sid, rec in dataset.samples.items()}
        obs["_unit"] = obs[SAMPLE_ID].map(donor_of).astype(str)
    else:
        obs["_unit"] = obs[SAMPLE_ID].astype(str)

    populations = sorted(set(sampling.labels_by_cell.values()))
    records: list[MetricRecord] = []

    for contrast in contrasts:
        common = dict(
            run_id=run_id, candidate_id=candidate_id, marker=contrast.marker,
            comparison_id=f"contrast:{contrast.contrast_id}", aggregation_weight=1.0,
        )
        if contrast.marker not in sampling.evaluation_markers:
            records.append(_record(sampling, status=MetricStatus.INELIGIBLE, value=None,
                reason=(
                    f"marker {contrast.marker!r} is not flagged use_for_evaluation; a "
                    "contrast is never computed on an undeclared marker"
                ),
                **common,
            ))
            continue
        col = int(dataset.feature_index([contrast.marker])[0])
        group_of = {sid: _group_level(rec, contrast.group_factor)
                    for sid, rec in dataset.samples.items()}
        obs["_group"] = obs[SAMPLE_ID].map(group_of)

        for pop in populations:
            rows = obs.loc[obs["_label"] == pop]
            rows_a = rows.loc[rows["_group"] == contrast.group_a]
            rows_b = rows.loc[rows["_group"] == contrast.group_b]
            local = dict(common, population=pop)
            n_a, n_b = len(rows_a), len(rows_b)

            if contrast.group_factor in design.confounded_factors:
                records.append(_record(sampling, status=MetricStatus.INELIGIBLE,
                    value=None,
                    reason=(
                        f"unassessable: the contrast factor {contrast.group_factor!r} "
                        f"is perfectly confounded with the correction factor "
                        f"{design.correction_factor!r}; the group contrast and the "
                        "batch contrast are the same contrast, so preservation cannot "
                        "be attributed to the correction and is never reported as "
                        "preserved"
                    ),
                    **local,
                ))
                continue
            batches_a = set(rows_a["_batch"])
            batches_b = set(rows_b["_batch"])
            if not batches_a & batches_b:
                records.append(_record(sampling, status=MetricStatus.INELIGIBLE,
                    value=None,
                    reason=(
                        f"unassessable: groups {contrast.group_a!r} and "
                        f"{contrast.group_b!r} are measured in disjoint batch sets "
                        f"{sorted(batches_a)} and {sorted(batches_b)} for population "
                        f"{pop!r}; the comparison is non-identifiable and is never "
                        "reported as preserved"
                    ),
                    **local,
                ))
                continue
            if n_a < MIN_GROUP_CELLS or n_b < MIN_GROUP_CELLS:
                records.append(_record(sampling, status=MetricStatus.INSUFFICIENT_DATA,
                    value=None,
                    reason=(
                        f"population {pop!r} has {n_a} and {n_b} cells in groups "
                        f"{contrast.group_a!r} and {contrast.group_b!r}; "
                        f"{MIN_GROUP_CELLS} cells per group are required"
                    ),
                    **local,
                ))
                continue

            units_a = sorted(set(rows_a["_unit"]))
            units_b = sorted(set(rows_b["_unit"]))
            if len(units_a) < MIN_UNITS_PER_GROUP or len(units_b) < MIN_UNITS_PER_GROUP:
                kind = "donor(s)" if donor_factor else "sample(s)"
                records.append(_record(sampling, status=MetricStatus.INSUFFICIENT_DATA,
                    value=None,
                    reason=(
                        f"population {pop!r} has {len(units_a)} and {len(units_b)} "
                        f"independent {kind} in the two groups; cells are not "
                        f"independent replicates and {MIN_UNITS_PER_GROUP} per group "
                        "are required"
                    ),
                    **local,
                ))
                continue

            def unit_means(matrix: np.ndarray, rows_grp, units: list[str], col=col) -> np.ndarray:
                v = matrix[rows_grp["_row"].to_numpy(), col]
                u = rows_grp["_unit"].to_numpy()
                return np.array([float(v[u == unit].mean()) for unit in units])

            m_a0 = unit_means(baseline, rows_a, units_a)
            m_b0 = unit_means(baseline, rows_b, units_b)
            m_a1 = unit_means(candidate_matrix, rows_a, units_a)
            m_b1 = unit_means(candidate_matrix, rows_b, units_b)

            shared = sorted(set(units_a) & set(units_b))
            if donor_factor and len(shared) >= 2:
                model = (
                    f"donor-paired difference over {len(shared)} donors measured in "
                    "both groups"
                )
                e0 = float(np.mean([
                    m_b0[units_b.index(d)] - m_a0[units_a.index(d)] for d in shared
                ]))
                e1 = float(np.mean([
                    m_b1[units_b.index(d)] - m_a1[units_a.index(d)] for d in shared
                ]))
            else:
                kind = "donor-level" if donor_factor else "sample-level"
                detail = (
                    "no donor spans both groups" if donor_factor
                    else "no donor factor is declared"
                )
                model = f"pooled {kind} means ({detail})"
                e0 = float(m_b0.mean() - m_a0.mean())
                e1 = float(m_b1.mean() - m_a1.mean())

            change = e1 - e0
            tol = contrast.preservation_tolerance
            preserved = abs(change) <= tol
            records.append(_record(sampling, status=MetricStatus.OK, value=change,
                n_cells=n_a + n_b,
                n_donors=len(set(units_a) | set(units_b)) if donor_factor else None,
                n_samples=len(set(rows_a[SAMPLE_ID]) | set(rows_b[SAMPLE_ID])),
                reason=(
                    f"effect before {e0:.4g}, after {e1:.4g} ({model}); change "
                    f"{change:.4g} is {'within' if preserved else 'beyond'} the "
                    f"predeclared preservation tolerance {tol:.4g} (origin: "
                    f"{contrast.preservation_tolerance_origin}) -> "
                    + ("preserved" if preserved else "not_preserved")
                ),
                **local,
            ))

    return records
