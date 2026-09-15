"""Derived comparisons between a candidate and the identity candidate.

Reductions use the *same* eligible comparisons in numerator and denominator: a
comparison that is missing, ineligible or insufficient for either candidate is
excluded from both.  Dropping an unavailable comparison can therefore never improve a
candidate's standing.
"""

from __future__ import annotations

from collections import defaultdict

from batchdx.cyto.metrics.registry import EMD_REDUCTION, EMD_WITHIN_POPULATION
from batchdx.cyto.schemas.contracts import MetricRecord, MetricResult, SamplingManifest
from batchdx.cyto.schemas.enums import MetricStatus

#: Denominators at or below this are treated as undefined rather than stabilised.
DENOMINATOR_FLOOR = 1e-9


def _key(record: MetricRecord) -> tuple:
    return (record.population, record.marker, record.comparison_id)


def compute_reductions(
    results: dict[str, MetricResult],
    sampling: SamplingManifest,
    *,
    run_id: str,
    identity_candidate_id: str,
) -> list[MetricRecord]:
    """Relative EMD reduction per candidate against the identity candidate."""

    if identity_candidate_id not in results:
        return [MetricRecord(
            run_id=run_id, candidate_id="*", metric_namespace="common",
            metric_name=EMD_REDUCTION.name, metric_version=EMD_REDUCTION.version,
            stage=EMD_REDUCTION.stage, status=MetricStatus.INELIGIBLE,
            direction=EMD_REDUCTION.direction, scale=EMD_REDUCTION.scale,
            reason=(
                f"the identity candidate {identity_candidate_id!r} produced no metrics, so "
                "no relative reduction is defined for any candidate"
            ),
        )]

    by_candidate: dict[str, dict[tuple, MetricRecord]] = defaultdict(dict)
    for candidate_id, result in results.items():
        for record in result.records:
            if record.metric_name != EMD_WITHIN_POPULATION.name:
                continue
            by_candidate[candidate_id][_key(record)] = record

    identity = by_candidate.get(identity_candidate_id, {})
    out: list[MetricRecord] = []

    for candidate_id in sorted(by_candidate):
        candidate = by_candidate[candidate_id]
        usable = [
            k for k in identity
            if k in candidate
            and identity[k].status is MetricStatus.OK
            and candidate[k].status is MetricStatus.OK
            and identity[k].value is not None
            and candidate[k].value is not None
        ]
        n_eligible = len(identity)
        coverage = (len(usable) / n_eligible) if n_eligible else 0.0
        base = dict(
            run_id=run_id, candidate_id=candidate_id, metric_namespace="common",
            metric_name=EMD_REDUCTION.name, metric_version=EMD_REDUCTION.version,
            stage=EMD_REDUCTION.stage, direction=EMD_REDUCTION.direction,
            scale=EMD_REDUCTION.scale, coverage=coverage,
            label_tier=sampling.label_source, label_confidence=sampling.label_confidence,
            representation=sampling.representation,
        )
        if not usable:
            out.append(MetricRecord(
                **base, status=MetricStatus.INSUFFICIENT_DATA,
                reason="no comparison succeeded for both this candidate and the identity "
                       "candidate",
            ))
            continue
        denominator = sum(float(identity[k].value) for k in usable)
        numerator = denominator - sum(float(candidate[k].value) for k in usable)
        if denominator <= DENOMINATOR_FLOOR:
            out.append(MetricRecord(
                **base, status=MetricStatus.UNDEFINED,
                reason=(
                    f"the identity candidate's total EMD over the shared comparisons is "
                    f"{denominator:.3g}; a relative reduction is undefined and is not "
                    "stabilised with an epsilon"
                ),
            ))
            continue
        out.append(MetricRecord(
            **base, value=float(numerator / denominator), status=MetricStatus.OK,
            reason=(
                f"over {len(usable)} of {n_eligible} comparisons that succeeded for both "
                "candidates"
            ),
        ))
    return out
