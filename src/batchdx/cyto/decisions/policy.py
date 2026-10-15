"""Eligibility and guardrails, not a universal weighted leaderboard.

A candidate reaches a recommendation only by clearing three independent gates:
sufficient comparison coverage, supported technical improvement, and no unacceptable
preservation failure.  There is no score that trades one against another, because the
weights that would do so are not calibrated.

Abstention and "no correction" are first-class outcomes.  Every threshold carries a
declared origin, and a threshold whose origin is *uncalibrated* can produce
``eligible_for_review`` but never ``recommended_under_policy``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from batchdx.cyto.io.hashing import canonical_hash
from batchdx.cyto.schemas.contracts import (
    CandidateAssessment,
    CorrectionResult,
    DecisionResult,
    DesignAudit,
    MetricRecord,
    MetricResult,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import (
    ComparisonSupport,
    CorrectionStatus,
    DecisionState,
    LabelTier,
    MetricStatus,
)


@dataclass(frozen=True)
class Policy:
    """A named, publishable set of gates."""

    name: str
    min_coverage: float = 0.6
    min_emd_reduction: float = 0.2
    max_mad_compression_fraction: float = 0.25
    max_joint_structure_change: float = 0.1
    #: EMD below which the identity candidate is taken to show no batch effect worth
    #: correcting.  Uncalibrated.
    no_effect_emd: float = 0.05
    require_reference_supported: bool = True
    threshold_origins: dict[str, str] = field(default_factory=lambda: {
        "min_coverage": "uncalibrated engineering default",
        "min_emd_reduction": "uncalibrated engineering default",
        "max_mad_compression_fraction": "uncalibrated engineering default",
        "max_joint_structure_change": "uncalibrated engineering default",
        "no_effect_emd": "uncalibrated engineering default",
    })
    #: Per-threshold calibration evidence: threshold name -> reference to a calibration
    #: artifact (path, URI or hash).  Origins are self-declared and the artifact is not
    #: re-verified here, but a policy that declares itself calibrated must at least name
    #: what it was calibrated against, for every threshold, by construction.
    calibration_artifacts: dict[str, str] = field(default_factory=dict)

    @property
    def any_threshold_uncalibrated(self) -> bool:
        for name, origin in self.threshold_origins.items():
            if "uncalibrated" in origin:
                return True
            if not self.calibration_artifacts.get(name):
                # An origin that merely avoids the word "uncalibrated" is not
                # calibration evidence; a calibrated policy must name its artifact.
                return True
        return False


POLICIES: dict[str, Policy] = {
    "conservative": Policy(name="conservative"),
    # Relaxes the evidence-quantity gates, never the attribution gates: clustering
    # labels can never support a recommendation under any policy, and neither can a
    # confounded design.  "Permissive" loosens how much evidence is required, not
    # what kind of evidence counts.
    "permissive": Policy(
        name="permissive",
        min_coverage=0.4,
        require_reference_supported=False,
    ),
}


def _summarise(records: list[MetricRecord], name: str) -> tuple[float | None, float]:
    """Median of the OK values for ``name`` and the fraction of records that were OK."""

    subset = [r for r in records if r.metric_name == name]
    if not subset:
        return None, 0.0
    ok = [r for r in subset if r.status is MetricStatus.OK and r.value is not None]
    if not ok:
        return None, 0.0
    weights = np.array([r.aggregation_weight for r in ok], dtype=float)
    values = np.array([float(r.value) for r in ok], dtype=float)
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights)
    median = float(values[int(np.searchsorted(cumulative, cumulative[-1] / 2.0))])
    return median, len(ok) / len(subset)


def _worst_marker_median(records: list[MetricRecord], name: str) -> float | None:
    """Median across comparisons within each marker, then the worst marker.

    A batch effect confined to one or two markers is still a batch effect. Taking a
    single median over every (population, marker, comparison) cell would average such
    an effect away against the many markers that are unaffected, so the summary is
    robust *within* a marker and sensitive *across* markers.
    """

    by_marker: dict[str, list[float]] = defaultdict(list)
    for r in records:
        if r.metric_name != name or r.status is not MetricStatus.OK or r.value is None:
            continue
        by_marker[r.marker or "*"].append(float(r.value))
    if not by_marker:
        return None
    return max(float(np.median(values)) for values in by_marker.values())


def _comparable_key(records: list[MetricRecord]) -> str:
    keys = sorted(
        f"{r.population}|{r.marker}|{r.comparison_id}"
        for r in records
        if r.metric_name == "emd_within_population" and r.status is MetricStatus.OK
    )
    return canonical_hash(keys)


def assess_candidates(
    design: DesignAudit,
    sampling: SamplingManifest,
    corrections: dict[str, CorrectionResult],
    metrics: dict[str, MetricResult],
    reductions: list[MetricRecord],
    *,
    run_id: str,
    policy_name: str = "conservative",
    identity_candidate_id: str = "uncorrected",
    overrides: list[dict[str, str]] | None = None,
    samples: dict[str, Any] | None = None,
) -> DecisionResult:
    """Assign one state per candidate and, only if every gate clears, a recommendation.

    ``samples`` is the sample-level record table, needed only when the frozen manifest
    carries a train/tune/validate split: the held-out check asks which batch contrasts
    the validation split can supply, and that is a question about the samples.  It is
    optional so that a caller with no split, and every unit test written before splits
    existed, keeps its previous behaviour exactly.
    """

    policy = POLICIES.get(policy_name)
    if policy is None:
        raise KeyError(f"unknown policy {policy_name!r}; available: {sorted(POLICIES)}")

    reduction_by_candidate = {r.candidate_id: r for r in reductions}
    warnings: list[str] = []
    assessments: list[CandidateAssessment] = []

    all_non_identifiable = bool(design.comparisons) and all(
        c.support is ComparisonSupport.NON_IDENTIFIABLE
        for c in design.comparisons if c.population is None
    )
    # Support status per population-free comparison, so the reference gate can be
    # evaluated on the comparisons a candidate's evidence actually rests on.
    support_by_comparison = {
        c.comparison_id: c.support for c in design.comparisons if c.population is None
    }

    identity_metrics = metrics.get(identity_candidate_id)
    baseline_emd = (
        _worst_marker_median(identity_metrics.records, "emd_within_population")
        if identity_metrics else None
    )
    no_effect = baseline_emd is not None and baseline_emd < policy.no_effect_emd

    # Which of the design's batch contrasts are estimable only on data that a
    # fit-and-apply candidate would have been fitted on.  Computed once, from the
    # frozen manifest, because it is a property of the design and the split rather
    # than of any one candidate.
    #
    # This is deliberately a *blanket* check over the design's comparisons rather than a
    # check on whichever comparison a given candidate happened to be scored on.  A
    # candidate that was fitted on the training split cannot be judged by a contrast
    # that only the training split can supply, no matter which row the number came
    # from, so the whole comparison set has to clear before any of it counts as
    # held-out evidence.  With no splits in the manifest (the pre-split case) there is
    # nothing to be held out from, so the check does not fire and the older behaviour
    # is unchanged.
    non_held_out = ""
    if sampling.splits and samples:
        from batchdx.cyto.sampling.splits import split_for_comparison

        blocked = []
        for comparison in design.comparisons:
            if comparison.population is not None:
                continue
            placement = split_for_comparison(comparison, sampling.splits, samples)
            if not placement.can_support_recommendation:
                blocked.append(
                    f"{placement.batch_a} vs {placement.batch_b} "
                    f"(estimable in {placement.evaluable_in or 'no split'})"
                )
        if blocked:
            non_held_out = "; ".join(blocked)

    for candidate_id in sorted(set(corrections) | set(metrics)):
        correction = corrections.get(candidate_id)
        result = metrics.get(candidate_id)
        records = result.records if result else []
        evidence: dict[str, float | None] = {}
        guardrails: list[str] = []

        if correction is not None and correction.status is CorrectionStatus.FAILED:
            assessments.append(CandidateAssessment(
                candidate_id=candidate_id, state=DecisionState.EXECUTION_FAILED,
                reason=correction.reason or "the backend failed",
            ))
            continue
        if correction is not None and correction.status in (
            CorrectionStatus.INELIGIBLE, CorrectionStatus.NOT_AVAILABLE
        ):
            assessments.append(CandidateAssessment(
                candidate_id=candidate_id, state=DecisionState.INSUFFICIENT_EVIDENCE,
                reason=correction.reason or "the backend did not run",
            ))
            continue

        # Which comparisons does this candidate's EMD evidence actually rest on, and is
        # every one of them reference-supported?  A global any() would let improvement
        # measured on conditional contrasts ride on an unrelated reference-supported pair.
        evidence_comparisons = {
            r.comparison_id for r in records
            if r.metric_name == "emd_within_population" and r.status is MetricStatus.OK
        }
        candidate_reference_supported = bool(evidence_comparisons) and all(
            support_by_comparison.get(cid) is ComparisonSupport.REFERENCE_SUPPORTED
            for cid in evidence_comparisons
        )

        reduction_record = reduction_by_candidate.get(candidate_id)
        reduction = (
            float(reduction_record.value)
            if reduction_record is not None
            and reduction_record.status is MetricStatus.OK
            and reduction_record.value is not None
            else None
        )
        coverage = (
            float(reduction_record.coverage)
            if reduction_record is not None and reduction_record.coverage is not None
            else None
        )
        emd_median, emd_ok = _summarise(records, "emd_within_population")
        joint_median, joint_ok = _summarise(records, "joint_structure_distance")
        compression, _ = _summarise(records, "mad_compression_fraction")

        evidence.update({
            "emd_within_population_median": emd_median,
            "emd_worst_marker_median": _worst_marker_median(records, "emd_within_population"),
            "emd_relative_reduction": reduction,
            "joint_structure_distance_median": joint_median,
            "mad_compression_fraction": compression,
        })
        coverage_map = {"emd_comparisons": coverage, "emd_records_ok": emd_ok,
                        "joint_structure_records_ok": joint_ok}

        if design.below_design_floor:
            state, reason = DecisionState.BELOW_DESIGN_FLOOR, (
                design.design_floor_reason or "the design is below the declared floor"
            )
        elif all_non_identifiable:
            state, reason = DecisionState.NON_IDENTIFIABLE, (
                "every batch contrast is non-identifiable given the protected factors"
            )
        elif candidate_id == identity_candidate_id and no_effect:
            state, reason = DecisionState.NOT_NEEDED, (
                f"the worst marker's median within-population EMD on the uncorrected "
                f"baseline is {baseline_emd:.3g}, below the declared no-effect threshold "
                f"{policy.no_effect_emd:.3g}; no marker shows a batch difference worth "
                "correcting"
            )
        elif candidate_id == identity_candidate_id:
            state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                "the identity candidate is always retained as a valid outcome"
            )
        else:
            if compression is not None and compression > policy.max_mad_compression_fraction:
                guardrails.append(
                    f"mad_compression_fraction {compression:.3g} exceeds "
                    f"{policy.max_mad_compression_fraction:.3g}"
                )
            if joint_median is not None and joint_median > policy.max_joint_structure_change:
                guardrails.append(
                    f"joint_structure_distance median {joint_median:.3g} exceeds "
                    f"{policy.max_joint_structure_change:.3g}"
                )
            if design.confounded_factors:
                # Design-level backstop, independent of per-backend eligibility rules: a
                # batch factor perfectly confounded with a protected biological factor
                # makes the batch contrast unattributable, whatever a measured reference
                # or a backend's own identifiability check says.
                state, reason = DecisionState.NON_IDENTIFIABLE, (
                    f"protected factor(s) {design.confounded_factors} are perfectly "
                    f"confounded with the batch factor {design.correction_factor!r}; no "
                    "automatic correction recommendation is possible for the affected "
                    "contrast, and adding the factor as a covariate does not help"
                )
            elif guardrails:
                state, reason = DecisionState.REJECTED_DISTORTION, "; ".join(guardrails)
            elif reduction is None:
                state, reason = DecisionState.INSUFFICIENT_EVIDENCE, (
                    "no relative reduction could be computed against the identity candidate"
                )
            elif coverage is None or coverage < policy.min_coverage:
                state, reason = DecisionState.INSUFFICIENT_EVIDENCE, (
                    "comparison coverage is unavailable for this candidate"
                    if coverage is None else
                    f"comparison coverage {coverage:.2f} is below the declared minimum "
                    f"{policy.min_coverage:.2f}"
                )
            elif compression is None or joint_median is None:
                # A guardrail metric with no usable record cannot be shown to be
                # satisfied; an unguarded candidate is insufficient evidence, not a
                # silent pass.
                missing = [name for name, v in (
                    ("mad_compression_fraction", compression),
                    ("joint_structure_distance", joint_median),
                ) if v is None]
                state, reason = DecisionState.INSUFFICIENT_EVIDENCE, (
                    f"preservation guardrail(s) {missing} have no usable record for this "
                    "candidate, so distortion cannot be ruled out"
                )
            elif reduction < policy.min_emd_reduction:
                state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                    f"technical improvement {reduction:.3g} is below the declared minimum "
                    f"{policy.min_emd_reduction:.3g}"
                )
            elif policy.require_reference_supported and not candidate_reference_supported:
                state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                    "no batch contrast is supported by a measured reference; the evidence "
                    "is conditional on comparable biological samples behaving identically"
                )
            elif sampling.label_source is LabelTier.CLUSTERING:
                # Invariant, not a policy option: labels derived from clustering of the
                # uncorrected baseline may themselves be batch-driven, so they can never
                # alone support a recommendation under any policy.
                state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                    "evaluation labels come from clustering of the uncorrected baseline and "
                    "may themselves be batch-driven; they cannot alone support a recommendation"
                )
            elif non_held_out:
                # A contrast estimated only on data the candidate was fitted on is not
                # held-out evidence, however good the number looks.  With a
                # fit-and-apply backend in play this is the difference between
                # measuring generalisation and measuring its own training set, so the
                # outcome is review rather than approval, and the reason names the
                # contrast so a reader can see which comparison is not held out.
                state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                    "the supporting batch contrast is not estimable on data held out from "
                    f"the fit: {non_held_out}"
                )
            elif policy.any_threshold_uncalibrated:
                state, reason = DecisionState.ELIGIBLE_FOR_REVIEW, (
                    "every gate passed, but the policy thresholds are uncalibrated, so the "
                    "outcome is review rather than automatic approval"
                )
            else:
                state, reason = DecisionState.RECOMMENDED_UNDER_POLICY, (
                    f"coverage {coverage:.2f}, improvement {reduction:.3g}, no guardrail "
                    "triggered" if coverage is not None else
                    f"improvement {reduction:.3g}, no guardrail triggered"
                )

        assessments.append(CandidateAssessment(
            candidate_id=candidate_id, state=state, reason=reason,
            guardrails_failed=guardrails, evidence=evidence, coverage=coverage_map,
            comparable_evidence_key=_comparable_key(records),
        ))

    recommended = [a.candidate_id for a in assessments
                   if a.state is DecisionState.RECOMMENDED_UNDER_POLICY]
    withheld = None
    if not recommended:
        withheld = (
            "no candidate cleared every gate under this policy; review the assessments and "
            "their reasons. Retaining the uncorrected data is a valid outcome."
        )
    elif len(recommended) > 1:
        # A tie is not a silent None: the report must say why no single candidate is named.
        withheld = (
            f"{len(recommended)} candidates cleared every gate under this policy "
            f"({', '.join(sorted(recommended))}); the policy does not break ties between "
            "fully qualified candidates. Choose on grounds the policy does not encode "
            "(cost, throughput, laboratory preference) and record the choice."
        )

    groups: dict[str, list[CandidateAssessment]] = defaultdict(list)
    for a in assessments:
        if a.comparable_evidence_key:
            groups[a.comparable_evidence_key].append(a)
    front: list[str] = []
    comparable = False
    largest = max(groups.values(), key=len, default=[])
    if len(largest) > 1:
        comparable = True
        front = _pareto_front(largest)
    elif groups:
        warnings.append(
            "candidates were not evaluated on comparable evidence, so no Pareto comparison "
            "is shown"
        )

    if policy.any_threshold_uncalibrated:
        warnings.append(
            "the policy thresholds in use are uncalibrated engineering defaults. They are "
            "not validated universal cutoffs and must not be cited as such."
        )

    return DecisionResult(
        run_id=run_id, policy=policy.name,
        policy_settings={
            "min_coverage": policy.min_coverage,
            "min_emd_reduction": policy.min_emd_reduction,
            "max_mad_compression_fraction": policy.max_mad_compression_fraction,
            "max_joint_structure_change": policy.max_joint_structure_change,
            "no_effect_emd": policy.no_effect_emd,
            "require_reference_supported": policy.require_reference_supported,
        },
        threshold_origins=dict(policy.threshold_origins),
        assessments=assessments,
        recommended_candidate_id=recommended[0] if len(recommended) == 1 else None,
        recommendation_withheld_reason=withheld,
        pareto_front=front, pareto_comparable=comparable,
        overrides=list(overrides or []), warnings=warnings,
    )


def _pareto_front(assessments: list[CandidateAssessment]) -> list[str]:
    """Non-dominated candidates over (reduction up, joint change down, compression down)."""

    def vector(a: CandidateAssessment) -> tuple[float, float, float] | None:
        r = a.evidence.get("emd_relative_reduction")
        j = a.evidence.get("joint_structure_distance_median")
        c = a.evidence.get("mad_compression_fraction")
        if r is None or j is None or c is None:
            return None
        return (float(r), -float(j), -float(c))

    points = {a.candidate_id: vector(a) for a in assessments}
    points = {k: v for k, v in points.items() if v is not None}
    front: list[str] = []
    for name, v in points.items():
        dominated = any(
            all(o >= s for o, s in zip(other, v, strict=True)) and any(
                o > s for o, s in zip(other, v, strict=True)
            )
            for key, other in points.items() if key != name
        )
        if not dominated:
            front.append(name)
    return sorted(front)
