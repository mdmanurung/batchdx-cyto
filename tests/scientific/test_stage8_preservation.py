"""Stage 8 (partial): joint-structure detection and distortion rejection.

Stage 8 as specified — prespecified biological contrasts with design-appropriate
models, donor-block bootstrap, and the OpenProblems functional-marker compatibility
component — is **not implemented** in this build (see
``docs/implementation_status.md``).  What *is* implemented and tested here is the
joint-structure metric the Stage 8 gate depends on, plus the distortion guardrails
that stop a collapsing correction from being recommended.  Nothing in this file should
be read as evidence that biological preservation has been assessed.
"""

from __future__ import annotations

import numpy as np
import pytest

from batchdx.cyto.adapters import BackendContext, get_backend
from batchdx.cyto.decisions import assess_candidates
from batchdx.cyto.metrics import compute_common_metrics, compute_reductions
from batchdx.cyto.schemas.enums import DecisionState, MetricStatus


def _run(backend_name, frozen, tmp_path, candidate_id):
    config, dataset, design, sampling = frozen
    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed)
    matrix, correction = get_backend(backend_name).run(dataset, candidate_id, {}, context)
    metrics = compute_common_metrics(dataset, matrix, sampling, design,
                                     run_id="t", candidate_id=candidate_id)
    return matrix, correction, metrics


def _median(result, name):
    values = [float(r.value) for r in result.records
              if r.metric_name == name and r.status is MetricStatus.OK and r.value is not None]
    return float(np.median(values)) if values else None


def _marker_median(result, name, marker):
    values = [float(r.value) for r in result.records
              if r.metric_name == name and r.marker == marker
              and r.status is MetricStatus.OK and r.value is not None]
    return float(np.median(values)) if values else None


def test_a_correct_offset_removal_improves_alignment_without_distortion(
    frozen, tmp_path
) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    _, correction, metrics = _run("fixture_offset_removal", frozen, tmp_path, "offset")
    assert correction.status.value == "success"

    # The fixture injects the technical offset on CD4 and HLA-DR only, so that is where
    # the improvement has to appear. A median over every marker would be dominated by
    # the markers that never carried an offset, where a finite reference sample can only
    # add estimation noise.
    for marker in ("CD4", "HLADR"):
        before = _marker_median(identity_metrics, "emd_within_population", marker)
        after = _marker_median(metrics, "emd_within_population", marker)
        assert after < 0.5 * before, f"{marker}: {after:.3g} is not well below {before:.3g}"

    reductions = compute_reductions(
        {"uncorrected": identity_metrics, "offset": metrics}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    reduction = next(r for r in reductions if r.candidate_id == "offset")
    assert reduction.status is MetricStatus.OK and float(reduction.value) > 0.5

    compression = _median(metrics, "mad_compression_fraction")
    assert compression is not None and compression < 0.25, "spread must survive"
    joint = _median(metrics, "joint_structure_distance")
    assert joint is not None and joint < 0.1, "co-expression must survive"


def test_a_collapsing_correction_is_rejected_despite_perfect_alignment(
    frozen, tmp_path
) -> None:
    config, dataset, design, sampling = frozen
    identity = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                      run_id="t", candidate_id="uncorrected")
    _, collapse_correction, collapse_metrics = _run(
        "fixture_collapse", frozen, tmp_path, "collapse"
    )

    # alignment is excellent...
    assert _median(collapse_metrics, "emd_within_population") < _median(
        identity, "emd_within_population"
    )
    # ...and the dispersion guardrail catches it anyway
    assert _median(collapse_metrics, "mad_compression_fraction") > 0.9

    reductions = compute_reductions(
        {"uncorrected": identity, "collapse": collapse_metrics}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    decision = assess_candidates(
        design, sampling,
        {"collapse": collapse_correction},
        {"uncorrected": identity, "collapse": collapse_metrics},
        reductions, run_id="t", policy_name="conservative",
        identity_candidate_id="uncorrected",
    )
    states = {a.candidate_id: a.state for a in decision.assessments}
    assert states["collapse"] is DecisionState.REJECTED_DISTORTION
    assert decision.recommended_candidate_id is None


def test_marginal_preserving_joint_destroying_corruption_is_detected(
    frozen, tmp_path
) -> None:
    """Only the joint-structure metric can see this corruption; it must see it."""

    _, dataset, design, sampling = frozen
    identity = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                      run_id="t", candidate_id="uncorrected")
    matrix, correction, metrics = _run("fixture_joint_scramble", frozen, tmp_path, "scramble")

    def by_key(result, name):
        return {
            (r.population, r.marker, r.comparison_id): float(r.value)
            for r in result.records
            if r.metric_name == name and r.status is MetricStatus.OK
        }

    # every marginal distribution is preserved exactly, so EMD cannot move
    before, after = by_key(identity, "emd_within_population"), by_key(metrics, "emd_within_population")
    shared = [k for k in before if k in after]
    assert shared
    assert max(abs(after[k] - before[k]) for k in shared) < 1e-9

    # ...and so is dispersion
    mad_before = by_key(identity, "mad_signed_change")
    mad_after = by_key(metrics, "mad_signed_change")
    assert all(abs(mad_after[k]) < 1e-9 for k in mad_after if k in mad_before)

    # the joint-structure metric is the only thing that sees the damage
    joint_identity = _median(identity, "joint_structure_distance")
    joint_scrambled = _median(metrics, "joint_structure_distance")
    assert joint_identity == pytest.approx(0.0, abs=1e-9)
    assert joint_scrambled > 0.1

    reductions = compute_reductions(
        {"uncorrected": identity, "scramble": metrics}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    decision = assess_candidates(
        design, sampling, {"scramble": correction},
        {"uncorrected": identity, "scramble": metrics}, reductions,
        run_id="t", policy_name="conservative", identity_candidate_id="uncorrected",
    )
    states = {a.candidate_id: a.state for a in decision.assessments}
    assert states["scramble"] is DecisionState.REJECTED_DISTORTION


def test_biology_only_differences_are_not_labelled_technical(frozen) -> None:
    """A difference that lives between groups, not between batches, stays a biological one."""

    _, dataset, design, sampling = frozen
    obs = dataset.obs
    group_of = {sid: rec.group for sid, rec in dataset.samples.items()}
    groups = obs["sample_id"].map(group_of).to_numpy()
    marker = "CD69"
    col = dataset.feature_index([marker])[0]

    matrix = dataset.expression()
    matrix[groups == "case", col] += 3.0  # a large, purely biological difference

    metrics = compute_common_metrics(dataset, matrix, sampling, design,
                                     run_id="t", candidate_id="biology")

    # The design is balanced: both groups appear in every batch, so a between-group
    # difference does not become a between-batch difference.
    per_batch = [
        float(r.value) for r in metrics.records
        if r.metric_name == "emd_within_population" and r.marker == marker
        and r.status is MetricStatus.OK
    ]
    assert per_batch
    assert np.median(per_batch) < 1.0, (
        "a purely biological shift present in every batch must not inflate the "
        "batch-pair EMD"
    )
    supports = {c.support.value for c in design.comparisons if c.population is None}
    assert supports == {"reference_supported"}


def test_stage8_is_marked_unimplemented_rather_than_quietly_missing() -> None:
    from batchdx.cyto.metrics.registry import DEFERRED_METRICS

    assert any("functional_marker_preservation" in key for key in DEFERRED_METRICS)
    for key, reason in DEFERRED_METRICS.items():
        if "functional_marker_preservation" in key:
            assert "not implemented" in reason
