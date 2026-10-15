"""Stage 9 acceptance gate: candidate states, guardrails and abstention."""

from __future__ import annotations

import numpy as np
import pytest

from batchdx.cyto.adapters import BackendContext, get_backend
from batchdx.cyto.decisions import POLICIES, assess_candidates
from batchdx.cyto.decisions.policy import Policy
from batchdx.cyto.metrics import compute_common_metrics, compute_reductions
from batchdx.cyto.schemas.contracts import MetricResult
from batchdx.cyto.schemas.enums import CorrectionStatus, DecisionState, LabelTier


def _score(frozen, backend_name, candidate_id, tmp_path):
    config, dataset, design, sampling = frozen
    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed)
    matrix, correction = get_backend(backend_name).run(dataset, candidate_id, {}, context)
    metrics = compute_common_metrics(dataset, matrix, sampling, design,
                                     run_id="t", candidate_id=candidate_id)
    return correction, metrics


def _decide(frozen, corrections, metrics, policy="conservative", **kwargs):
    _, _, design, sampling = frozen
    reductions = compute_reductions(metrics, sampling, run_id="t",
                                    identity_candidate_id="uncorrected")
    return assess_candidates(design, sampling, corrections, metrics, reductions,
                             run_id="t", policy_name=policy,
                             identity_candidate_id="uncorrected", **kwargs)


def _states(decision):
    return {a.candidate_id: a.state for a in decision.assessments}


def test_absent_batch_effect_yields_not_needed(frozen, tmp_path) -> None:
    config, dataset, design, sampling = frozen
    rng = np.random.default_rng(7)
    flat = np.tile(rng.normal(size=(1, dataset.adata.n_vars)), (dataset.n_cells, 1))
    flat += rng.normal(scale=1e-4, size=flat.shape)

    identity = compute_common_metrics(dataset, flat, sampling, design,
                                      run_id="t", candidate_id="uncorrected")
    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed)
    _, correction = get_backend("none").run(dataset, "uncorrected", {}, context)

    decision = _decide(frozen, {"uncorrected": correction}, {"uncorrected": identity})
    assert _states(decision)["uncorrected"] is DecisionState.NOT_NEEDED
    assert decision.recommended_candidate_id is None


def test_successful_correction_reaches_review_but_not_automatic_approval(
    frozen, tmp_path
) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    decision = _decide(
        frozen,
        {"uncorrected": identity_correction, "offset": correction},
        {"uncorrected": identity_metrics, "offset": metrics},
    )
    states = _states(decision)
    assert states["offset"] is DecisionState.ELIGIBLE_FOR_REVIEW
    assert "uncalibrated" in next(
        a.reason for a in decision.assessments if a.candidate_id == "offset"
    )
    assert decision.recommended_candidate_id is None
    assert any("uncalibrated" in w for w in decision.warnings)


def test_a_calibrated_policy_can_recommend(frozen, tmp_path) -> None:
    """With thresholds that have a declared, non-uncalibrated origin, approval is possible."""

    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    POLICIES["calibrated_for_test"] = Policy(
        name="calibrated_for_test",
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
        calibration_artifacts={
            k: "tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        decision = _decide(
            frozen,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics},
            policy="calibrated_for_test",
        )
        assert _states(decision)["offset"] is DecisionState.RECOMMENDED_UNDER_POLICY
        assert decision.recommended_candidate_id == "offset"
    finally:
        POLICIES.pop("calibrated_for_test")


def test_information_destroying_correction_is_rejected(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    correction, metrics = _score(frozen, "fixture_collapse", "collapse", tmp_path)
    decision = _decide(frozen, {"collapse": correction},
                       {"uncorrected": identity_metrics, "collapse": metrics})
    assessment = next(a for a in decision.assessments if a.candidate_id == "collapse")
    assert assessment.state is DecisionState.REJECTED_DISTORTION
    assert assessment.guardrails_failed


def test_sparse_validation_coverage_yields_insufficient_evidence(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)
    sparse = MetricResult(
        run_id="t", candidate_id="offset", metric_group="common",
        records=[r for r in metrics.records if r.comparison_id == "B1__vs__B2"],
        status="partial", reason="only one comparison retained for this test",
    )
    decision = _decide(frozen, {"offset": correction},
                       {"uncorrected": identity_metrics, "offset": sparse})
    assessment = next(a for a in decision.assessments if a.candidate_id == "offset")
    assert assessment.state is DecisionState.INSUFFICIENT_EVIDENCE
    assert "coverage" in assessment.reason


def test_confounded_design_is_non_identifiable_for_every_candidate(
    confounded_project, tmp_path
) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels

    config = load_config(confounded_project / "config" / "config.yaml")
    dataset = load_dataset(config, confounded_project)
    design = audit_design(
        list(dataset.samples.values()), dataset_id="c", correction_factor="batch",
        protected=["group"], donor_factor="donor_id",
        minimum_batches=config.design.minimum_batches,
        require_shared_reference=config.design.require_shared_reference,
    )
    labels, tier, confidence = resolve_labels(dataset, config, confounded_project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    metrics = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                     run_id="t", candidate_id="uncorrected")
    reductions = compute_reductions({"uncorrected": metrics}, sampling, run_id="t",
                                    identity_candidate_id="uncorrected")
    decision = assess_candidates(design, sampling, {}, {"uncorrected": metrics},
                                 reductions, run_id="t",
                                 identity_candidate_id="uncorrected")
    assert decision.assessments
    assert all(
        a.state in (DecisionState.BELOW_DESIGN_FLOOR, DecisionState.NON_IDENTIFIABLE)
        for a in decision.assessments
    )
    assert decision.recommended_candidate_id is None


def test_below_floor_design_blocks_every_recommendation(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    design.below_design_floor = True
    design.design_floor_reason = "forced for this test"
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    decision = _decide(frozen, {"offset": correction},
                       {"uncorrected": identity_metrics, "offset": metrics})
    assert all(a.state is DecisionState.BELOW_DESIGN_FLOOR for a in decision.assessments)


def test_failed_and_unavailable_backends_get_distinct_states(frozen, tmp_path) -> None:
    from batchdx.cyto.adapters.base import BackendError, failure_result, make_spec

    config, dataset, design, sampling = frozen
    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed)
    with pytest.raises(BackendError):
        get_backend("fixture_failure").run(dataset, "boom", {}, context)

    spec = make_spec("boom", "fixture_failure", None, {}, 473, dataset, [])
    failed = failure_result(spec, CorrectionStatus.FAILED, "it exploded")
    unavailable = failure_result(
        make_spec("absent", "cytonorm", None, {}, 473, dataset, []),
        CorrectionStatus.NOT_AVAILABLE, "not installed",
    )
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    decision = _decide(frozen, {"boom": failed, "absent": unavailable},
                       {"uncorrected": identity_metrics})
    states = _states(decision)
    assert states["boom"] is DecisionState.EXECUTION_FAILED
    assert states["absent"] is DecisionState.INSUFFICIENT_EVIDENCE


def test_dropping_a_failed_metric_cannot_improve_a_recommendation(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    correction, metrics = _score(frozen, "fixture_collapse", "collapse", tmp_path)

    full = _decide(frozen, {"collapse": correction},
                   {"uncorrected": identity_metrics, "collapse": metrics})
    assert _states(full)["collapse"] is DecisionState.REJECTED_DISTORTION

    # Remove exactly the guardrail metric that rejected it.  The candidate must not be
    # promoted by the absence of its own failing evidence.
    pruned = MetricResult(
        run_id="t", candidate_id="collapse", metric_group="common",
        records=[r for r in metrics.records
                 if r.metric_name not in ("mad_compression_fraction", "mad_signed_change")],
    )
    after = _decide(frozen, {"collapse": correction},
                    {"uncorrected": identity_metrics, "collapse": pruned})
    state = _states(after)["collapse"]
    assert state is not DecisionState.RECOMMENDED_UNDER_POLICY
    assert state in (
        DecisionState.REJECTED_DISTORTION,
        DecisionState.ELIGIBLE_FOR_REVIEW,
        DecisionState.INSUFFICIENT_EVIDENCE,
    )


def test_clustering_labels_alone_cannot_support_a_recommendation(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    sampling = sampling.model_copy(update={"label_source": LabelTier.CLUSTERING,
                                           "label_confidence": "low"})
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    POLICIES["calibrated_clustering_test"] = Policy(
        name="calibrated_clustering_test",
        threshold_origins={k: "calibrated for this test"
                           for k in Policy(name="x").threshold_origins},
    )
    try:
        reductions = compute_reductions(
            {"uncorrected": identity_metrics, "offset": metrics}, sampling,
            run_id="t", identity_candidate_id="uncorrected",
        )
        decision = assess_candidates(
            design, sampling, {"offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics}, reductions,
            run_id="t", policy_name="calibrated_clustering_test",
            identity_candidate_id="uncorrected",
        )
        assessment = next(a for a in decision.assessments if a.candidate_id == "offset")
        assert assessment.state is DecisionState.ELIGIBLE_FOR_REVIEW
        assert "batch-driven" in assessment.reason
    finally:
        POLICIES.pop("calibrated_clustering_test")


def test_policy_settings_and_threshold_origins_are_published(frozen, tmp_path) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    decision = _decide(frozen, {}, {"uncorrected": identity_metrics})
    assert set(decision.threshold_origins) >= {
        "min_coverage", "min_emd_reduction", "max_mad_compression_fraction",
        "max_joint_structure_change", "no_effect_emd",
    }
    assert decision.policy_settings["min_coverage"] == POLICIES["conservative"].min_coverage


def test_overrides_are_recorded_beside_the_original_assessment(frozen) -> None:
    _, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    override = {"candidate_id": "offset", "reason": "operator accepted the risk"}
    decision = _decide(frozen, {}, {"uncorrected": identity_metrics}, overrides=[override])
    assert decision.overrides == [override]
    assert decision.assessments, "the original assessment is still displayed"


def test_confounded_design_blocks_recommendation_even_when_reference_supported(
    frozen, tmp_path
) -> None:
    """Regression for the policy-layer confounding gate.

    A design whose batch factor is perfectly confounded with a protected biological
    factor must never reach ``recommended_under_policy`` — not via per-backend
    eligibility alone.  Before this gate existed, a shared measured reference made
    every comparison ``reference_supported`` and a calibrated policy could recommend a
    correction on a batch-by-group design.
    """
    _, dataset, design, sampling = frozen
    design.confounded_factors = ["group"]

    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    POLICIES["calibrated_for_test"] = Policy(
        name="calibrated_for_test",
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
        calibration_artifacts={
            k: "tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        decision = _decide(
            frozen,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics},
            policy="calibrated_for_test",
        )
        states = _states(decision)
        assert states["offset"] is DecisionState.NON_IDENTIFIABLE
        assert "confounded" in next(
            a.reason for a in decision.assessments if a.candidate_id == "offset"
        )
        assert decision.recommended_candidate_id is None
    finally:
        POLICIES.pop("calibrated_for_test")


def test_permissive_policy_cannot_recommend_on_clustering_labels(frozen, tmp_path) -> None:
    """The clustering-label gate is an invariant, not a policy option: no shipped or
    user policy may reach recommended_under_policy on clustering-derived labels."""
    _, dataset, design, sampling = frozen
    object.__setattr__(sampling, "label_source", LabelTier.CLUSTERING)

    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    POLICIES["calibrated_for_test"] = Policy(
        name="calibrated_for_test",
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
        calibration_artifacts={
            k: "tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        decision = _decide(
            frozen,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics},
            policy="permissive",
        )
        assert _states(decision)["offset"] is not DecisionState.RECOMMENDED_UNDER_POLICY
        assert decision.recommended_candidate_id is None
    finally:
        POLICIES.pop("calibrated_for_test")


def test_reference_gate_is_evaluated_on_the_evidence_actually_used(
    frozen, tmp_path
) -> None:
    """A candidate whose EMD evidence comes from conditional contrasts must not ride on
    an unrelated reference-supported pair.  The gate is scoped to the comparisons that
    produced the candidate's evidence."""
    _, dataset, design, sampling = frozen
    # Downgrade every comparison except one pair to conditional; the candidate's
    # evidence spans all pairs, so the gate must now fail.
    from batchdx.cyto.schemas.enums import ComparisonSupport

    for c in design.comparisons:
        if c.population is None and c.support is ComparisonSupport.REFERENCE_SUPPORTED:
            c.support = ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL
    design.comparisons[0].support = ComparisonSupport.REFERENCE_SUPPORTED

    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    POLICIES["calibrated_for_test"] = Policy(
        name="calibrated_for_test",
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
        calibration_artifacts={
            k: "tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        decision = _decide(
            frozen,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics},
            policy="calibrated_for_test",
        )
        assert _states(decision)["offset"] is not DecisionState.RECOMMENDED_UNDER_POLICY
    finally:
        POLICIES.pop("calibrated_for_test")


def test_a_contrast_only_the_training_split_can_estimate_blocks_approval(frozen, tmp_path):
    """The held-out rule is enforced by the decision layer, not merely available.

    `split_for_comparison` was implemented and unit-tested in isolation, which left the
    guarantee unenforced: a candidate fitted on the training split could still be
    recommended on a contrast that only the training split can supply.  This drives
    `assess_candidates` with a manifest whose validate split holds one batch only, and
    asserts the calibrated policy stops at review and says why.
    """

    config, dataset, design, sampling = frozen
    identity_metrics = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    identity_correction, _ = _score(frozen, "none", "uncorrected", tmp_path)
    correction, metrics = _score(frozen, "fixture_offset_removal", "offset", tmp_path)

    batches = sorted({rec.batch for rec in dataset.samples.values() if rec.batch})
    assert len(batches) >= 2, "the paired-reference fixture must have two batches"
    all_samples = sorted(dataset.samples)

    # Every sample of every batch except the last goes to train, so no split but train
    # can compare those two batches.  The last batch's samples go to validate, which
    # leaves validate holding exactly one batch and therefore no batch contrast at all.
    leaky = sampling.model_copy(update={"splits": {
        "split_unit": ["sample_id"],
        "train": [s for s in all_samples if dataset.samples[s].batch != batches[-1]],
        "tune": [],
        "validate": [s for s in all_samples if dataset.samples[s].batch == batches[-1]],
    }})

    POLICIES["calibrated_for_test"] = Policy(
        name="calibrated_for_test",
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
        calibration_artifacts={
            k: "tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        # First: with the real (well-stratified) manifest this same policy does approve,
        # so the refusal below is caused by the split and not by the policy.
        approved = _decide(
            frozen,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics},
            policy="calibrated_for_test", samples=dict(dataset.samples),
        )
        assert _states(approved)["offset"] is DecisionState.RECOMMENDED_UNDER_POLICY

        # Then: with validate holding one batch, approval must be withheld.
        reductions = compute_reductions(
            {"uncorrected": identity_metrics, "offset": metrics}, leaky,
            run_id="t", identity_candidate_id="uncorrected",
        )
        blocked = assess_candidates(
            design, leaky,
            {"uncorrected": identity_correction, "offset": correction},
            {"uncorrected": identity_metrics, "offset": metrics}, reductions,
            run_id="t", policy_name="calibrated_for_test",
            identity_candidate_id="uncorrected", samples=dict(dataset.samples),
        )
    finally:
        POLICIES.pop("calibrated_for_test")

    assessment = next(a for a in blocked.assessments if a.candidate_id == "offset")
    assert assessment.state is DecisionState.ELIGIBLE_FOR_REVIEW, assessment.state
    assert assessment.state is not DecisionState.RECOMMENDED_UNDER_POLICY
    assert blocked.recommended_candidate_id is None
    assert "not estimable on data held out" in assessment.reason
    assert "estimable in" in assessment.reason, (
        "the reason must name where the contrast could be estimated, so a reader can see "
        "it was train-only rather than unestimable anywhere"
    )
