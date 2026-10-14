"""Stage 8: biological preservation and uncertainty.

Prespecified biological contrasts (with donor-paired models where the design declares
donor structure), donor-block bootstrap intervals and restricted permutation are
implemented and tested here, alongside the joint-structure metric and the distortion
guardrails the Stage 8 gate already depended on.  The OpenProblems functional-marker
compatibility component remains deferred (see ``docs/implementation_status.md``);
nothing in this file claims otherwise.
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


# --------------------------------------------------------------------------------------
# Stage 8 proper: biological preservation and uncertainty
# --------------------------------------------------------------------------------------


def _contrast(marker="CD69", tolerance=0.2):
    from batchdx.cyto.config import BiologicalContrastConfig

    return BiologicalContrastConfig(
        contrast_id=f"{marker}_case_vs_control", marker=marker, group_a="case",
        group_b="control", preservation_tolerance=tolerance,
        preservation_tolerance_origin=(
            "test-declared: half of the 0.9 biological shift injected by the fixture, "
            "an order of magnitude above the sampling noise of a 30-cell group mean"
        ),
    )


def _make_project(spec, tmp_path, contrasts=()):
    """A simulated dataset with annotation labels and a fully declared design."""
    from batchdx.cyto._version import CONFIG_SCHEMA_VERSION
    from batchdx.cyto.config import (
        CandidateConfig,
        Config,
        CorrectionConfig,
        DesignConfig,
        EvaluationConfig,
        InputConfig,
        PreprocessingConfig,
        ProjectConfig,
    )
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
    from batchdx.cyto.simulate import write_simulation

    paths = write_simulation(spec, tmp_path)
    config = Config(
        schema_version=CONFIG_SCHEMA_VERSION,
        project=ProjectConfig(id="stage8"),
        input=InputConfig(
            format="h5ad", path=str(paths["h5ad"]), expression_layer="preprocessed",
            samples=str(paths["samples"]), markers=str(paths["markers"]),
            processing_manifest=str(paths["processing"]),
        ),
        preprocessing=PreprocessingConfig(),
        design=DesignConfig(
            correction_batch="batch", protect=["group"], donor="donor_id",
            require_shared_reference=False, biological_contrasts=list(contrasts),
        ),
        evaluation=EvaluationConfig(label_source="annotation", labels="true_population"),
        correction=CorrectionConfig(
            candidates=[CandidateConfig(id="uncorrected", backend="none")]
        ),
    )
    dataset = load_dataset(config, tmp_path)
    design = audit_design(
        list(dataset.samples.values()), dataset_id=dataset.manifest.dataset_id,
        correction_factor="batch", protected=["group"], donor_factor="donor_id",
        minimum_batches=3, require_shared_reference=False,
    )
    labels, tier, confidence = resolve_labels(dataset, config, tmp_path)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    return config, dataset, design, sampling


def test_a_correct_offset_removal_preserves_the_prespecified_biological_contrast(
    tmp_path
) -> None:
    from batchdx.cyto.adapters import BackendContext, get_backend
    from batchdx.cyto.metrics.biology import compute_biological_contrasts
    from batchdx.cyto.schemas.enums import MetricStatus
    from batchdx.cyto.simulate import paired_reference_spec

    contrast = _contrast()
    config, dataset, design, sampling = _make_project(
        paired_reference_spec(), tmp_path, [contrast]
    )
    baseline = dataset.expression()

    def verdicts(matrix, candidate_id):
        records = compute_biological_contrasts(
            dataset, baseline, matrix, sampling, design, [contrast],
            run_id="t", candidate_id=candidate_id, donor_factor="donor_id",
        )
        ok = [r for r in records if r.status is MetricStatus.OK]
        assert ok, f"no computable contrast record for {candidate_id}: {records}"
        return ok

    # before vs before: the change is zero by construction and must read as preserved
    for record in verdicts(baseline, "uncorrected"):
        assert record.reason.endswith("-> preserved")
        assert float(record.value) == 0.0

    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed)
    matrix, correction = get_backend("fixture_offset_removal").run(
        dataset, "offset", {}, context
    )
    assert correction.status.value == "success"

    for record in verdicts(matrix, "offset"):
        assert record.reason.endswith("-> preserved"), record.reason
        assert abs(float(record.value)) <= contrast.preservation_tolerance
        # donors in this fixture never span both groups, so the model must say pooled
        assert "pooled donor-level means" in record.reason, record.reason


def test_a_donor_paired_model_is_used_when_donors_span_groups(tmp_path) -> None:
    from batchdx.cyto.metrics.biology import compute_biological_contrasts
    from batchdx.cyto.schemas.enums import MetricStatus
    from batchdx.cyto.simulate import SampleSpec, SimulationSpec

    samples = []
    for batch in ("B1", "B2"):
        for k in range(2):
            for group in ("case", "control"):
                samples.append(SampleSpec(
                    sample_id=f"P{k}_{group}_{batch}", donor_id=f"P{k}", batch=batch,
                    group=group, n_cells=800,
                ))
    spec = SimulationSpec(
        samples=samples, seed=473, batch_shift={"B2": {"CD4": 0.5}},
        group_shift={"case": {"CD69": 0.9}},
    )
    config, dataset, design, sampling = _make_project(spec, tmp_path, [_contrast()])

    records = compute_biological_contrasts(
        dataset, dataset.expression(), dataset.expression(), sampling, design,
        [_contrast()], run_id="t", candidate_id="uncorrected", donor_factor="donor_id",
    )
    ok = [r for r in records if r.status is MetricStatus.OK]
    assert ok
    assert all("donor-paired difference" in r.reason for r in ok), [
        r.reason for r in ok
    ]

    metrics = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                     run_id="t", candidate_id="uncorrected")
    with_ci = [
        r for r in metrics.records
        if r.metric_name == "emd_within_population"
        and r.uncertainty_kind == "donor_block_bootstrap_ci"
    ]
    assert with_ci, "two donors per (population, batch) must yield bootstrap intervals"


def test_a_contrast_confounded_with_batch_is_unassessable_never_preserved(tmp_path) -> None:
    from batchdx.cyto.metrics.biology import compute_biological_contrasts
    from batchdx.cyto.schemas.enums import MetricStatus
    from batchdx.cyto.simulate import confounded_spec

    config, dataset, design, sampling = _make_project(
        confounded_spec(), tmp_path, [_contrast()]
    )
    assert design.confounded_factors == ["group"]

    records = compute_biological_contrasts(
        dataset, dataset.expression(), dataset.expression(), sampling, design,
        [_contrast()], run_id="t", candidate_id="any", donor_factor="donor_id",
    )
    assert records
    assert all(r.status is MetricStatus.INELIGIBLE for r in records)
    assert all(r.value is None for r in records)
    assert all("unassessable" in (r.reason or "") for r in records)
    assert not any((r.reason or "").endswith("-> preserved") for r in records)


def test_marginal_preserving_corruption_leaves_contrasts_intact_but_is_still_rejected(
    frozen, tmp_path
) -> None:
    """The contrast layer alone cannot see this corruption; the guardrail still is."""

    from batchdx.cyto.config import BiologicalContrastConfig
    from batchdx.cyto.decisions import assess_candidates
    from batchdx.cyto.metrics.biology import compute_biological_contrasts
    from batchdx.cyto.schemas.enums import DecisionState

    config, dataset, design, sampling = frozen
    contrast = BiologicalContrastConfig(
        contrast_id="case_vs_control_CD69", marker="CD69", group_a="case",
        group_b="control", preservation_tolerance=0.2,
        preservation_tolerance_origin="test-declared: half of the injected 0.9 shift",
    )
    baseline = dataset.expression()
    matrix, correction, scramble_metrics = _run(
        "fixture_joint_scramble", frozen, tmp_path, "scramble"
    )

    def preserved(mtrx):
        records = compute_biological_contrasts(
            dataset, baseline, mtrx, sampling, design, [contrast],
            run_id="t", candidate_id="x", donor_factor=config.design.donor,
        )
        ok = [r for r in records if r.status is MetricStatus.OK]
        assert ok
        return all(r.reason.endswith("-> preserved") for r in ok)

    assert preserved(baseline)
    assert not preserved(matrix), (
        "the within-(population, batch) scramble permutes values across groups and "
        "donors, so donor-level group means are destroyed and the contrast correctly "
        "reads not_preserved; the joint-structure guardrail catches the same "
        "corruption independently"
    )

    identity = compute_common_metrics(dataset, baseline, sampling, design,
                                      run_id="t", candidate_id="uncorrected")
    reductions = compute_reductions(
        {"uncorrected": identity, "scramble": scramble_metrics}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    decision = assess_candidates(
        design, sampling, {"scramble": correction},
        {"uncorrected": identity, "scramble": scramble_metrics}, reductions,
        run_id="t", policy_name="conservative", identity_candidate_id="uncorrected",
    )
    states = {a.candidate_id: a.state for a in decision.assessments}
    assert states["scramble"] is DecisionState.REJECTED_DISTORTION


def test_rare_population_yields_insufficient_data_not_zeros(tmp_path) -> None:
    from batchdx.cyto.metrics.biology import compute_biological_contrasts
    from batchdx.cyto.schemas.enums import MetricStatus
    from batchdx.cyto.simulate import rare_population_spec

    config, dataset, design, sampling = _make_project(rare_population_spec(), tmp_path)
    metrics = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                     run_id="t", candidate_id="uncorrected")

    rare = [r for r in metrics.records if r.population == "Rare"]
    assert rare, "the rare population must be reported, never silently dropped"
    insufficient = [r for r in rare if r.status is MetricStatus.INSUFFICIENT_DATA]
    assert insufficient, "rare groups below the cell floor must be insufficient_data"
    assert all(r.value is None for r in insufficient)
    assert not any(
        r.value == 0.0 for r in rare if r.status is not MetricStatus.OK
    ), "a missing value must never be reported as a zero"

    records = compute_biological_contrasts(
        dataset, dataset.expression(), dataset.expression(), sampling, design,
        [_contrast(marker="CD19")], run_id="t", candidate_id="uncorrected",
        donor_factor="donor_id",
    )
    rare_contrasts = [r for r in records if r.population == "Rare"]
    assert rare_contrasts
    assert all(
        r.status is MetricStatus.INSUFFICIENT_DATA and r.value is None
        for r in rare_contrasts
    ), [r.reason for r in rare_contrasts]


def test_bootstrap_ci_is_deterministic_under_the_project_seed(tmp_path) -> None:
    from batchdx.cyto.schemas.enums import MetricStatus
    from batchdx.cyto.simulate import paired_reference_spec

    config, dataset, design, sampling = _make_project(paired_reference_spec(), tmp_path)

    def intervals():
        metrics = compute_common_metrics(dataset, dataset.expression(), sampling, design,
                                         run_id="t", candidate_id="uncorrected")
        return metrics, {
            (r.population, r.marker, r.comparison_id):
                (r.uncertainty_low, r.uncertainty_high, r.uncertainty_kind)
            for r in metrics.records if r.metric_name == "emd_within_population"
        }

    metrics, first = intervals()
    _, second = intervals()
    assert first == second, "the bootstrap must be deterministic given the seed"

    with_ci = [v for v in first.values() if v[2] is not None]
    assert with_ci
    assert all(
        v[2] == "donor_block_bootstrap_ci" and v[0] <= v[1] for v in with_ci
    )
    ok_records = [
        r for r in metrics.records
        if r.metric_name == "emd_within_population" and r.status is MetricStatus.OK
    ]
    assert all(
        r.uncertainty_kind == "donor_block_bootstrap_ci" or "bootstrap" in (r.reason or "")
        for r in ok_records
    ), "every OK EMD record carries an interval or says why it has none"


def test_restricted_permutation_permutes_only_within_strata_and_is_deterministic() -> None:
    from batchdx.cyto.metrics.uncertainty import restricted_permutation_test

    rng = np.random.default_rng(7)
    n_donors = 12
    values, labels, strata = [], [], []
    for d in range(n_donors):
        for batch, shift in (("A", 0.0), ("B", 0.8)):
            values.extend(rng.normal(shift, 1.0, size=6).tolist())
            labels.extend([batch] * 6)
            strata.extend([f"donor{d}"] * 6)
    values = np.array(values)
    labels = np.array(labels, dtype=object)
    strata = np.array(strata, dtype=object)

    result = restricted_permutation_test(values, labels, strata, n_permutations=500,
                                         seed=473)
    assert result.p_value is not None and result.p_value < 0.05
    assert result.n_exchangeable_strata == n_donors

    null = values - 0.8 * (labels == "B")
    null_result = restricted_permutation_test(null, labels, strata, n_permutations=500,
                                              seed=473)
    assert null_result.p_value is not None and null_result.p_value > 0.3

    again = restricted_permutation_test(values, labels, strata, n_permutations=500,
                                        seed=473)
    assert again.p_value == result.p_value, "same seed must give the same p-value"

    # no stratum contains both batch labels: no valid exchange, no p-value
    blocked = restricted_permutation_test(
        values, labels, np.arange(len(labels)).astype(object), n_permutations=100,
        seed=473,
    )
    assert blocked.p_value is None and blocked.reason
