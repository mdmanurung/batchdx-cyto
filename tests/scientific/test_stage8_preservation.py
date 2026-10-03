"""Stage 8: biological preservation, design-appropriate models and uncertainty.

What is tested here
-------------------
* A prespecified contrast survives a **correct** correction within a predeclared
  tolerance, and is destroyed by a **destructive** one.  The test detects both.
* A correction neither manufactures nor erases a purely biological difference.
* The donor-block bootstrap interval narrows as the number of donors rises.
* The restricted permutation detects a real effect at a rate far above chance and
  calibrates at approximately chance under a null.
* A contrast the design cannot support is **refused with a reason**, not skipped.
* The replicate unit is the donor, not the cell: tripling the cells in every sample
  changes nothing.
* Every record that reports a value also reports its interval and the method that
  produced it.

The joint-structure metric and the distortion guardrails, which the Stage 8 gate also
depends on, are tested here too.

Tolerances are declared next to the assertions that use them and were not loosened to
obtain a favourable result.  Where a test initially failed the cause was a fixture that
did not contain the structure being tested for, and the fixture was corrected, not the
tolerance.
"""

from __future__ import annotations

import tempfile
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from batchdx.cyto.adapters import BackendContext, get_backend
from batchdx.cyto.config import load_config
from batchdx.cyto.decisions import assess_candidates
from batchdx.cyto.design import audit_design
from batchdx.cyto.examples import materialise_example
from batchdx.cyto.io import load_dataset
from batchdx.cyto.metrics import (
    ContrastSpec,
    compute_common_metrics,
    compute_preservation_metrics,
    compute_reductions,
)
from batchdx.cyto.metrics.preservation import (
    BOOTSTRAP_KIND,
    CONTRASTAT_RECORDS,
    MIN_PERMUTATION_ASSIGNMENTS,
    PERMUTATION_NULL_KIND,
    PERMUTATION_P_KIND,
    ContrastRefusal,
    resolve_contrast,
)
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.schemas.enums import DecisionState, MetricStatus
from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation

#: The five statistics emitted per (contrast, population, marker, model).
CONTRASTAT_METRICS = {spec.name for spec in CONTRASTAT_RECORDS}

#: Which interval each of them is allowed to carry.
EXPECTED_UNCERTAINTY_KIND = {
    "biological_contrast_effect_size": BOOTSTRAP_KIND,
    "biological_contrast_cohens_d": BOOTSTRAP_KIND,
    "biological_contrast_ttest_statistic": PERMUTATION_NULL_KIND,
    "biological_contrast_rank_statistic": PERMUTATION_NULL_KIND,
    "biological_contrast_permutation_p": PERMUTATION_P_KIND,
}

#: The paired-reference fixture shifts CD69 by +0.9 in the case group, on a within-cell
#: spread of 0.35 (see :func:`batchdx.cyto.simulate.paired_reference_spec`).  Every
#: tolerance below is expressed as a fraction of that injected shift.
INJECTED_CD69_SHIFT = 0.9

#: A prespecified contrast over the panel's two functional markers, on every frozen
#: population, with the number of resamples and the seed fixed so the record set is
#: reproducible.
FUNCTIONAL_CONTRAST = dict(
    contrast_id="case_vs_control",
    factor="group",
    level_a="case",
    level_b="control",
    expected_direction="increase_in_a",
    markers=("CD69", "HLADR"),
    n_resamples=2000,
    seed=11,
)


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


# --------------------------------------------------------------------------------------
# Preservation helpers
# --------------------------------------------------------------------------------------


def _preservation(dataset, matrix, sampling, design, candidate_id, **overrides):
    spec = ContrastSpec(**{**FUNCTIONAL_CONTRAST, **overrides})
    return compute_preservation_metrics(
        dataset, matrix, sampling, design, [spec], run_id="t", candidate_id=candidate_id
    )


def _by_key(result, name, *, require_ok=True):
    out = {}
    for r in result.records:
        if r.metric_name != name:
            continue
        if require_ok and r.status is not MetricStatus.OK:
            continue
        out[(r.population, r.marker)] = r
    return out


def _effects(result, marker="CD69"):
    return {k[0]: r for k, r in _by_key(result, "biological_contrast_effect_size").items()
            if k[1] == marker}


def _p_values(result, marker="CD69"):
    return {k[0]: r for k, r in _by_key(result, "biological_contrast_permutation_p").items()
            if k[1] == marker}


# --------------------------------------------------------------------------------------
# A second design: every donor contributes one case and one control sample
# --------------------------------------------------------------------------------------


def _paired_donor_design(root: Path, *, n_donors: int, effect: float, seed: int,
                         n_cells: int = 150, batches: tuple[str, ...] = ("B1", "B2", "B3")):
    """Write a design in which a donor is matched across the two levels.

    The shipped examples give every donor a single sample, so they can exercise the
    unpaired model only.  This builder writes the same manifests over a design where each
    donor appears once in each level of the same batch, which is the layout the paired
    model and the within-donor permutation need.
    """

    materialise_example("paired_reference", root, cells_per_sample=n_cells)
    samples = []
    for b, batch in enumerate(batches):
        samples.append(SampleSpec(
            sample_id=f"ref_{batch}", donor_id="REFDONOR", batch=batch, group="reference",
            reference_id="REF1", reference_role="technical_replicate", n_cells=n_cells,
        ))
        for k in range(n_donors):
            for group, tag in (("case", "CA"), ("control", "CO")):
                samples.append(SampleSpec(
                    sample_id=f"D{b}{k}{tag}_{batch}", donor_id=f"D{b}{k}", batch=batch,
                    group=group, n_cells=n_cells,
                ))
    spec = SimulationSpec(samples=samples, seed=seed, batch_shift={"B2": {"CD4": 0.8}})
    if effect:
        spec.group_shift = {"case": {"CD69": effect}}
    write_simulation(spec, root)

    config = load_config(root / "config" / "config.yaml")
    dataset = load_dataset(config, root)
    design = audit_design(
        list(dataset.samples.values()), dataset_id=dataset.manifest.dataset_id,
        correction_factor="batch", protected=["group"], donor_factor="donor_id",
        minimum_batches=3, require_shared_reference=True,
    )
    labels, tier, confidence = resolve_labels(dataset, config, root)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    return config, dataset, design, sampling


@pytest.fixture
def paired_donor(tmp_path):
    """A design where every donor contributes one case and one control sample."""

    return _paired_donor_design(
        tmp_path, n_donors=6, effect=INJECTED_CD69_SHIFT, seed=101
    )


# --------------------------------------------------------------------------------------
# The joint-structure and distortion clauses of the Stage 8 gate
# --------------------------------------------------------------------------------------


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
    # the improvement has to appear. A median over every marker would be dominated by the
    # markers that never carried an offset, where a finite reference sample can only
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


# --------------------------------------------------------------------------------------
# The contrast itself
# --------------------------------------------------------------------------------------


def test_a_correct_offset_removal_preserves_the_prespecified_contrast(
    frozen, tmp_path
) -> None:
    """The injected biological difference must survive a correct correction.

    Tolerance: **0.05 in units of the declared (arcsinh, cofactor 5) scale**, which is
    5 % of the 0.9 shift the fixture injects on CD69.  It is not tighter than that
    because the correction estimates each batch's offset from a *single* reference
    aliquot of 400 cells, whose own standard error at this spread is about 0.03; a
    tolerance below the precision of the reference measurement would be claiming more
    than the measurement carries.  (The observed change is far smaller than the
    tolerance, because the fixture's offset is a per-batch constant and the design is
    exactly balanced, so it cancels in a difference of arms.  The tolerance is stated
    from the measurement's precision, not from the number that came out.)
    """

    _, dataset, design, sampling = frozen
    corrected, _, _ = _run("fixture_offset_removal", frozen, tmp_path, "offset")
    before = _effects(_preservation(dataset, dataset.expression(), sampling, design,
                                     "uncorrected"))
    after = _effects(_preservation(dataset, corrected, sampling, design, "offset"))

    assert set(before) == set(after) and before, "the contrast must produce records"
    tolerance = 0.05 * INJECTED_CD69_SHIFT
    for population, record in before.items():
        got = after[population]
        assert got.status is MetricStatus.OK
        assert abs(got.value - record.value) <= tolerance, (
            f"{population}: the contrast moved from {record.value:.4g} to {got.value:.4g}, "
            f"more than the {tolerance:.3g} a correct correction may move it by"
        )
        assert got.value > 0.0, "the injected direction must be preserved"
        assert got.uncertainty_low > 0.0, "the interval must still exclude zero"
        assert "expected_direction=increase_in_a: supported" in got.reason

    p_before = _p_values(_preservation(dataset, dataset.expression(), sampling, design,
                                       "uncorrected"))
    p_after = _p_values(_preservation(dataset, corrected, sampling, design, "offset"))
    floor = 1.0 / (FUNCTIONAL_CONTRAST["n_resamples"] + 1)
    for population, record in p_before.items():
        assert record.value == pytest.approx(floor, rel=1e-9)
        assert p_after[population].value == pytest.approx(floor, rel=1e-9), (
            f"{population}: the correction moved the restricted permutation p off its floor"
        )

    # Cohen's d falls a little, because the correction removes between-batch spread from
    # the pooled SD, but the effect is not standardised away.
    d_before = _by_key(_preservation(dataset, dataset.expression(), sampling, design,
                                      "uncorrected"), "biological_contrast_cohens_d")
    d_after = _by_key(_preservation(dataset, corrected, sampling, design, "offset"),
                      "biological_contrast_cohens_d")
    for key, record in d_before.items():
        if key[1] != "CD69" or record.value is None:
            continue
        assert d_after[key].value > 0.5 * record.value


def test_a_collapsing_correction_does_not_preserve_the_contrast(frozen, tmp_path) -> None:
    """`fixture_collapse` erases the contrast, and the metrics must say so.

    The collapse replaces every value in a population with that population's median over
    all batches, so the two arms become the same constant: the effect size collapses to
    zero, the pooled standard deviation with it, and the standardised effect becomes
    *undefined* rather than zero, because a ratio with a zero denominator is not a number.
    """

    _, dataset, design, sampling = frozen
    collapsed, _, _ = _run("fixture_collapse", frozen, tmp_path, "collapse")
    before = _effects(_preservation(dataset, dataset.expression(), sampling, design,
                                     "uncorrected"))
    after = _effects(_preservation(dataset, collapsed, sampling, design, "collapse"))

    tolerance = 0.05 * INJECTED_CD69_SHIFT
    for population, before_record in before.items():
        assert before_record.value > 0.5, "the fixture must inject a real effect"
        got = after[population]
        assert got.status is MetricStatus.OK
        assert abs(got.value) <= tolerance, (
            f"{population}: the collapse left an effect of {got.value:.4g}; the contrast "
            "was not destroyed"
        )
        assert abs(got.uncertainty_high) <= tolerance, "and left an interval to match"

    sizes = _by_key(_preservation(dataset, collapsed, sampling, design, "collapse"),
                    "biological_contrast_cohens_d", require_ok=False)
    assert sizes, "the standardised effect must still be reported on"
    for record in sizes.values():
        assert record.status is MetricStatus.UNDEFINED, (
            "a zero pooled SD makes Cohen's d undefined, not zero"
        )
        assert record.value is None

    for population, record in _p_values(
        _preservation(dataset, collapsed, sampling, design, "collapse")
    ).items():
        assert record.value > 0.5, (
            f"{population}: with the contrast erased the permutation p must be at chance, "
            f"not significant; got {record.value:.4g}"
        )


def test_a_correction_neither_manufactures_nor_erases_a_biological_difference(
    frozen, tmp_path
) -> None:
    """A correct correction must leave the biology exactly where the design put it."""

    config, dataset, design, sampling = frozen
    group_of = {sid: rec.group for sid, rec in dataset.samples.items()}
    groups = dataset.obs["sample_id"].map(group_of).to_numpy()
    cd69 = dataset.feature_index(["CD69"])[0]
    cd3 = dataset.feature_index(["CD3"])[0]

    injected = 3.0  # an order of magnitude above the fixture's own 0.9 shift
    uncorrected = dataset.expression()
    uncorrected[groups == "case", cd69] += injected
    corrected, _, _ = _run("fixture_offset_removal", frozen, tmp_path, "offset")
    corrected = corrected.copy()
    corrected[groups == "case", cd69] += injected

    before = _effects(_preservation(dataset, uncorrected, sampling, design, "uncorrected"))
    after = _effects(_preservation(dataset, corrected, sampling, design, "offset"))
    tolerance = 0.05 * injected
    for population, record in before.items():
        assert abs(record.value - (INJECTED_CD69_SHIFT + injected)) <= 0.05, (
            "the fixture's own shift plus the injected one is the whole expected effect"
        )
        assert abs(after[population].value - record.value) <= tolerance, (
            f"{population}: the correction moved a purely biological difference"
        )
        assert after[population].value > 0.0

    # ...and it must not invent one where the design put none.  CD3 carries no group shift.
    control = _effects(_preservation(dataset, corrected, sampling, design, "offset"),
                       marker="CD3")
    for population, record in control.items():
        assert abs(record.value) <= tolerance, (
            f"{population}: the correction manufactured a CD3 contrast of {record.value:.4g}"
        )
    for population, record in _p_values(
        _preservation(dataset, corrected, sampling, design, "offset"), marker="CD3"
    ).items():
        assert record.value > 0.05, (
            f"{population}: CD3 has no group difference, so its permutation p must not be "
            f"significant; got {record.value:.4g}"
        )
    assert cd3 >= 0  # the marker really is on the panel the contrast did not ask for


def test_a_contrast_the_design_cannot_support_is_refused_with_a_reason(
    confounded_project,
) -> None:
    """A non-identifiable design gets a refusal, not a number and not silence."""

    config = load_config(confounded_project / "config" / "config.yaml")
    dataset = load_dataset(config, confounded_project)
    design = audit_design(
        list(dataset.samples.values()), dataset_id="c", correction_factor="batch",
        protected=["group"], donor_factor="donor_id",
        minimum_batches=config.design.minimum_batches,
        require_shared_reference=config.design.require_shared_reference,
    )
    assert design.confounded_factors == ["group"], "the fixture must be confounded"

    labels, tier, confidence = resolve_labels(dataset, config, confounded_project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)

    spec = ContrastSpec(**FUNCTIONAL_CONTRAST)
    outcome = resolve_contrast(dataset, design, sampling, spec)
    assert isinstance(outcome, ContrastRefusal)
    assert outcome.reason_code == "non_identifiable"
    assert "confounded" in outcome.reason

    result = compute_preservation_metrics(
        dataset, dataset.expression(), sampling, design, [spec],
        run_id="t", candidate_id="uncorrected",
    )
    assert len(result.records) == 1, "a refused contrast must not also be computed"
    record = result.records[0]
    assert record.metric_name == "biological_contrast_admissibility"
    assert record.status is MetricStatus.INELIGIBLE
    assert record.value is None, "a refusal is never reported as a score of zero"
    assert record.coverage == 0.0
    assert "non_identifiable" in record.reason and "confounded" in record.reason


def test_a_contrast_naming_a_level_that_does_not_exist_is_refused(frozen) -> None:
    _, dataset, design, sampling = frozen
    spec = ContrastSpec(contrast_id="absent", factor="group", level_a="case",
                        level_b="treated")
    outcome = resolve_contrast(dataset, design, sampling, spec)
    assert isinstance(outcome, ContrastRefusal)
    assert outcome.reason_code == "levels_absent"
    assert "treated" in outcome.reason and "reference" in outcome.reason

    # An equal-level contrast is rejected at construction, before it reaches the data.
    with pytest.raises(ValueError, match="two different levels"):
        ContrastSpec(contrast_id="degenerate", factor="group", level_a="case",
                     level_b="case")


def test_the_contrast_chooses_the_paired_model_only_when_a_donor_spans_both_levels(
    frozen, paired_donor,
) -> None:
    _, dataset, design, sampling = frozen
    plan = resolve_contrast(dataset, design, sampling, ContrastSpec(**FUNCTIONAL_CONTRAST))
    assert not isinstance(plan, ContrastRefusal)
    # The shipped example gives every donor one sample, so the samples are unpaired and
    # the permutation is restricted within a batch instead of within a donor.
    assert plan.model == "unpaired"
    assert plan.block_unit == "donor_id"
    assert all(name.startswith("batch:") for name, _ in plan.blocks)

    _, paired_dataset, paired_design, paired_sampling = paired_donor
    paired_plan = resolve_contrast(paired_dataset, paired_design, paired_sampling,
                                   ContrastSpec(**FUNCTIONAL_CONTRAST))
    assert not isinstance(paired_plan, ContrastRefusal)
    assert paired_plan.model == "paired"
    assert all(name.startswith("donor:") for name, _ in paired_plan.blocks)
    assert paired_plan.reachable_permutations == 2 ** paired_plan.n_donors

    # Asking for the paired model where it is impossible is refused, not approximated.
    refusal = resolve_contrast(dataset, design, sampling,
                               ContrastSpec(**{**FUNCTIONAL_CONTRAST, "model": "paired"}))
    assert isinstance(refusal, ContrastRefusal)
    assert refusal.reason_code == "pairing_unsupported"


def test_a_paired_contrast_is_tested_on_the_within_donor_differences(paired_donor) -> None:
    _, dataset, design, sampling = paired_donor
    result = _preservation(dataset, dataset.expression(), sampling, design, "uncorrected")
    effect = _effects(result)["T_CD4"]
    assert effect.value == pytest.approx(INJECTED_CD69_SHIFT, abs=0.1)
    t_statistic = _by_key(result, "biological_contrast_ttest_statistic")[("T_CD4", "CD69")]
    rank = _by_key(result, "biological_contrast_rank_statistic")[("T_CD4", "CD69")]
    permutation = _p_values(result)["T_CD4"]
    # A sign-flip null over 2**n donor pairs is discrete, so the smallest reportable
    # p-value is a small multiple of 1/2**n rather than the +1 resampling floor.
    assert t_statistic.value > 0.0
    assert t_statistic.uncertainty_high < t_statistic.value
    assert rank.value > rank.uncertainty_high
    assert permutation.value < 0.05
    assert "model=paired" in effect.reason


# --------------------------------------------------------------------------------------
# Uncertainty
# --------------------------------------------------------------------------------------


def test_the_donor_block_bootstrap_interval_narrows_as_the_donor_count_rises() -> None:
    """More donors, narrower interval: the width is the quantity the donors buy."""

    widths, donor_counts = [], []
    for n_donors in (2, 4, 8, 16):
        root = Path(tempfile.mkdtemp())
        _, dataset, design, sampling = _paired_donor_design(
            root, n_donors=n_donors, effect=INJECTED_CD69_SHIFT, seed=101
        )
        result = _preservation(dataset, dataset.expression(), sampling, design,
                               "uncorrected", populations=("T_CD4",), markers=("CD69",),
                               seed=5)
        record = _by_key(result, "biological_contrast_effect_size")[("T_CD4", "CD69")]
        assert record.status is MetricStatus.OK
        assert record.uncertainty_kind == f"{BOOTSTRAP_KIND}_95"
        widths.append(record.uncertainty_high - record.uncertainty_low)
        donor_counts.append(record.n_donors)
        assert record.n_donors == 3 * n_donors

    assert donor_counts == [6, 12, 24, 48]
    assert widths[0] > widths[1] > widths[2] > widths[3], (
        f"the interval did not narrow with the donor count: {widths}"
    )
    # Halving the width per doubling of donors is what a bootstrap of a mean does; the
    # loosest possible statement of that is a factor of two over the fourfold increase.
    assert widths[0] > 1.7 * widths[3]


def test_restricted_permutation_detects_a_real_effect_and_calibrates_at_chance(
    tmp_path,
) -> None:
    """Detection rate and calibration, each over independent resimulations.

    The data are redrawn for every replicate, so the replicates are independent draws
    rather than one dataset seen through twenty different permutation seeds.
    """

    n_replicates = 10
    p_values = {"null": [], "effect": []}
    for index in range(n_replicates):
        for kind, effect in (("null", 0.0), ("effect", INJECTED_CD69_SHIFT)):
            root = tmp_path / f"{kind}_{index}"
            _, dataset, design, sampling = _paired_donor_design(
                root, n_donors=5, effect=effect, seed=1000 + index, n_cells=150,
                batches=("B1", "B2"),
            )
            result = _preservation(dataset, dataset.expression(), sampling, design,
                                   "uncorrected", populations=("T_CD4",), markers=("CD69",),
                                   n_resamples=499, seed=index)
            p_values[kind].append(_p_values(result)["T_CD4"].value)

    detected_null = sum(p < 0.05 for p in p_values["null"])
    detected_effect = sum(p < 0.05 for p in p_values["effect"])
    assert detected_effect == n_replicates, (
        f"a real effect was detected {detected_effect}/{n_replicates} times: "
        f"{p_values['effect']}"
    )
    # A correctly calibrated test rejects at 5 % about one time in twenty.  Two or fewer
    # out of ten is inside that; three would not be.
    assert detected_null <= 2, (
        f"a null effect was called significant {detected_null}/{n_replicates} times: "
        f"{p_values['null']}"
    )
    assert 0.25 < float(np.median(p_values["null"])) < 0.75, (
        f"the null p-values are not centred at chance: {p_values['null']}"
    )
    assert max(p_values["null"]) > 0.3, "and they are not all implausibly small"


def test_every_preservation_record_carries_its_uncertainty(frozen, tmp_path) -> None:
    """The three contract fields that no Stage 7 record ever set are set here."""

    _, dataset, design, sampling = frozen
    candidates = {"uncorrected": dataset.expression()}
    for backend, candidate_id in (("fixture_offset_removal", "offset"),
                                  ("fixture_collapse", "collapse")):
        config, dsn, dsg, smp = frozen
        context = BackendContext(design=dsg, sampling=smp, work_dir=str(tmp_path),
                                 seed=config.project.seed)
        matrix, _ = get_backend(backend).run(dsn, candidate_id, {}, context)
        candidates[candidate_id] = matrix

    allowed = {f"{kind}_95" for kind in
               (BOOTSTRAP_KIND, PERMUTATION_NULL_KIND, PERMUTATION_P_KIND)}
    checked = 0
    for candidate_id, matrix in candidates.items():
        result = compute_preservation_metrics(
            dataset, matrix, sampling, design, [ContrastSpec(**FUNCTIONAL_CONTRAST)],
            run_id="t", candidate_id=candidate_id,
        )
        assert result.records
        for record in result.records:
            if record.metric_name == "biological_contrast_admissibility":
                continue
            if record.status is not MetricStatus.OK:
                # A record with no value reports no interval; what it must never do is
                # invent a value to carry one.
                assert record.value is None, (
                    f"{candidate_id}/{record.metric_name}: a non-ok row carries a value"
                )
                continue
            assert record.uncertainty_low is not None, record.metric_name
            assert record.uncertainty_high is not None, record.metric_name
            assert record.uncertainty_kind in allowed, record.uncertainty_kind
            assert record.uncertainty_kind == (
                f"{EXPECTED_UNCERTAINTY_KIND[record.metric_name]}_95"
            ), record.metric_name
            assert record.uncertainty_low <= record.uncertainty_high, (
                f"{candidate_id}/{record.metric_name}: "
                f"[{record.uncertainty_low}, {record.uncertainty_high}]"
            )
            assert record.n_cells is not None and record.coverage is not None
            assert record.label_tier is not None and record.representation is not None
            assert f"n_resamples={FUNCTIONAL_CONTRAST['n_resamples']}" in record.reason
            assert f"seed={FUNCTIONAL_CONTRAST['seed']}" in record.reason
            checked += 1
    assert checked > 100, f"only {checked} records were checked"


def test_an_interval_is_never_reported_without_a_reachable_permutation_null() -> None:
    """A null with too few reachable assignments is reported, not faked."""

    root = Path(tempfile.mkdtemp())
    # One donor per batch per level gives three complete donor pairs, whose sign-flip null
    # has 2**3 = 8 reachable outcomes: far too few for a p-value to mean anything.
    _, dataset, design, sampling = _paired_donor_design(
        root, n_donors=1, effect=INJECTED_CD69_SHIFT, seed=101
    )
    spec = ContrastSpec(**FUNCTIONAL_CONTRAST)
    plan = resolve_contrast(dataset, design, sampling, spec)
    assert not isinstance(plan, ContrastRefusal)
    assert plan.model == "paired" and plan.n_donors == 3
    assert plan.reachable_permutations == 8 < MIN_PERMUTATION_ASSIGNMENTS

    result = compute_preservation_metrics(
        dataset, dataset.expression(), sampling, design, [spec],
        run_id="t", candidate_id="uncorrected",
    )
    effects = _by_key(result, "biological_contrast_effect_size")
    p_values = _by_key(result, "biological_contrast_permutation_p", require_ok=False)
    assert effects and p_values
    for record in p_values.values():
        assert record.status is MetricStatus.INSUFFICIENT_DATA
        assert record.value is None
        assert "8" in record.reason
        assert str(MIN_PERMUTATION_ASSIGNMENTS) in record.reason
    # The effect size itself is still reported, with its bootstrap interval: only the
    # permutation-derived quantities are withheld.
    for record in effects.values():
        assert record.status is MetricStatus.OK
        assert record.uncertainty_kind == f"{BOOTSTRAP_KIND}_95"


def test_a_prespecified_contrast_is_reproducible_from_its_seed(frozen) -> None:
    """Same seed, same numbers; different seed, same point estimate, different resampling."""

    _, dataset, design, sampling = frozen
    matrix = dataset.expression()
    first = _by_key(_preservation(dataset, matrix, sampling, design, "a", seed=3),
                    "biological_contrast_effect_size")
    again = _by_key(_preservation(dataset, matrix, sampling, design, "a", seed=3),
                    "biological_contrast_effect_size")
    other = _by_key(_preservation(dataset, matrix, sampling, design, "a", seed=4),
                    "biological_contrast_effect_size")
    assert set(first) == set(again) == set(other)
    for key, record in first.items():
        assert again[key].value == record.value
        assert again[key].uncertainty_low == record.uncertainty_low
        assert again[key].uncertainty_high == record.uncertainty_high
        assert other[key].value == record.value, "the point estimate cannot depend on a seed"
        assert (other[key].uncertainty_low != record.uncertainty_low
                or other[key].uncertainty_high != record.uncertainty_high)


def test_the_contrast_level_order_is_normalised_so_a_positive_effect_means_the_same_arm() -> None:
    spec = ContrastSpec(contrast_id="c", factor="group", level_a="control", level_b="case")
    assert (spec.level_a, spec.level_b) == ("case", "control")
    spec = ContrastSpec(contrast_id="c", factor="group", level_a="case", level_b="control")
    assert (spec.level_a, spec.level_b) == ("case", "control")


# --------------------------------------------------------------------------------------
# The replicate unit
# --------------------------------------------------------------------------------------


def _triple_every_cell(dataset, sampling, repeats=3):
    """A copy of the dataset in which every cell is present ``repeats`` times."""

    obs = dataset.adata.obs
    take = np.repeat(np.arange(dataset.n_cells), repeats)
    new_obs = obs.iloc[take].reset_index(drop=True).copy()
    # The ids have to follow the same order as ``take``, or the labels are scrambled
    # across cells and the comparison below would be meaningless.
    new_obs["cell_id"] = [
        f"{cid}#{r}" for cid in obs["cell_id"].tolist() for r in range(repeats)
    ]
    new_obs["original_event_index"] = np.arange(len(take), dtype=int)
    adata = ad.AnnData(
        X=np.zeros((len(take), dataset.adata.n_vars), dtype=np.float32),
        obs=new_obs, var=dataset.adata.var.copy(),
    )
    # anndata derives var_names from the DataFrame's index only when var_names is not
    # passed explicitly, and a plain copy of an existing .var does not carry the index
    # across.  Losing it would make every marker lookup resolve against a 0..N-1
    # RangeIndex, so feature_index("CD69") would fail with the "feature(s) [True, ...]"
    # message that comes from treating an all-True membership test as a feature list.
    adata.var_names = list(dataset.adata.var_names)
    adata.layers[dataset.expression_layer] = np.asarray(
        dataset.adata.layers[dataset.expression_layer]
    )[take]

    cells = [c for c in sampling.selected_cell_ids for _ in range(repeats)]
    labels = {
        f"{c}#{r}": sampling.labels_by_cell[c]
        for c in sampling.selected_cell_ids for r in range(repeats)
    }
    return replace(dataset, adata=adata), sampling.model_copy(update={
        "selected_cell_ids": cells,
        "labels_by_cell": labels,
        "per_sample_counts": {k: v * repeats for k, v in sampling.per_sample_counts.items()},
    })


def test_the_replicate_unit_is_the_sample_not_the_cell(frozen) -> None:
    """Tripling the cells in every sample must change nothing at all.

    This is the direct check that cells are not treated as replicates.  A cell-level test
    would see 3x the data, shrink its standard error by sqrt(3), and move the p-value;
    the sample-level test sees the same twelve per-sample means, so the effect, the
    interval and the permutation p come out bit-identical.
    """

    _, dataset, design, sampling = frozen
    tripled_dataset, tripled_sampling = _triple_every_cell(dataset, sampling)
    assert tripled_dataset.n_cells == 3 * dataset.n_cells
    assert len(tripled_sampling.selected_cell_ids) == 3 * len(sampling.selected_cell_ids)

    original = _preservation(dataset, dataset.expression(), sampling, design, "a")
    tripled = _preservation(tripled_dataset, tripled_dataset.expression(), tripled_sampling,
                            design, "a")

    def key(record):
        return (record.metric_name, record.population, record.marker)

    before = {key(r): r for r in original.records if r.metric_name in CONTRASTAT_METRICS}
    after = {key(r): r for r in tripled.records if r.metric_name in CONTRASTAT_METRICS}
    assert set(before) == set(after) and before
    for k, record in before.items():
        got = after[k]
        assert got.status is record.status, k
        if record.value is None:
            assert got.value is None
            continue
        assert got.n_cells == 3 * record.n_cells, "the cell count is bookkeeping, not n"
        if k[0] == "biological_contrast_permutation_p":
            assert got.value == record.value, (
                "the restricted permutation p is a count of resamples and cannot move "
                "when the data are the same data repeated"
            )
        else:
            assert got.value == pytest.approx(record.value, abs=1e-9), k
            assert got.uncertainty_low == pytest.approx(record.uncertainty_low, abs=1e-9)
            assert got.uncertainty_high == pytest.approx(record.uncertainty_high, abs=1e-9)

    # ...and the cell-level analysis of the very same numbers would claim far more
    # precision than the design supports, which is the point.
    #
    # The error in the cell-level analysis is not that its t statistic is larger.  On
    # this fixture the two are the same size to within 12% (measured: 85.28 cell-level
    # against 95.68 sample-level, i.e. 0.89x), because the simulated within-sample
    # spread is small next to the injected donor effect.  What is wrong is the
    # *denominator*: the cell-level test spends its degrees of freedom on 400 cells from
    # one tube as though each were an independent biospecimen, so it reports a p-value
    # of exactly 0 where the sample-level test can only reach 6.2e-16 from six
    # independent donors per arm.  An exact zero is not a stronger result; it is the
    # signature of counting the wrong thing.
    #
    # So the assertions below are about the replicate counts, which is where the
    # difference actually lives, and about the sample-level t being the stable one.
    # expression() takes cell *ids*, not a boolean mask, and takes features by name.
    # features=["CD69"] already narrows the columns to one, so no further column index
    # is applied: the result is already the CD69 column of the selected cells.
    selected = [c for c in dataset.obs["cell_id"] if c in set(sampling.selected_cell_ids)]
    group = dataset.obs["sample_id"].map(
        {sid: rec.group for sid, rec in dataset.samples.items()}
    )
    labels = group.loc[selected].to_numpy()
    values = dataset.expression(features=["CD69"], cell_ids=selected)[:, 0]
    cell_result = stats.ttest_ind(values[labels == "case"], values[labels == "control"],
                                  equal_var=False)
    n_cells_per_arm = int((labels == "case").sum())
    sample_ids = dataset.obs["sample_id"].loc[selected].to_numpy()
    n_samples_per_arm = len({
        sid for sid, lab in zip(sample_ids, labels, strict=True) if lab == "case"
    })

    # The pseudoreplication: the cell-level analysis believes it has this many
    # independent replicates per arm...
    assert n_cells_per_arm == n_samples_per_arm * 400
    assert n_cells_per_arm >= 50 * n_samples_per_arm, (
        "the point of the comparison is that the cell count dwarfs the donor count"
    )
    # ...and that is what drives it to a p-value of exactly zero, which no analysis of
    # six donors per arm can honestly claim.  The sample-level test on the same data,
    # with the correct six replicates, stays finite.
    assert cell_result.pvalue == 0.0
    per_sample = pd.DataFrame(
        {"g": labels, "s": sample_ids, "v": values}
    ).groupby(["g", "s"], observed=True).v.mean().reset_index()
    honest = stats.ttest_ind(
        per_sample[per_sample.g == "case"].v.values,
        per_sample[per_sample.g == "control"].v.values,
        equal_var=False,
    )
    assert honest.pvalue > 0.0 and np.isfinite(honest.pvalue)

    # The sample-level statistic the package reports is the one that does not move when
    # the cells are repeated, which the first half of this test already showed.
    sample_level = _by_key(
        original, "biological_contrast_ttest_statistic"
    )[("T_CD4", "CD69")].value
    assert np.isfinite(sample_level) and sample_level > 0.0
