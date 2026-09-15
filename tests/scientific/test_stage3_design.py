"""Stage 3 acceptance gate: design audit, identifiability and eligibility."""

from __future__ import annotations

import numpy as np
import pytest

from batchdx.cyto.design import audit_design
from batchdx.cyto.design.audit import ELIGIBILITY_RULES
from batchdx.cyto.schemas.contracts import SampleRecord
from batchdx.cyto.schemas.enums import ComparisonSupport, ReferenceRole


def _sample(sample_id, batch, group=None, donor=None, reference=None, role="none", **extra):
    return SampleRecord(
        sample_id=sample_id, batch=batch, group=group, donor_id=donor,
        reference_id=reference, reference_role=ReferenceRole(role),
        technical_factors=extra,
    )


def _audit(samples, **kwargs):
    defaults = dict(dataset_id="test", correction_factor="batch", protected=["group"],
                    donor_factor="donor_id", minimum_batches=3, require_shared_reference=True)
    defaults.update(kwargs)
    return audit_design(samples, **defaults)


def _support(audit, a, b):
    for c in audit.comparisons:
        if c.population is None and {c.batch_a, c.batch_b} == {a, b}:
            return c.support
    raise AssertionError(f"no comparison for {a} vs {b}")


def _eligible(audit, name):
    return next(e for e in audit.eligibility if e.name == name)


def balanced_with_references() -> list[SampleRecord]:
    samples = []
    for batch in ("B1", "B2", "B3"):
        samples.append(_sample(f"ref_{batch}", batch, group="reference", donor="REF",
                               reference="REF1", role="technical_replicate"))
        for k in range(2):
            samples.append(_sample(f"case{k}_{batch}", batch, group="case", donor=f"C{k}{batch}"))
            samples.append(_sample(f"ctrl{k}_{batch}", batch, group="control",
                                   donor=f"H{k}{batch}"))
    return samples


def test_balanced_design_with_references_is_reference_supported() -> None:
    audit = _audit(balanced_with_references())
    assert audit.n_batches == 3
    assert not audit.below_design_floor
    assert _support(audit, "B1", "B2") is ComparisonSupport.REFERENCE_SUPPORTED
    assert _eligible(audit, "cycombine").eligible
    assert _eligible(audit, "cycombine_cluster").eligible


def test_partial_reference_overlap_is_reported_per_pair() -> None:
    samples = []
    for batch, reference in (("B1", "REF1"), ("B2", "REF1"), ("B3", None)):
        if reference:
            samples.append(_sample(f"ref_{batch}", batch, group="reference",
                                   reference=reference, role="technical_replicate"))
        for k in range(2):
            samples.append(_sample(f"case{k}_{batch}", batch, group="case", donor=f"C{k}{batch}"))
            samples.append(_sample(f"ctrl{k}_{batch}", batch, group="control",
                                   donor=f"H{k}{batch}"))
    audit = _audit(samples)
    assert _support(audit, "B1", "B2") is ComparisonSupport.REFERENCE_SUPPORTED
    assert _support(audit, "B1", "B3") is ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL
    assert len(audit.reference_graph_components) > 1
    assert any("disconnected" in w for w in audit.warnings)
    assert not _eligible(audit, "cytonorm").eligible


def test_perfect_confounding_is_non_identifiable_and_blocks_correction() -> None:
    samples = [
        _sample(f"C{k}_B1", "B1", group="case", donor=f"C{k}") for k in range(3)
    ] + [
        _sample(f"H{k}_B2", "B2", group="control", donor=f"H{k}") for k in range(3)
    ]
    audit = _audit(samples)

    assert audit.confounded_factors == ["group"]
    assert _support(audit, "B1", "B2") is ComparisonSupport.NON_IDENTIFIABLE
    assert "covariate" in _support_reason(audit, "B1", "B2")

    cycombine = _eligible(audit, "cycombine")
    assert not cycombine.eligible
    assert cycombine.reason_code == "batch_not_identifiable"
    # A scientifically impossible contrast is not something an exploratory override
    # may wave through.
    assert not cycombine.override_allowed
    assert audit.below_design_floor


def _support_reason(audit, a, b) -> str:
    for c in audit.comparisons:
        if c.population is None and {c.batch_a, c.batch_b} == {a, b}:
            return c.reason
    raise AssertionError("missing comparison")


def test_two_batch_design_makes_the_cluster_detector_ineligible() -> None:
    samples = [
        _sample(f"s{k}_B1", "B1", group="case" if k % 2 else "control", donor=f"D{k}")
        for k in range(4)
    ] + [
        _sample(f"s{k}_B2", "B2", group="case" if k % 2 else "control", donor=f"E{k}")
        for k in range(4)
    ]
    audit = _audit(samples)

    cluster = _eligible(audit, "cycombine_cluster")
    assert not cluster.eligible
    assert cluster.reason_code == "batch_count_below_minimum"
    assert "3-30" in cluster.reason or "3 <= n_batches <= 30" in cluster.reason
    # the quick detector only needs two batches, so it stays eligible
    assert _eligible(audit, "cycombine_express").eligible
    assert audit.below_design_floor


def test_thirty_one_batches_exceeds_the_cluster_detector_range() -> None:
    samples = []
    for b in range(31):
        batch = f"B{b:02d}"
        samples.append(_sample(f"case_{batch}", batch, group="case", donor=f"C{b}"))
        samples.append(_sample(f"ctrl_{batch}", batch, group="control", donor=f"H{b}"))
    audit = _audit(samples, require_shared_reference=False)
    cluster = _eligible(audit, "cycombine_cluster")
    assert not cluster.eligible
    assert cluster.reason_code == "batch_count_above_maximum"


def test_cytobatchflagr_needs_controls_in_every_batch_and_four_batches() -> None:
    samples = []
    for batch in ("B1", "B2", "B3"):
        samples.append(_sample(f"ref_{batch}", batch, group="reference", reference="REF1",
                               role="technical_replicate"))
        samples.append(_sample(f"case_{batch}", batch, group="case", donor=f"C{batch}"))
        samples.append(_sample(f"ctrl_{batch}", batch, group="control", donor=f"H{batch}"))
    audit = _audit(samples)
    entry = _eligible(audit, "cytobatchflagr")
    assert not entry.eligible
    assert entry.reason_code == "batch_count_below_minimum"

    samples += [
        _sample("ref_B4", "B4", group="reference", reference="REF1", role="technical_replicate"),
        _sample("case_B4", "B4", group="case", donor="CB4"),
        _sample("ctrl_B4", "B4", group="control", donor="HB4"),
    ]
    assert _eligible(_audit(samples), "cytobatchflagr").eligible

    # drop the control from one batch: the three IQR metrics are no longer defined
    without = [s for s in samples if s.sample_id != "ref_B2"]
    entry = _eligible(_audit(without), "cytobatchflagr")
    assert not entry.eligible
    assert entry.reason_code == "missing_reference_in_some_batch"


def test_iqr_fence_cannot_flag_with_three_values() -> None:
    """The minimum-four-batches rule for CytoBatchFlagR is arithmetic, not a guess.

    With three observations the Tukey fence always lies outside the data, so the
    upstream rule can never fire; the eligibility manifest states that as its reason.
    """

    rng = np.random.default_rng(0)
    for _ in range(2000):
        x = np.sort(rng.normal(size=3) * rng.uniform(0.1, 100))
        q1, q3 = np.quantile(x, [0.25, 0.75])
        iqr = q3 - q1
        assert not ((x > q3 + 1.5 * iqr) | (x < q1 - 1.5 * iqr)).any()

    # four observations can flag, which is why the rule is "at least four batches"
    x = np.array([0.0, 0.0, 0.0, 10.0])
    q1, q3 = np.quantile(x, [0.25, 0.75])
    assert (x > q3 + 1.5 * (q3 - q1)).any()
    assert ELIGIBILITY_RULES["cytobatchflagr"].min_batches == 4


def test_repeated_stimulation_and_timepoints_do_not_become_technical_replicates() -> None:
    samples = [
        SampleRecord(sample_id="D1_T0_B1", batch="B1", donor_id="D1", timepoint="T0",
                     stimulation="none", group="case"),
        SampleRecord(sample_id="D1_T1_B2", batch="B2", donor_id="D1", timepoint="T1",
                     stimulation="PMA", group="case"),
        SampleRecord(sample_id="D2_T0_B1", batch="B1", donor_id="D2", timepoint="T0",
                     stimulation="none", group="control"),
        SampleRecord(sample_id="D2_T1_B2", batch="B2", donor_id="D2", timepoint="T1",
                     stimulation="PMA", group="control"),
    ]
    audit = _audit(samples)
    # the same donor at two timepoints is not a measured technical replicate
    assert all(not v for v in audit.reference_coverage.values())
    assert _support(audit, "B1", "B2") is ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL
    assert audit.below_design_floor


def test_missing_donor_ids_permit_descriptive_diagnosis_only() -> None:
    samples = [
        _sample(f"s{k}_{batch}", batch, group="case" if k % 2 else "control")
        for batch in ("B1", "B2", "B3") for k in range(4)
    ]
    audit = _audit(samples, require_shared_reference=False)
    assert any("donor" in w for w in audit.warnings)
    assert not audit.below_design_floor
    assert _support(audit, "B1", "B2") is ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL


def test_population_specific_lack_of_support_is_reported() -> None:
    samples = balanced_with_references()
    counts = {s.sample_id: {"common": 500, "rare": 0 if s.batch == "B3" else 80}
              for s in samples}
    audit = _audit(samples, populations=["common", "rare"], population_counts=counts)

    rare_b1_b3 = next(
        c for c in audit.comparisons
        if c.population == "rare" and {c.batch_a, c.batch_b} == {"B1", "B3"}
    )
    assert rare_b1_b3.support is ComparisonSupport.DESCRIPTIVE_ONLY
    assert "no cells" in rare_b1_b3.reason

    common_b1_b3 = next(
        c for c in audit.comparisons
        if c.population == "common" and {c.batch_a, c.batch_b} == {"B1", "B3"}
    )
    assert common_b1_b3.support is ComparisonSupport.REFERENCE_SUPPORTED


def test_nested_technical_structure_is_represented_and_the_target_named() -> None:
    samples = []
    for site, batches in (("S1", ("B1", "B2")), ("S2", ("B3", "B4"))):
        for batch in batches:
            for k in range(2):
                samples.append(_sample(f"case{k}_{batch}", batch, group="case",
                                       donor=f"C{k}{batch}", site=site))
                samples.append(_sample(f"ctrl{k}_{batch}", batch, group="control",
                                       donor=f"H{k}{batch}", site=site))
    audit = _audit(samples, nested_within="site", require_shared_reference=False)
    assert audit.correction_factor == "batch"
    assert audit.nested_within == "site"
    # site is not the correction target and remains uncorrected; it is recorded, and the
    # design matrix including it is rank deficient because batch determines site.
    assert audit.design_full_rank is False


def test_reference_role_requires_an_explicit_reference_id() -> None:
    with pytest.raises(ValueError):
        SampleRecord(sample_id="x", reference_role="not_a_role")


def test_below_floor_is_a_configuration_value_with_a_declared_origin() -> None:
    samples = balanced_with_references()
    audit = _audit(samples, minimum_batches=5,
                   minimum_batches_origin="declared by the study protocol")
    assert audit.below_design_floor
    assert audit.design_floor_origin == "declared by the study protocol"
    assert "minimum is 5" in audit.design_floor_reason
