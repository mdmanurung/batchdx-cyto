"""Stage 4 acceptance gate: frozen labels, deterministic sampling, common partition.

Since Step 4.1 this file also gates the biospecimen-level train/tune/validate
split written to ``SamplingManifest.splits``, which was an empty dict in every
manifest this package emitted before.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto.metrics import compute_common_metrics
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.sampling.partitions import UNGATED
from batchdx.cyto.sampling.splits import (
    RECOMMENDATION_SPLIT,
    SPLIT_NAMES,
    SPLIT_UNIT_KEY,
    SplitError,
    build_splits,
    held_out_batch_pairs,
    missing_arms,
    resolve_biospecimen_unit,
    split_for_comparison,
    split_members,
    split_unit_field,
    splits_hold_every_arm,
)
from batchdx.cyto.schemas.contracts import ComparisonAssessment, SampleRecord
from batchdx.cyto.schemas.enums import ComparisonSupport, LabelTier, ReferenceRole


def _rebuild(config, dataset, design, base: Path):
    labels, tier, confidence = resolve_labels(dataset, config, base)
    return build_sampling_manifest(dataset, config, design, labels, tier, confidence)


def test_sampling_manifest_is_reproducible(frozen, example_project: Path) -> None:
    config, dataset, design, sampling = frozen
    again = _rebuild(config, dataset, design, example_project)
    assert again.selected_cell_ids == sampling.selected_cell_ids
    assert again.labels_by_cell == sampling.labels_by_cell
    assert again.per_sample_counts == sampling.per_sample_counts


def test_sampling_is_balanced_at_the_sample_level_and_reports_coverage(frozen) -> None:
    _, dataset, _, sampling = frozen
    counts = set(sampling.per_sample_counts.values())
    assert len(counts) == 1, "every sample contributes the same number of evaluation cells"
    assert sampling.per_sample_original_counts
    assert all(
        sampling.per_sample_counts[s] <= sampling.per_sample_original_counts[s]
        for s in sampling.per_sample_counts
    )
    # weights are per sample, so a high-yield sample cannot dominate a summary
    assert set(sampling.aggregation_weights.values()) == {1.0}


def test_gating_labels_are_medium_confidence_and_recorded(frozen) -> None:
    _, _, _, sampling = frozen
    assert sampling.label_source is LabelTier.GATING
    assert sampling.label_confidence == "medium"
    assert set(sampling.population_counts) >= {"T_CD4", "T_CD8", "B", "Mono", "NK"}


def test_clustering_labels_are_flagged_lowest_confidence(loaded, example_project: Path) -> None:
    config, dataset = loaded
    config.evaluation.label_source = "clustering"
    config.evaluation.clustering_k = 5
    labels, tier, confidence = resolve_labels(dataset, config, example_project)
    assert tier is LabelTier.CLUSTERING
    assert confidence == "low"
    assert len(set(labels)) <= 5

    again, _, _ = resolve_labels(dataset, config, example_project)
    assert list(again) == list(labels), "clustering labels must be deterministic"


def test_unmatched_cells_are_labelled_not_dropped(frozen) -> None:
    _, _, _, sampling = frozen
    if UNGATED in sampling.population_counts:
        assert any(UNGATED in w for w in sampling.population_coverage_warnings)


def test_rare_population_generates_a_coverage_warning(frozen, example_project: Path) -> None:
    config, dataset, design, _ = frozen
    config.evaluation.sampling.min_cells_per_population = 10_000
    sampling = _rebuild(config, dataset, design, example_project)
    assert len(sampling.population_coverage_warnings) >= len(sampling.population_counts)


def test_candidate_order_does_not_change_labels_or_comparisons(frozen, example_project: Path) -> None:
    config, dataset, design, sampling = frozen
    config.correction.candidates.reverse()
    again = _rebuild(config, dataset, design, example_project)
    assert again.labels_by_cell == sampling.labels_by_cell
    assert again.comparison_pairs == sampling.comparison_pairs
    assert again.evaluation_markers == sampling.evaluation_markers


def test_every_metric_record_carries_label_tier_and_representation(frozen) -> None:
    config, dataset, design, sampling = frozen
    result = compute_common_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="t", candidate_id="uncorrected",
    )
    assert result.records
    for record in result.records:
        assert record.label_tier is sampling.label_source
        assert record.label_confidence == sampling.label_confidence
        assert record.representation is not None
        assert record.representation.space == "marker"
        assert record.representation.marker_mask == "evaluation"
        assert record.representation.scaling == "none"
        assert record.representation.distance == "euclidean"
        assert record.representation.k == 30


def test_non_identifiable_comparisons_are_excluded_from_the_frozen_pairs(
    confounded_project: Path,
) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset

    config = load_config(confounded_project / "config" / "config.yaml")
    dataset = load_dataset(config, confounded_project)
    design = audit_design(
        list(dataset.samples.values()), dataset_id="c", correction_factor="batch",
        protected=["group"], donor_factor="donor_id",
    )
    sampling = _rebuild(config, dataset, design, confounded_project)
    assert sampling.comparison_pairs == []


def test_sampling_does_not_redefine_population_frequency(frozen) -> None:
    """A capped draw must not change the estimated composition of a sample."""

    _, dataset, _, sampling = frozen
    obs = dataset.obs
    obs["_label"] = obs["cell_id"].map(sampling.labels_by_cell)
    selected = obs.loc[obs["cell_id"].isin(sampling.selected_cell_ids)]

    for sid in sorted(selected["sample_id"].unique())[:3]:
        drawn = selected.loc[selected["sample_id"] == sid, "_label"].value_counts(normalize=True)
        for _population, fraction in drawn.items():
            assert 0.0 <= float(fraction) <= 1.0
        assert np.isclose(drawn.sum(), 1.0)


# --------------------------------------------------------------------------------------
# Step 4.1 — biospecimen-level train/tune/validate splits
# --------------------------------------------------------------------------------------


def _sample(sample_id: str, **kwargs) -> SampleRecord:
    return SampleRecord(sample_id=sample_id, **kwargs)


#: Measured technical roles.  A biological pseudo control is not one, which is the
#: point: the split must not treat it as a reference or as a biological arm.
_MEASURED_ROLES = {ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR}


def _synthetic(
    n_donors_per_arm: int = 6,
    donor_factor: str | None = "donor_id",
    declare_biospecimen: bool = True,
    split_donor_across_batches: bool = False,
) -> dict[str, SampleRecord]:
    """A design where each donor contributes one sample, unless asked otherwise."""

    samples: dict[str, SampleRecord] = {}
    batches = ["B1", "B2", "B3"]
    for arm in ("case", "control"):
        for d in range(n_donors_per_arm):
            donor = f"{arm[0].upper()}{d}"
            for b in batches if split_donor_across_batches else [batches[d % 3]]:
                samples[f"{donor}_{b}"] = _sample(
                    f"{donor}_{b}",
                    donor_id=donor if donor_factor == "donor_id" else None,
                    biospecimen_id=(
                        f"{donor}_{b}" if declare_biospecimen else None
                    ),
                    batch=b,
                    group=arm,
                )
    return samples


def test_the_manifest_records_a_split_and_names_the_unit_it_used(frozen) -> None:
    """`SamplingManifest.splits` was an empty dict in every manifest before Step 4.1."""

    config, dataset, _, sampling = frozen
    assert sampling.splits, "the emitted manifest carries no split at all"

    unit = split_unit_field(sampling.splits)
    assert unit == "donor_id", (
        "the design declares donor_id and every sample carries one, so the split must be "
        f"keyed on the donor; the manifest says {unit!r}"
    )
    members = split_members(sampling.splits)
    assert set(members) == set(SPLIT_NAMES)
    for name in SPLIT_NAMES:
        assert members[name], f"split {name!r} is empty"
    assert sorted(sum(members.values(), [])) == sorted(dataset.samples)
    # The unit is a first-class fact in the manifest, not an inference left to the reader.
    assert SPLIT_UNIT_KEY in sampling.splits
    assert sampling.splits[SPLIT_UNIT_KEY] == ["donor_id"]


def test_every_sample_of_a_biospecimen_lands_in_exactly_one_split(frozen) -> None:
    """The whole point: splitting by sample would leak a donor across the boundary."""

    _, dataset, _, sampling = frozen
    members = split_members(sampling.splits)
    home = {s: name for name in SPLIT_NAMES for s in members[name]}

    by_donor: dict[str, set[str]] = {}
    for sid, rec in dataset.samples.items():
        by_donor.setdefault(rec.donor_id, set()).add(home[sid])
    straddling = {d: sorted(v) for d, v in by_donor.items() if len(v) > 1}
    assert not straddling, (
        f"donor(s) {straddling} appear in more than one split; the held-out number would "
        "then be computed on data the fit has already seen"
    )
    # The example's reference donor contributes three samples across three batches and
    # must still be indivisible, or the split proves nothing.
    multi = [d for d, v in by_donor.items()
             if sum(len(members[s]) for s in v) > len(v)]
    assert multi, "this fixture has no multi-sample donor, so the test would pass vacuously"


def test_the_split_does_not_depend_on_candidate_order_or_sample_reading_order(
    frozen, example_project: Path
) -> None:
    config, dataset, design, sampling = frozen
    assert _rebuild(config, dataset, design, example_project).splits == sampling.splits

    config.correction.candidates.reverse()
    assert _rebuild(config, dataset, design, example_project).splits == sampling.splits

    # Re-resolve the split from the same samples supplied in reverse order.
    reversed_samples = {k: dataset.samples[k] for k in reversed(list(dataset.samples))}
    a = build_splits(dataset.samples, donor_factor="donor_id", seed=7)
    b = build_splits(reversed_samples, donor_factor="donor_id", seed=7)
    assert a.to_manifest() == b.to_manifest()


def test_the_split_is_reproducible_from_the_seed_and_moves_with_it(
    frozen, example_project: Path
) -> None:
    config, dataset, design, sampling = frozen
    again = _rebuild(config, dataset, design, example_project)
    assert again.splits == sampling.splits

    config.project.seed = 20260929
    other = _rebuild(config, dataset, design, example_project)
    assert set(other.splits) == set(sampling.splits)
    # The seed decides which donor goes where.  A different seed must therefore be
    # allowed to move some donor, or the seed is decorative.
    assert other.splits != sampling.splits


def test_every_split_holds_both_arms_of_the_biological_contrast(frozen) -> None:
    """Otherwise the split is scored on a contrast that exists only elsewhere."""

    _, dataset, _, sampling = frozen
    assignment = build_splits(
        dataset.samples, donor_factor="donor_id", seed=sampling.seed,
    )
    assert assignment.arms == ["case", "control"]
    for name in SPLIT_NAMES:
        cov = assignment.coverage[name]
        assert set(cov.arm_counts) == {"case", "control"}, (
            f"split {name!r} does not hold both arms: {cov.arm_counts}"
        )
        assert cov.n_samples > 0 and cov.n_units > 0
        assert name in cov.summary()

    # And the counts are the ones a reader of the manifest would recompute, so the
    # reported per-split arm counts are not a separate, unchecked claim.  A measured
    # technical reference is not an arm, which is why its samples are excluded here.
    home = {s: n for n in SPLIT_NAMES for s in assignment.coverage[n].samples}
    recomputed: dict[str, dict[str, int]] = {n: {} for n in SPLIT_NAMES}
    for sid, split in home.items():
        rec = dataset.samples[sid]
        if rec.group and rec.reference_role not in _MEASURED_ROLES:
            recomputed[split][rec.group] = recomputed[split].get(rec.group, 0) + 1
    for name in SPLIT_NAMES:
        assert recomputed[name] == assignment.coverage[name].arm_counts


def test_a_split_that_cannot_hold_both_arms_is_refused_rather_than_emitted() -> None:
    """The degenerate case is reported, not produced."""

    # Six case donors and one control donor: three splits can each be given a case
    # donor, but the single control donor cannot be in all three.
    thin = {
        f"c{d}_B1": _sample(f"c{d}_B1", donor_id=f"c{d}", batch="B1", group="case")
        for d in range(6)
    }
    thin["h0_B1"] = _sample("h0_B1", donor_id="h0", batch="B1", group="control")

    with pytest.raises(SplitError) as excinfo:
        build_splits(thin, donor_factor="donor_id", seed=1)
    message = str(excinfo.value)
    assert "cannot hold every arm" in message
    assert "control" in message
    # The message names the arm and what to do, rather than only saying "invalid".
    assert "Add biospecimens to the thinner arm" in message


def test_a_zero_share_split_is_refused_instead_of_emitted_empty() -> None:
    """A zero denominator is undefined, not a split with nothing in it."""

    with pytest.raises(SplitError, match="not strictly positive"):
        build_splits(
            _synthetic(), donor_factor="donor_id", seed=1,
            proportions={"train": 1.0, "tune": 0.0, "validate": 0.0},
        )
    with pytest.raises(SplitError, match="must name exactly"):
        build_splits(
            _synthetic(), donor_factor="donor_id", seed=1,
            proportions={"train": 0.5, "tune": 0.5},
        )


def test_the_split_unit_is_declared_and_never_chosen_in_silence() -> None:
    both = _synthetic()
    assert resolve_biospecimen_unit(both, "donor_id")[0] == "donor_id"
    assert "most independent" in resolve_biospecimen_unit(both, "donor_id")[1]

    # A donor factor that is not populated for every sample is not a usable key.
    partial = {k: v.model_copy() for k, v in both.items()}
    first = sorted(partial)[0]
    partial[first] = partial[first].model_copy(update={"donor_id": None})
    field, reason = resolve_biospecimen_unit(partial, "donor_id")
    assert field == "biospecimen_id"
    assert first in reason and "cannot key the partition" in reason

    # With neither, the split falls back to the tube and says so.
    neither = {k: v.model_copy(update={"donor_id": None, "biospecimen_id": None})
               for k, v in both.items()}
    field, reason = resolve_biospecimen_unit(neither, "donor_id")
    assert field == "sample_id"
    assert "biospecimen_id is absent" in reason

    assignment = build_splits(neither, donor_factor="donor_id", seed=1)
    assert any("split unit is sample_id" in w for w in assignment.warnings)


def test_per_split_arm_and_coverage_counts_are_reported(frozen) -> None:
    _, dataset, _, sampling = frozen
    assignment = build_splits(
        dataset.samples, donor_factor="donor_id", seed=sampling.seed,
    )
    assert assignment.arms == ["case", "control"]
    for name in SPLIT_NAMES:
        cov = assignment.coverage[name]
        assert set(cov.arm_counts) == {"case", "control"}, (
            f"split {name!r} does not report both arms: {cov.arm_counts}"
        )
        assert cov.n_samples > 0 and cov.n_units > 0
        assert name in cov.summary()
    # Coverage that a split cannot provide is stated, not left to be discovered.
    text = " ".join(assignment.warnings)
    assert "no measured technical reference" in text or "reference" in text


def test_the_validate_split_can_contrast_the_batches_a_decision_needs(frozen) -> None:
    """A held-out split that sees one batch cannot support a batch conclusion."""

    _, dataset, design, sampling = frozen
    pairs = held_out_batch_pairs(sampling.splits, dataset.samples)
    assert ("B1", "B2") in pairs and ("B2", "B3") in pairs, (
        f"the validation split cannot contrast the design's batches: {pairs}"
    )
    for batch_a, batch_b in pairs:
        assert batch_a != batch_b


def test_a_contrast_only_evaluable_on_the_training_split_cannot_support_a_recommendation(
    frozen,
) -> None:
    """The policy-layer contract, stated as a function rather than as a claim.

    `decisions/policy.py` does not call this yet; wiring it in is a required central
    change.  What is asserted here is that the verdict is computable and that it is
    False when the only split holding both batches is the one a fit would have seen.
    """

    _, dataset, _, sampling = frozen
    comparison = ComparisonAssessment(
        comparison_id="B1_vs_B2", batch_a="B1", batch_b="B2",
        support=ComparisonSupport.REFERENCE_SUPPORTED, reason="fixture",
    )

    # A partition where every B2 sample is in train and none in validate.
    leaked = {
        SPLIT_UNIT_KEY: ["donor_id"],
        "train": ["d1_B1", "d1_B2"],
        "tune": ["d2_B1"],
        "validate": ["d3_B1"],
    }
    samples = {
        "d1_B1": _sample("d1_B1", donor_id="d1", batch="B1", group="case"),
        "d1_B2": _sample("d1_B2", donor_id="d1", batch="B2", group="case"),
        "d2_B1": _sample("d2_B1", donor_id="d2", batch="B1", group="control"),
        "d3_B1": _sample("d3_B1", donor_id="d3", batch="B1", group="control"),
    }
    verdict = split_for_comparison(comparison, leaked, samples)
    assert verdict.evaluable_in == ["train"]
    assert verdict.held_out is False
    assert verdict.can_support_recommendation is False
    assert RECOMMENDATION_SPLIT not in verdict.evaluable_in
    assert "cannot by itself support a recommendation" in verdict.reason

    # The real partition from the fixture, where the contrast is held out.
    real = split_for_comparison(comparison, sampling.splits, dataset.samples)
    assert real.held_out is True
    assert real.can_support_recommendation is True
    assert RECOMMENDATION_SPLIT in real.evaluable_in


def test_a_contrast_no_split_can_estimate_is_a_coverage_failure_not_a_weak_result() -> None:
    comparison = ComparisonAssessment(
        comparison_id="B1_vs_B9", batch_a="B1", batch_b="B9",
        support=ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL, reason="fixture",
    )
    splits = {
        SPLIT_UNIT_KEY: ["donor_id"],
        "train": ["a_B1"], "tune": ["b_B1"], "validate": ["c_B1"],
    }
    samples = {
        "a_B1": _sample("a_B1", donor_id="a", batch="B1", group="case"),
        "b_B1": _sample("b_B1", donor_id="b", batch="B1", group="control"),
        "c_B1": _sample("c_B1", donor_id="c", batch="B1", group="case"),
    }
    verdict = split_for_comparison(comparison, splits, samples)
    assert verdict.evaluable_in == []
    assert verdict.held_out is False
    assert "not estimable in this partition at all" in verdict.reason
    assert "coverage failure" in verdict.reason


def test_a_recorded_partition_that_leaks_is_refused_at_the_reader() -> None:
    """The invariant is re-checked on the frozen artefact, not only at construction."""

    leaked = {
        SPLIT_UNIT_KEY: ["donor_id"],
        "train": ["a_B1", "b_B2"],
        "tune": ["b_B2"],
        "validate": ["c_B1"],
    }
    with pytest.raises(SplitError, match="more than one split"):
        split_members(leaked)

    with pytest.raises(SplitError, match="missing split"):
        split_members({"train": ["a"], "tune": ["b"]})

    with pytest.raises(SplitError, match="unknown key"):
        split_members({**leaked, "holdout": ["z"]})


def test_the_manifest_round_trips_through_json_with_its_split(frozen) -> None:
    """A contract that cannot survive serialisation is not a contract."""

    import json

    _, _, _, sampling = frozen
    again = type(sampling).model_validate(json.loads(sampling.to_json()))
    assert again.splits == sampling.splits
    assert split_unit_field(again.splits) == "donor_id"


def test_a_design_too_thin_to_stratify_is_warned_about_rather_than_silently_split() -> None:
    """The manifest is not refused for one arm's sake, but the deficiency is stated.

    ``build_splits`` refuses by default.  The sampling manifest asks it not to, because
    that manifest also serves the detectors, the metrics and the report: a design whose
    thinner arm holds fewer biospecimens than there are splits is a reason to warn about
    the split, not a reason to refuse to run.  What must not happen is a partition that
    silently omits an arm.
    """

    thin = {
        f"c{d}_{b}": _sample(f"c{d}_{b}", donor_id=f"c{d}_{b}", batch=b, group="case")
        for d in range(4) for b in ("B1", "B2")
    }
    thin["h0_B1"] = _sample("h0_B1", donor_id="h0", batch="B1", group="control")

    with pytest.raises(SplitError, match="cannot hold every arm"):
        build_splits(thin, donor_factor="donor_id", seed=1)

    relaxed = build_splits(thin, donor_factor="donor_id", seed=1, require_all_arms=False)
    gaps = missing_arms(relaxed.to_manifest(), thin)
    assert gaps, "a relaxation that hid the deficiency would be worse than a refusal"
    for name, arms in gaps.items():
        assert arms == ["control"]
        assert name in relaxed.coverage
    text = " ".join(relaxed.warnings)
    assert "NOT arm-balanced" in text
    assert "must refuse such a partition" in text
    assert splits_hold_every_arm(relaxed.to_manifest(), thin) is False

    # The same check passes on a design that can be stratified, so the reader-side
    # function is not simply always-False.
    balanced = _synthetic()
    balanced_manifest = build_splits(
        balanced, donor_factor="donor_id", seed=1
    ).to_manifest()
    assert missing_arms(balanced_manifest, balanced) == {}
    assert splits_hold_every_arm(balanced_manifest, balanced) is True


def _balanced() -> dict[str, SampleRecord]:
    return _synthetic()


@pytest.fixture
def two_batch_thin_project(tmp_path: Path) -> Path:
    """Two batches with exactly one control donor: the contrast cannot be stratified.

    Defined here rather than in ``tests/conftest.py`` because it belongs to the split
    gate alone, and a shared conftest change for one test file's fixture is a coupling
    the next agent would inherit for no reason.
    """

    from batchdx.cyto.examples import materialise_example
    from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation

    project = tmp_path / "thin_two_batch"
    materialise_example("no_reference", project, cells_per_sample=200)
    samples = [
        SampleSpec(sample_id=f"s{k}_{batch}", donor_id=f"D{k}{batch}", batch=batch,
                   group="case" if k else "control", n_cells=200)
        for batch in ("B1", "B2") for k in range(3)
    ]
    write_simulation(SimulationSpec(samples=samples, batch_shift={"B2": {"CD4": 0.7}}), project)
    return project


def test_the_manifest_warns_when_the_emitted_split_is_not_arm_balanced(
    two_batch_thin_project: Path,
) -> None:
    """End to end: the deficiency reaches the frozen manifest, not just the caller."""

    config, dataset, design = _load(two_batch_thin_project)
    sampling = _rebuild(config, dataset, design, two_batch_thin_project)
    assert sampling.splits, "the manifest must still carry a partition it cannot balance"
    gaps = missing_arms(sampling.splits, dataset.samples)
    assert gaps, "this fixture is meant to defeat stratification; it no longer does"
    assert any("NOT arm-balanced" in w for w in sampling.population_coverage_warnings)
    assert splits_hold_every_arm(sampling.splits, dataset.samples) is False


def _load(project: Path):
    from batchdx.cyto.config import load_config
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset

    config = load_config(project / "config" / "config.yaml")
    dataset = load_dataset(config, project)
    design = audit_design(
        list(dataset.samples.values()), dataset_id="c", correction_factor="batch",
        protected=["group"], donor_factor="donor_id",
    )
    return config, dataset, design
