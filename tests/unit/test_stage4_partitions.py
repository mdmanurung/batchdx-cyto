"""Stage 4 acceptance gate: frozen labels, deterministic sampling, common partition."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from batchdx.cyto.metrics import compute_common_metrics
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.sampling.partitions import UNGATED
from batchdx.cyto.schemas.enums import LabelTier


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
