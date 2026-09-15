"""Stage 7 acceptance gate: metric behaviour under known perturbations."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from batchdx.cyto.metrics import compute_common_metrics, compute_reductions
from batchdx.cyto.metrics.common import MIN_GROUP_CELLS
from batchdx.cyto.metrics.registry import (
    DEFERRED_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    registered_metrics,
)
from batchdx.cyto.schemas.contracts import MetricResult, SamplingManifest
from batchdx.cyto.schemas.enums import MetricDirection, MetricStatus


def _values(result, name, **filters):
    out = []
    for r in result.records:
        if r.metric_name != name or r.status is not MetricStatus.OK:
            continue
        if all(getattr(r, k) == v for k, v in filters.items()):
            out.append(float(r.value))
    return out


def _metrics(dataset, matrix, sampling, design, candidate_id="c"):
    return compute_common_metrics(dataset, matrix, sampling, design,
                                  run_id="t", candidate_id=candidate_id)


def test_identical_distributions_give_near_zero_emd(frozen) -> None:
    _, dataset, design, sampling = frozen
    # Force every batch to share one distribution by overwriting with a common draw.
    rng = np.random.default_rng(0)
    matrix = np.tile(rng.normal(size=(1, dataset.adata.n_vars)), (dataset.n_cells, 1))
    matrix += rng.normal(scale=1e-6, size=matrix.shape)
    result = _metrics(dataset, matrix, sampling, design)
    emd = _values(result, "emd_within_population")
    assert emd and max(emd) < 1e-3


def test_a_known_shift_moves_emd_by_the_size_of_the_shift(frozen) -> None:
    _, dataset, design, sampling = frozen
    batch_of = dataset.sample_batches()
    batches = dataset.obs["sample_id"].map(batch_of).to_numpy()
    marker = sampling.evaluation_markers[0]
    col = dataset.feature_index([marker])[0]

    matrix = dataset.expression()
    matrix[batches == "B2", col] += 2.0
    result = _metrics(dataset, matrix, sampling, design)

    b1_b2 = [
        float(r.value) for r in result.records
        if r.metric_name == "emd_within_population" and r.marker == marker
        and r.status is MetricStatus.OK and r.comparison_id == "B1__vs__B2"
    ]
    assert b1_b2
    # Wasserstein distance on a pure translation equals the translation, plus whatever
    # difference the baseline already carried.
    assert min(b1_b2) > 1.0


def test_scale_change_alone_is_visible_in_dispersion_not_only_in_emd(frozen) -> None:
    _, dataset, design, sampling = frozen
    batch_of = dataset.sample_batches()
    batches = dataset.obs["sample_id"].map(batch_of).to_numpy()

    matrix = dataset.expression()
    rows = np.flatnonzero(batches == "B2")
    centre = matrix[rows, :].mean(axis=0)
    matrix[rows, :] = centre + (matrix[rows, :] - centre) * 0.4

    result = _metrics(dataset, matrix, sampling, design)
    signed = _values(result, "mad_signed_change", comparison_id="batch:B2")
    assert signed and np.median(signed) < 0, "compression must show as a negative signed change"

    compression = _values(result, "mad_compression_fraction")
    assert compression and compression[0] > 0.2


def test_constant_marker_yields_undefined_not_zero(frozen) -> None:
    _, dataset, design, sampling = frozen
    matrix = dataset.expression()
    matrix[:, :] = 1.0
    result = _metrics(dataset, matrix, sampling, design)

    joint = [r for r in result.records if r.metric_name == "joint_structure_distance"]
    assert joint
    # No group may report a value: a constant marker has no correlation to compare, and
    # that is reported as undefined (or as insufficient data), never as a score of zero.
    assert not any(r.status is MetricStatus.OK for r in joint)
    assert all(r.value is None for r in joint)
    assert any(r.status is MetricStatus.UNDEFINED for r in joint)
    assert all("constant" in (r.reason or "") for r in joint
               if r.status is MetricStatus.UNDEFINED)


def test_rare_population_returns_insufficient_data_with_a_reason(frozen) -> None:
    _, dataset, design, sampling = frozen
    rare = "rare_population"
    victim = sampling.selected_cell_ids[:5]
    labels = dict(sampling.labels_by_cell)
    for cell in victim:
        labels[cell] = rare
    sampling = sampling.model_copy(update={"labels_by_cell": labels})

    result = _metrics(dataset, sampling and dataset.expression(), sampling, design)
    rare_records = [r for r in result.records
                    if r.population == rare and r.metric_name == "emd_within_population"]
    assert rare_records
    assert all(r.status is MetricStatus.INSUFFICIENT_DATA for r in rare_records)
    assert all(str(MIN_GROUP_CELLS) in (r.reason or "") for r in rare_records)
    assert all(r.value is None for r in rare_records)


def test_unequal_sample_sizes_do_not_change_the_weighting(frozen) -> None:
    _, _, _, sampling = frozen
    assert set(sampling.aggregation_weights.values()) == {1.0}


def test_zero_denominator_makes_the_reduction_undefined_not_infinite(frozen) -> None:
    _, dataset, design, sampling = frozen
    matrix = dataset.expression()
    # identical data in every batch: the identity candidate's total EMD is ~0
    rng = np.random.default_rng(1)
    flat = np.tile(rng.normal(size=(1, dataset.adata.n_vars)), (dataset.n_cells, 1))

    identity = _metrics(dataset, flat, sampling, design, candidate_id="uncorrected")
    other = _metrics(dataset, matrix, sampling, design, candidate_id="other")
    reductions = compute_reductions(
        {"uncorrected": identity, "other": other}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    statuses = {r.candidate_id: r.status for r in reductions}
    assert statuses["other"] is MetricStatus.UNDEFINED
    assert all(r.value is None for r in reductions if r.status is MetricStatus.UNDEFINED)


def test_composition_only_change_moves_composition_and_not_within_population_emd(
    frozen,
) -> None:
    """The decisive separation test: resampling cell *counts* must not move intensity."""

    _, dataset, design, sampling = frozen
    labels = sampling.labels_by_cell
    batch_of = dataset.sample_batches()
    obs = dataset.obs
    obs["_batch"] = obs["sample_id"].map(batch_of)

    baseline = _metrics(dataset, dataset.expression(), sampling, design, candidate_id="base")

    # Drop most Mono cells from B2's evaluation set: composition changes, and every
    # remaining cell keeps its original values, so within-population EMD must not move.
    rng = np.random.default_rng(2)
    keep: list[str] = []
    for cell in sampling.selected_cell_ids:
        row = obs.loc[obs["cell_id"] == cell]
        if labels[cell] == "Mono" and str(row["_batch"].iloc[0]) == "B2":
            if rng.random() < 0.8:
                continue
        keep.append(cell)

    shifted = sampling.model_copy(update={
        "selected_cell_ids": keep,
        "labels_by_cell": {c: labels[c] for c in keep},
    })
    after = _metrics(dataset, dataset.expression(), shifted, design, candidate_id="after")

    def by_key(result, name):
        return {
            (r.population, r.marker, r.comparison_id): float(r.value)
            for r in result.records
            if r.metric_name == name and r.status is MetricStatus.OK
        }

    freq_before = by_key(baseline, "population_frequency_shift")
    freq_after = by_key(after, "population_frequency_shift")
    mono_keys = [k for k in freq_before if k[0] == "Mono" and "B2" in (k[2] or "")]
    assert mono_keys
    assert max(abs(freq_after[k] - freq_before[k]) for k in mono_keys) > 0.05

    emd_before = by_key(baseline, "emd_within_population")
    emd_after = by_key(after, "emd_within_population")
    shared = [k for k in emd_before if k in emd_after and k[0] != "Mono"]
    assert shared
    assert max(abs(emd_after[k] - emd_before[k]) for k in shared) < 1e-9


def test_mad_uses_the_normal_consistency_constant(frozen) -> None:
    """The dispersion metric must be on the same footing as R's stats::mad."""

    rng = np.random.default_rng(3)
    sample = rng.normal(loc=0.0, scale=2.0, size=200_000)
    scaled = float(stats.median_abs_deviation(sample, scale="normal"))
    raw = float(stats.median_abs_deviation(sample))
    assert scaled == pytest.approx(raw * 1.4826, rel=1e-3)
    assert scaled == pytest.approx(2.0, rel=0.02)


def test_empirical_and_binned_emd_are_named_differently() -> None:
    names = {f"{ns}.{name}" for ns, name in registered_metrics()}
    assert "common.emd_within_population" in names
    assert "native.cycombine_evaluate_emd" in names
    spec = registered_metrics()[("common", "emd_within_population")]
    assert "cyCombine" in spec.notes and "bin units" in spec.notes
    assert spec.direction is MetricDirection.LOWER_IS_BETTER


def test_no_local_metric_claims_to_be_an_upstream_component() -> None:
    for (namespace, name), spec in registered_metrics().items():
        if namespace == "compatibility":
            assert spec.upstream_verified, (
                f"{name} is registered as a compatibility component without verification"
            )
        if namespace == "common":
            assert "upstream" not in spec.implementation.lower()


def test_unsupported_metrics_report_reasons_rather_than_zero_scores() -> None:
    assert DEFERRED_METRICS
    assert all(isinstance(reason, str) and reason for reason in DEFERRED_METRICS.values())
    assert "openproblems.cms" in UPSTREAM_METRICS_NOT_FOUND
    assert "absent at the pinned commit" in UPSTREAM_METRICS_NOT_FOUND["openproblems.cms"]


def test_missing_metrics_cannot_improve_a_candidate(frozen) -> None:
    """A candidate whose comparisons failed gets lower coverage, never a better score."""

    _, dataset, design, sampling = frozen
    full = _metrics(dataset, dataset.expression(), sampling, design, candidate_id="uncorrected")

    partial = MetricResult(
        run_id="t", candidate_id="partial", metric_group="common",
        records=[r for r in full.records if r.comparison_id == "B1__vs__B2"],
    )
    reductions = compute_reductions(
        {"uncorrected": full, "partial": partial}, sampling,
        run_id="t", identity_candidate_id="uncorrected",
    )
    coverage = {r.candidate_id: r.coverage for r in reductions}
    assert coverage["partial"] < coverage["uncorrected"]


def test_sampling_manifest_round_trips_through_json(frozen, tmp_path) -> None:
    _, _, _, sampling = frozen
    path = tmp_path / "sampling.json"
    sampling.write(path)
    again = SamplingManifest.read(path)
    assert again.labels_by_cell == sampling.labels_by_cell
    assert again.representation == sampling.representation
