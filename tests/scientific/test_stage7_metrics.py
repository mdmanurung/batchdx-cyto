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


# --------------------------------------------------------------------------------------
# native namespace: the two specs are executed by cyCombine itself
# --------------------------------------------------------------------------------------


def _native_context(frozen, tmp_path, **kwargs):
    from batchdx.cyto.adapters import BackendContext

    config, _, design, sampling = frozen
    return BackendContext(
        design=design, sampling=sampling, work_dir=str(tmp_path),
        seed=config.project.seed, **kwargs
    )


@pytest.mark.requires_cycombine
def test_native_metrics_are_produced_by_calling_cycombine_not_by_reimplementing_it(
    frozen, tmp_path
) -> None:
    """The two `native` specs must run the upstream function, not a local copy.

    This is the clause Stage 7 could not previously satisfy: the specs were registered
    but had no producer, so nothing executed `cyCombine::evaluate_emd` or
    `evaluate_mad`.  They are called here on the identity "correction" (uncorrected
    against itself), where the honest answer is a reduction of zero, and the returned
    record must carry that zero as a real value from upstream rather than as a locally
    computed stand-in.
    """

    from batchdx.cyto.adapters import get_backend
    from batchdx.cyto.adapters.r_runner import r_package_version
    from batchdx.cyto.metrics.native import compute_native_metrics

    if not r_package_version("cyCombine").available:
        pytest.skip("cyCombine is not installed; see docs/blocked.md")

    config, dataset, design, sampling = frozen
    baseline = dataset.expression()
    corrected, _ = get_backend("none").run(dataset, "uncorrected", {}, _native_context(
        frozen, tmp_path
    ))

    records = compute_native_metrics(
        dataset, baseline, corrected, sampling, design,
        run_id="t", candidate_id="uncorrected",
        context=_native_context(frozen, tmp_path),
        work_dir=str(tmp_path / "native"),
    )
    by_name = {r.metric_name: r for r in records}
    assert set(by_name) == {"cycombine_evaluate_emd", "cycombine_evaluate_mad"}

    for name, record in by_name.items():
        assert record.metric_namespace == "native", name
        assert record.status is MetricStatus.OK, (
            f"{name}: {record.status} — {record.reason}"
        )
        assert record.value is not None, name
        # Identical inputs in, so upstream's own reduction is zero.  A tolerance is
        # used rather than ==0 because upstream sums and rounds internally.
        assert record.value == pytest.approx(0.0, abs=1e-9), (
            f"{name}: the identity correction cannot move an upstream reduction"
        )
        assert "cyCombine" in record.reason and "NOT comparable" in record.reason
        assert record.n_cells == len(sampling.selected_cell_ids)
        assert record.representation is not None


@pytest.mark.requires_cycombine
def test_a_native_metric_that_cannot_run_reports_a_reason_and_no_value(
    frozen, tmp_path, monkeypatch
) -> None:
    """An unavailable native metric is a status, never a zero score.

    A zero here would be indistinguishable from a genuine finding of 'no batch effect',
    which is the specific failure mode the registry's `native` namespace exists to
    prevent.  Forcing the R call to fail is the only honest way to see the branch.
    """

    from batchdx.cyto.adapters.r_runner import RScriptFailed, r_package_version
    from batchdx.cyto.metrics import native

    if not r_package_version("cyCombine").available:
        pytest.skip("cyCombine is not installed; see docs/blocked.md")

    config, dataset, design, sampling = frozen
    baseline = dataset.expression()

    def boom(*args, **kwargs):
        raise RScriptFailed("cycombine_evaluate.R", 1, "synthetic failure for the test")

    monkeypatch.setattr(native, "run_r_script", boom)
    records = native.compute_native_metrics(
        dataset, baseline, baseline, sampling, design,
        run_id="t", candidate_id="uncorrected",
        context=_native_context(frozen, tmp_path),
        work_dir=str(tmp_path / "native_fail"),
    )
    for record in records:
        assert record.status is MetricStatus.FAILED, record.metric_name
        assert record.value is None, (
            f"{record.metric_name}: a failed evaluation must not carry a value"
        )
        assert "synthetic failure" in record.reason

    # An absent interpreter is not a failure of the method; it is not_implemented here
    # and must be distinguishable from a broken evaluation.
    def absent(*args, **kwargs):
        raise native.RNotAvailable("no Rscript on PATH")

    monkeypatch.setattr(native, "run_r_script", absent)
    records = native.compute_native_metrics(
        dataset, baseline, baseline, sampling, design,
        run_id="t", candidate_id="uncorrected",
        context=_native_context(frozen, tmp_path),
        work_dir=str(tmp_path / "native_absent"),
    )
    statuses = {r.metric_name: r.status for r in records}
    assert set(statuses.values()) == {MetricStatus.NOT_IMPLEMENTED}
    assert all(r.value is None for r in records)


@pytest.mark.requires_cycombine
def test_the_native_r_script_refuses_a_misaligned_pair_of_frames(frozen, tmp_path) -> None:
    """Event identity is checked before the two frames are compared.

    `evaluate_emd` aligns on the `id` column.  If the corrected frame carried different
    events under the same labels, upstream would still return a plausible reduction,
    and the number would be meaningless.  The check is therefore in the R script, and
    this exercises the script rather than the Python wrapper.
    """

    import json
    import subprocess

    from batchdx.cyto.io.binary import write_matrix
    from batchdx.cyto.paths import workflow_root

    script = workflow_root() / "scripts" / "cycombine_evaluate.R"
    if not script.exists():
        pytest.skip("the packaged R script is not reachable from here")

    work = tmp_path / "misalign"
    work.mkdir(parents=True, exist_ok=True)
    markers = ["CD3", "CD4"]
    unc = np.arange(8, dtype=float).reshape(4, 2)
    cor = np.arange(8, dtype=float).reshape(4, 2)
    # Same shape, different events.
    write_matrix(work / "u.binmat", unc, row_ids=[f"c{i}" for i in range(4)],
                 col_names=markers)
    write_matrix(work / "c.binmat", cor, row_ids=[f"d{i}" for i in range(4)],
                 col_names=markers)
    params = {
        "uncorrected_matrix": str(work / "u.binmat"),
        "corrected_matrix": str(work / "c.binmat"),
        "markers": markers,
        "batch": ["B1", "B1", "B2", "B2"],
        "sample": ["s1", "s1", "s2", "s2"],
        "label": ["a", "a", "b", "b"],
        "cell_col": "label",
        "bin_size": 0.1,
    }
    (work / "params.json").write_text(json.dumps(params), encoding="utf-8")
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(script), str(work / "params.json"),
         str(work / "result.json")],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode != 0, "a misaligned pair must not be accepted"
    assert "same event ids" in proc.stdout + proc.stderr
