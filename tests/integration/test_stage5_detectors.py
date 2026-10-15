"""Stage 5 acceptance gate: detectors run headlessly, and skips are explicit.

The two native detection systems need R packages that could not be installed in the
environment this build was developed in (see ``docs/blocked.md``).  The tests that need
them skip with a stated reason; what they must *never* do is pass by producing flags
from an absent dependency, and that is asserted here directly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto.adapters.r_runner import r_package_version
from batchdx.cyto.design import audit_design
from batchdx.cyto.diagnostics import DetectorContext, get_detector
from batchdx.cyto.schemas.enums import MetricStatus
from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation


def _context(frozen, tmp_path, **kwargs):
    config, dataset, design, sampling = frozen
    return DetectorContext(
        run_id="t", design=design, sampling=sampling, out_dir=str(tmp_path),
        seed=config.project.seed, **kwargs,
    )


def test_drift_detector_finds_an_injected_within_file_trend(tmp_path: Path) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.examples import materialise_example
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels

    project = tmp_path / "drift"
    materialise_example("paired_reference", project, cells_per_sample=400)

    # Re-simulate with a monotonic within-file trend on CD8 in exactly one sample.
    from batchdx.cyto.simulate import paired_reference_spec

    spec = paired_reference_spec(n_cells=400)
    target = spec.samples[1].sample_id
    spec.drift = {target: {"CD8": 2.0}}
    write_simulation(spec, project)

    config = load_config(project / "config" / "config.yaml")
    dataset = load_dataset(config, project)
    design = audit_design(list(dataset.samples.values()), dataset_id="d",
                          correction_factor="batch", protected=["group"])
    labels, tier, confidence = resolve_labels(dataset, config, project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)

    result = get_detector("acquisition_drift").run(
        dataset,
        DetectorContext(run_id="t", design=design, sampling=sampling,
                        out_dir=str(tmp_path / "out"), seed=473),
    )
    assert result.status == "ok"

    flagged = {
        (a.sample_id, a.marker) for a in result.assessments
        if a.statistic == "spearman_rho_vs_acquisition_order" and a.flag
    }
    assert (target, "CD8") in flagged
    assert any("not addressed by batch-level correction" in w for w in result.warnings)

    # a sample with no injected drift must not be flagged for the same marker
    untouched = spec.samples[2].sample_id
    assert (untouched, "CD8") not in flagged


def test_drift_detector_is_ineligible_without_an_acquisition_order(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    dataset.adata.obs = dataset.adata.obs.drop(
        columns=["acquisition_order", "original_event_index"]
    )
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    assert result.status == "ineligible"
    assert "acquisition order" in result.reason
    assert not result.assessments


def test_composition_detector_separates_frequency_from_intensity(frozen, tmp_path) -> None:
    result = get_detector("composition").run(frozen[1], _context(frozen, tmp_path))
    assert result.status == "ok"

    frequencies = [a for a in result.assessments if a.statistic == "population_frequency"]
    assert frequencies
    assert all(a.evidence_class == "descriptive" for a in frequencies)

    tests = [a for a in result.assessments if a.statistic == "kruskal_batch_frequency_pvalue"]
    assert tests
    assert all(
        "not evidence that marker" in (a.reason or "")
        for a in tests if a.status is MetricStatus.OK
    )
    assert result.requirements["independent_unit"] == "sample"
    assert result.parameters["candidate_explanations"]


def test_composition_detector_reports_insufficient_data_rather_than_a_p_value(
    frozen, tmp_path
) -> None:
    """One sample per batch is not enough: the sample is the independent unit."""

    config, dataset, design, sampling = frozen
    keep = []
    seen = set()
    batch_of = dataset.sample_batches()
    for sid in sorted(dataset.samples):
        batch = batch_of[sid]
        if batch not in seen:
            seen.add(batch)
            keep.append(sid)
    obs = dataset.obs
    allowed = set(obs.loc[obs["sample_id"].isin(keep), "cell_id"])
    trimmed = sampling.model_copy(update={
        "labels_by_cell": {k: v for k, v in sampling.labels_by_cell.items() if k in allowed},
    })
    context = DetectorContext(run_id="t", design=design, sampling=trimmed,
                              out_dir=str(tmp_path), seed=473)
    result = get_detector("composition").run(dataset, context)
    tests = [a for a in result.assessments if a.statistic == "kruskal_batch_frequency_pvalue"]
    assert tests
    assert all(a.status is MetricStatus.INSUFFICIENT_DATA for a in tests)
    assert all(a.value is None for a in tests)


def test_two_batch_fixture_makes_the_cluster_detector_ineligible_not_uninformative(
    tmp_path: Path,
) -> None:
    samples = []
    for batch in ("B1", "B2"):
        for k in range(3):
            samples.append(SampleSpec(sample_id=f"s{k}_{batch}", donor_id=f"D{k}{batch}",
                                      batch=batch, group="case" if k else "control",
                                      n_cells=200))
    spec = SimulationSpec(samples=samples, batch_shift={"B2": {"CD4": 0.7}})
    project = tmp_path / "twobatch"
    from batchdx.cyto.examples import materialise_example

    materialise_example("no_reference", project, cells_per_sample=200)
    write_simulation(spec, project)

    from batchdx.cyto.config import load_config
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels

    config = load_config(project / "config" / "config.yaml")
    dataset = load_dataset(config, project)
    design = audit_design(list(dataset.samples.values()), dataset_id="two",
                          correction_factor="batch", protected=["group"])
    labels, tier, confidence = resolve_labels(dataset, config, project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    context = DetectorContext(run_id="t", design=design, sampling=sampling,
                              out_dir=str(tmp_path / "out"), seed=473)

    result = get_detector("cycombine_cluster").run(dataset, context)
    assert result.status == "ineligible"
    assert result.eligibility is not None
    assert result.eligibility.reason_code == "batch_count_below_minimum"
    assert "3-30" in result.reason or "3 <= n_batches <= 30" in result.reason
    assert not result.assessments, "an ineligible detector produces no flags at all"


def test_absent_r_package_yields_not_available_never_an_absence_of_batch_effects(
    frozen, tmp_path
) -> None:
    # Reproduce the unavailable case rather than waiting for an environment that lacks
    # cyCombine.  Gating on whether the package is installed meant this test was a
    # no-op in every environment where the native detector actually runs, which is how
    # the most important failure path ended up with no coverage at all.
    config, dataset, design, sampling = frozen
    # three batches, so the detector is eligible and the only obstacle is the package
    result = get_detector("cycombine_cluster").run(
        dataset, _context(frozen, tmp_path, rscript="Rscript-does-not-exist")
    )
    assert result.status == "not_available"
    assert "not installed" in result.reason
    assert "never reported as an absence of batch effects" in result.reason
    assert not result.assessments


def test_cytobatchflagr_without_a_checkout_is_not_available_with_a_stated_remedy(
    frozen, tmp_path
) -> None:
    result = get_detector("cytobatchflagr").run(
        frozen[1], _context(frozen, tmp_path, upstream_dir=str(tmp_path / "nowhere"))
    )
    assert result.status in ("not_available", "ineligible")
    assert not result.assessments
    if result.status == "not_available":
        assert "fetch_upstream.sh" in result.reason
        assert "never a fabricated flag" in result.reason


def test_cytobatchflagr_refuses_to_choose_thresholds_itself(frozen, tmp_path, repo_root) -> None:
    """Analyst-to-analyst threshold variation is the documented failure mode."""

    config, dataset, design, sampling = frozen
    # four batches with controls everywhere would make it eligible; here we only check
    # that the missing-threshold path is a stated limitation rather than a guess.
    context = _context(frozen, tmp_path, upstream_dir=str(repo_root / ".upstream"),
                       threshold_path=None)
    result = get_detector("cytobatchflagr").run(dataset, context)
    assert result.status in ("ineligible", "not_available")
    assert not result.assessments


@pytest.mark.requires_cycombine
def test_native_cycombine_detector_matches_a_direct_upstream_call(frozen, tmp_path) -> None:
    if not r_package_version("cyCombine").available:
        pytest.skip(
            "cyCombine is not installed; the native detector comparison cannot run. "
            "See docs/blocked.md."
        )
    config, dataset, design, sampling = frozen
    result = get_detector("cycombine_express").run(dataset, _context(frozen, tmp_path))
    assert result.status == "ok"
    assert any(a.statistic == "cycombine_emd" for a in result.assessments)


def test_detector_results_round_trip_through_json(frozen, tmp_path) -> None:
    from batchdx.cyto.schemas.contracts import DetectorResult

    result = get_detector("composition").run(frozen[1], _context(frozen, tmp_path))
    path = tmp_path / "detector.json"
    result.write(path)
    again = DetectorResult.read(path)
    assert len(again.assessments) == len(result.assessments)
    assert again.status == result.status


def test_flag_counts_are_not_presented_as_probabilities(frozen, tmp_path) -> None:
    result = get_detector("composition").run(frozen[1], _context(frozen, tmp_path))
    for assessment in result.assessments:
        if assessment.flag is not None:
            assert assessment.threshold is not None, (
                "a flag must carry the threshold that produced it"
            )
        if assessment.value is not None:
            assert np.isfinite(assessment.value)


def test_drift_detector_reports_undefined_never_zero_for_a_constant_marker(
    frozen, tmp_path
) -> None:
    """Regression: a constant marker has no defined rank trend.  The detector used to
    coerce the NaN Spearman to 0.0 with status OK — a missing value recorded as a
    measured zero, violating the package's own contract."""
    config, dataset, design, sampling = frozen

    # Force one marker to be constant within one sample of the baseline.
    marker = "CD8"
    sid = sorted(dataset.samples)[0]
    obs = dataset.adata.obs
    mask = (obs["sample_id"] == sid).to_numpy()
    layer = dataset.adata.layers["preprocessed"]
    col = list(dataset.adata.var_names).index(marker)
    layer[mask, col] = 1.0

    result = get_detector("acquisition_drift").run(
        dataset, _context(frozen, tmp_path / "const")
    )
    assert result.status == "ok"
    target = [
        a for a in result.assessments
        if a.sample_id == sid and a.marker == marker
        and a.statistic == "spearman_rho_vs_acquisition_order"
    ]
    assert target, "no drift assessment for the constant sample-marker pair"
    assert all(a.status is MetricStatus.UNDEFINED and a.value is None for a in target)
    assert not any(a.flag for a in target)
