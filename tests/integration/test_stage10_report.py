"""Stage 10 acceptance gate: offline report, traceable scores, safe exports."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto.io.export import export_fcs_per_sample
from batchdx.cyto.io.fcs import read_fcs
from batchdx.cyto.pipeline import (
    RunLayout,
    stage_correct,
    stage_decide,
    stage_design,
    stage_diagnose,
    stage_metrics,
    stage_partition,
    stage_report,
    stage_status,
    stage_validate,
)

EXTERNAL_URL = re.compile(r'(?:src|href)\s*=\s*["\'](https?:)?//', re.I)


@pytest.fixture
def completed_run(example_project: Path, tmp_path: Path):
    layout = RunLayout(tmp_path / "run")
    config_path = example_project / "config" / "config.yaml"
    stage_validate(config_path, example_project, layout)
    stage_design(config_path, example_project, layout)
    stage_partition(config_path, example_project, layout)
    for detector in ("acquisition_drift", "composition"):
        stage_diagnose(detector, config_path, example_project, layout)
    for candidate in ("uncorrected", "cycombine_default", "cytonorm_reference"):
        stage_correct(candidate, config_path, example_project, layout)
        stage_metrics(candidate, config_path, example_project, layout)
    stage_decide(config_path, example_project, layout)
    stage_status(config_path, layout)
    stage_report(config_path, example_project, layout)
    return config_path, example_project, layout


def test_report_is_self_contained_and_offline(completed_run) -> None:
    _, _, layout = completed_run
    html = layout.report.read_text(encoding="utf-8")
    assert html.startswith("<!doctype html>")
    external = EXTERNAL_URL.findall(html)
    assert not external, f"the report references external assets: {external[:3]}"
    assert "<svg" in html, "figures are embedded as inline SVG"
    assert "data:text/csv;base64," in html, "metric records are downloadable from the page"
    assert 'id="metricFilter"' in html, "the metric table is filterable"


def test_every_headline_number_is_traceable_to_a_record(completed_run) -> None:
    _, _, layout = completed_run
    decision = json.loads(layout.decision.read_text())
    html = layout.report.read_text(encoding="utf-8")

    for assessment in decision["assessments"]:
        assert assessment["candidate_id"] in html
        assert assessment["state"] in html

    # the metric table in the page and the records on disk agree on how many rows exist
    records = 0
    for candidate in ("uncorrected", "cycombine_default", "cytonorm_reference"):
        path = layout.metric_result(candidate)
        records += len(json.loads(path.read_text())["records"])
    assert records > 0
    assert html.count("emd_within_population") >= 2

    # The narrative is built from the records, so its numbers are traceable too: the
    # metric-record count on the page is the count in the artifacts -- the per-candidate
    # results plus the reduction records, which is exactly what the page's table renders.
    page_records = records + sum(
        len(json.loads(path.read_text())["records"])
        for path in [layout.reductions] if path.exists()
    )
    assert f"{page_records} metric record" in html
    detector_assessments = sum(
        len(json.loads(path.read_text())["assessments"])
        for path in layout.root.glob("diagnostics/*/detector_result.json")
    )
    assert f"The detectors produced {detector_assessments} assessments" in html

    # The embedding panels state what they project: a deterministic PCA over the
    # evaluation markers, one point per sample, and whether a UMAP exists.
    assert "Sample embeddings" in html
    assert "deterministic PCA (full SVD)" in html
    assert "umap package is not installed" in html or "UMAP panel is drawn" in html


def test_failed_and_unavailable_candidates_stay_visible(completed_run) -> None:
    _, _, layout = completed_run
    html = layout.report.read_text(encoding="utf-8")
    assert "Failures and limitations" in html
    assert "cytonorm_reference" in html
    assert "not_available" in html or "insufficient_evidence" in html
    # metrics named in the plan but not produced are listed with a reason
    assert "openproblems.cms" in html
    assert "absent at the pinned commit" in html


def test_missing_preservation_evidence_is_not_presented_as_preservation(completed_run) -> None:
    _, _, layout = completed_run
    html = layout.report.read_text(encoding="utf-8")
    assert "Missing preservation evidence is not evidence of preservation" in html
    assert "not implemented" in html


def test_identifiers_are_redacted_and_the_mapping_is_stored_separately(completed_run) -> None:
    _, project, layout = completed_run
    html = layout.report.read_text(encoding="utf-8")

    import pandas as pd

    samples = pd.read_csv(project / "config" / "samples.tsv", sep="\t")
    for donor in samples["donor_id"].astype(str):
        assert donor not in html, f"donor {donor} leaked into the report"

    mapping = layout.report.parent / "identifier_map.json"
    assert mapping.exists()
    payload = json.loads(mapping.read_text())
    assert payload["mapping"]
    assert "separate" in payload["note"]


def test_regeneration_triggers_no_correction(completed_run) -> None:
    config_path, project, layout = completed_run
    corrected = layout.corrected_h5ad("uncorrected")
    before = corrected.stat().st_mtime_ns
    result_before = layout.correction_result("uncorrected").read_text()

    stage_report(config_path, project, layout)

    assert corrected.stat().st_mtime_ns == before
    assert layout.correction_result("uncorrected").read_text() == result_before


def test_diagnosis_only_report_renders_before_any_correction(
    example_project: Path, tmp_path: Path
) -> None:
    layout = RunLayout(tmp_path / "diagnosis")
    config_path = example_project / "config" / "config.yaml"
    stage_validate(config_path, example_project, layout)
    stage_design(config_path, example_project, layout)
    stage_partition(config_path, example_project, layout)
    for detector in ("acquisition_drift", "composition"):
        stage_diagnose(detector, config_path, example_project, layout)

    stage_report(config_path, example_project, layout, diagnosis_only=True)
    html = layout.report.read_text(encoding="utf-8")
    assert "diagnosis only, before any correction" in html
    assert "No decision layer ran" in html
    assert not layout.decision.exists()


def test_exports_round_trip_with_preserved_ids_and_declared_scale(
    completed_run, tmp_path: Path
) -> None:
    import anndata as ad

    config_path, project, layout = completed_run
    from batchdx.cyto.config import load_config
    from batchdx.cyto.io import load_dataset

    config = load_config(config_path)
    dataset = load_dataset(config, project)

    stored = ad.read_h5ad(layout.corrected_h5ad("uncorrected"))
    assert list(map(str, stored.obs["cell_id"])) == list(map(str, dataset.obs["cell_id"]))
    assert list(map(str, stored.var_names)) == list(map(str, dataset.adata.var_names))
    matrix = np.asarray(stored.layers[str(stored.uns["expression_layer"])])
    assert np.allclose(matrix, dataset.expression())
    assert stored.uns["output_scale"] == dataset.processing.scale

    written = export_fcs_per_sample(dataset, matrix, tmp_path / "fcs",
                                    candidate_id="uncorrected")
    sample_id = sorted(written)[0]
    contents = read_fcs(written[sample_id])
    mask = (dataset.obs["sample_id"] == sample_id).to_numpy()
    assert contents.n_events == int(mask.sum())
    expected = np.sinh(dataset.expression()[mask, :]) * 5.0
    assert np.allclose(contents.events, expected, rtol=1e-5, atol=1e-5)

    keywords = {k.upper(): v for k, v in contents.keywords.items()}
    processing = next(v for k, v in keywords.items() if k.endswith("BATCHDX_PROCESSING"))
    assert "do NOT compensate or transform these values again" in processing
    first_cell = next(v for k, v in keywords.items() if k.endswith("BATCHDX_FIRST_CELL_ID"))
    assert first_cell == str(dataset.obs.loc[mask, "cell_id"].iloc[0])


def test_status_file_reports_every_planned_unit_not_just_the_ones_that_ran(
    completed_run,
) -> None:
    _, _, layout = completed_run
    payload = json.loads(layout.status.read_text())
    units = {(e["stage"], e["unit"]): e for e in payload["entries"]}

    assert ("correct", "cytonorm_reference") in units
    assert units[("correct", "cytonorm_reference")]["status"] in ("blocked", "failed")
    assert units[("correct", "cytonorm_reference")]["reason"]

    # detectors that were not requested are reported as skipped, with the reason
    assert units[("diagnose", "cytobatchflagr")]["status"] == "skipped"
    assert "not requested" in units[("diagnose", "cytobatchflagr")]["reason"]
    assert payload["analytical_fingerprint"]


def test_report_settings_do_not_change_the_analytical_fingerprint(
    example_project: Path,
) -> None:
    from batchdx.cyto.config import load_config

    config_path = example_project / "config" / "config.yaml"
    before = load_config(config_path).analytical_fingerprint()
    config_path.write_text(
        config_path.read_text().replace("title: batchdx-cyto", "title: a different title")
    )
    after = load_config(config_path).analytical_fingerprint()
    assert before == after
