"""Candidate failures are reportable without scoring stale or invalid outputs."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml

from batchdx.cyto import pipeline
from batchdx.cyto.adapters import BackendError, BackendNotAvailable, get_backend
from batchdx.cyto.adapters.validate import OutputValidationError
from batchdx.cyto.config import ConfigError
from batchdx.cyto.io.dataset import DatasetValidationError
from batchdx.cyto.schemas.contracts import CorrectionResult, MetricResult
from batchdx.cyto.schemas.enums import CorrectionStatus


@pytest.fixture
def candidate_run(example_project: Path, tmp_path: Path, monkeypatch):
    config_path = example_project / "config" / "config.yaml"
    raw = yaml.safe_load(config_path.read_text())
    raw["correction"]["candidates"] = [
        {"id": "uncorrected", "backend": "none"},
        {"id": "boom", "backend": "cycombine"},
        {"id": "working", "backend": "fixture_offset_removal"},
    ]
    raw["upstream"]["rscript"] = "Rscript-does-not-exist"
    config_path.write_text(yaml.safe_dump(raw))
    layout = pipeline.RunLayout(tmp_path / "run")
    pipeline.stage_validate(config_path, example_project, layout)
    pipeline.stage_design(config_path, example_project, layout)
    pipeline.stage_partition(config_path, example_project, layout)
    monkeypatch.setattr(get_backend("cycombine"), "availability", lambda context: [])
    return config_path, example_project, layout


@pytest.mark.parametrize(("error", "status"), [
    (OutputValidationError("invalid corrected output"), CorrectionStatus.FAILED),
    (FileNotFoundError("expected output missing"), CorrectionStatus.FAILED),
    (json.JSONDecodeError("invalid result JSON", "{", 1), CorrectionStatus.FAILED),
    (subprocess.TimeoutExpired("correction", 1), CorrectionStatus.FAILED),
    (BackendError("backend execution failed"), CorrectionStatus.FAILED),
    (BackendNotAvailable("backend is not installed"), CorrectionStatus.NOT_AVAILABLE),
])
def test_candidate_errors_become_records_and_keep_tracebacks(candidate_run, monkeypatch, error, status):
    config_path, project, layout = candidate_run

    def fail(*args):
        raise error

    monkeypatch.setattr(get_backend("cycombine"), "run", fail)
    pipeline.stage_correct("boom", config_path, project, layout)
    result = CorrectionResult.read(layout.correction_result("boom"))
    assert result.status is status
    assert str(error) in result.reason
    assert result.output_path is None
    assert result.validation["produced_output"] is False
    # The Snakemake rule declares corrected.h5ad as an output, so a failed candidate
    # writes the unmodified baseline in its place; the failure lives in the result,
    # which is what every consumer reads before the file.
    assert layout.corrected_h5ad("boom").exists()
    log = (layout.root / "logs" / "correct_boom.log").read_text()
    assert type(error).__name__ in log
    assert str(error) in log


def test_invalid_candidate_leaves_other_candidates_and_report_usable(candidate_run, monkeypatch):
    config_path, project, layout = candidate_run

    def invalid_output(dataset, candidate_id, parameters, context):
        from batchdx.cyto.adapters.validate import validate_output

        matrix = dataset.expression()
        matrix[0, 0] = np.inf
        validate_output(dataset, matrix, dataset.marker_ids("correction"), candidate_id=candidate_id)

    monkeypatch.setattr(get_backend("cycombine"), "run", invalid_output)
    for candidate in ("boom", "uncorrected", "working"):
        pipeline.stage_correct(candidate, config_path, project, layout)
        pipeline.stage_metrics(candidate, config_path, project, layout)
    for candidate in ("uncorrected", "working"):
        assert CorrectionResult.read(layout.correction_result(candidate)).status is CorrectionStatus.SUCCESS
        assert MetricResult.read(layout.metric_result(candidate)).records
    metrics = MetricResult.read(layout.metric_result("boom"))
    assert metrics.status == "ineligible"
    assert metrics.records == []
    pipeline.stage_decide(config_path, project, layout)
    pipeline.stage_status(config_path, layout)
    pipeline.stage_report(config_path, project, layout)
    html = layout.report.read_text()
    assert "boom" in html and "non-finite" in html
    entries = json.loads(layout.status.read_text())["entries"]
    assert any(e["unit"] == "boom" and e["status"] == "failed" for e in entries)


def test_failed_retry_invalidates_previous_output_and_metrics(candidate_run, monkeypatch):
    config_path, project, layout = candidate_run
    backend = get_backend("cycombine")
    monkeypatch.setattr(backend, "run", get_backend("none").run)
    pipeline.stage_correct("boom", config_path, project, layout)
    pipeline.stage_metrics("boom", config_path, project, layout)
    assert layout.corrected_h5ad("boom").exists()
    assert MetricResult.read(layout.metric_result("boom")).records

    def fail(*args):
        raise OutputValidationError("retry output invalid")

    monkeypatch.setattr(backend, "run", fail)
    pipeline.stage_correct("boom", config_path, project, layout)
    assert CorrectionResult.read(layout.correction_result("boom")).status is CorrectionStatus.FAILED
    assert layout.corrected_h5ad("boom").exists()
    assert not layout.metric_result("boom").exists()
    pipeline.stage_metrics("boom", config_path, project, layout)
    assert MetricResult.read(layout.metric_result("boom")).records == []


def test_partial_export_is_removed_and_failure_is_recorded(candidate_run, monkeypatch):
    config_path, project, layout = candidate_run
    monkeypatch.setattr(get_backend("cycombine"), "run", get_backend("none").run)
    real_write = pipeline.write_candidate_output

    def fail_write(dataset, matrix, path, **kwargs):
        # Fail only the correction export (the call that carries output_scale); the
        # failure candidate's baseline write must succeed so the declared output exists.
        if "output_scale" in kwargs:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"incomplete H5AD")
            raise OSError("candidate export failed")
        return real_write(dataset, matrix, path, **kwargs)

    monkeypatch.setattr(pipeline, "write_candidate_output", fail_write)
    pipeline.stage_correct("boom", config_path, project, layout)
    assert CorrectionResult.read(layout.correction_result("boom")).status is CorrectionStatus.FAILED
    assert layout.corrected_h5ad("boom").read_bytes() != b"incomplete H5AD"


@pytest.mark.parametrize("invalidity", ["shape", "non-finite"])
def test_unvalidated_backend_output_is_contained_at_export(candidate_run, monkeypatch, invalidity):
    config_path, project, layout = candidate_run

    def bad_output(dataset, candidate_id, parameters, context):
        matrix, result = get_backend("none").run(dataset, candidate_id, {}, context)
        if invalidity == "shape":
            matrix = matrix[:1]
        else:
            matrix[0, 0] = np.nan
        return matrix, result

    monkeypatch.setattr(get_backend("cycombine"), "run", bad_output)
    pipeline.stage_correct("boom", config_path, project, layout)
    result = CorrectionResult.read(layout.correction_result("boom"))
    assert result.status is CorrectionStatus.FAILED
    assert invalidity in result.reason
    assert layout.corrected_h5ad("boom").exists()


def test_non_success_return_never_exports_a_matrix(candidate_run, monkeypatch):
    from batchdx.cyto.adapters import failure_result, make_spec

    config_path, project, layout = candidate_run

    def failed_result(dataset, candidate_id, parameters, context):
        spec = make_spec(candidate_id, "cycombine", None, parameters, context.seed, dataset, [])
        return None, failure_result(spec, CorrectionStatus.FAILED, "no usable output")

    monkeypatch.setattr(get_backend("cycombine"), "run", failed_result)
    pipeline.stage_correct("boom", config_path, project, layout)
    assert CorrectionResult.read(layout.correction_result("boom")).status is CorrectionStatus.FAILED
    assert layout.corrected_h5ad("boom").exists()


def test_availability_failure_cannot_mask_original_candidate_error(candidate_run, monkeypatch):
    config_path, project, layout = candidate_run

    def fail(*args):
        raise BackendError("original backend failure")

    def unavailable(context):
        raise subprocess.TimeoutExpired("software probe", 1)

    backend = get_backend("cycombine")
    monkeypatch.setattr(backend, "run", fail)
    monkeypatch.setattr(backend, "availability", unavailable)
    pipeline.stage_correct("boom", config_path, project, layout)
    result = CorrectionResult.read(layout.correction_result("boom"))
    assert result.status is CorrectionStatus.FAILED
    assert "original backend failure" in result.reason
    assert result.software == []
    assert any("software probe" in w for w in result.warnings)


@pytest.mark.parametrize("error", [KeyboardInterrupt(), SystemExit(2), TypeError("programming defect")])
def test_cancellation_and_unexpected_defects_still_abort(candidate_run, monkeypatch, error):
    config_path, project, layout = candidate_run

    def fail(*args):
        raise error

    monkeypatch.setattr(get_backend("cycombine"), "run", fail)
    with pytest.raises(type(error)):
        pipeline.stage_correct("boom", config_path, project, layout)
    assert not layout.correction_result("boom").exists()


def test_invalid_shared_input_still_aborts(candidate_run, monkeypatch):
    config_path, project, layout = candidate_run

    def fail(*args):
        raise DatasetValidationError("invalid shared input")

    monkeypatch.setattr(pipeline, "_load", fail)
    with pytest.raises(DatasetValidationError, match="invalid shared input"):
        pipeline.stage_correct("boom", config_path, project, layout)
    assert not layout.correction_result("boom").exists()


def test_unknown_backend_is_rejected_before_dataset_validation(candidate_run):
    config_path, project, layout = candidate_run
    raw = yaml.safe_load(config_path.read_text())
    raw["correction"]["candidates"][1]["backend"] = "misspelled_backend"
    config_path.write_text(yaml.safe_dump(raw))
    layout = pipeline.RunLayout(layout.root / "unknown")
    with pytest.raises(ConfigError, match="unknown correction backend"):
        pipeline.stage_validate(config_path, project, layout)
    assert not layout.dataset_manifest.exists()


@pytest.mark.parametrize("candidate_id", ["", ".", "..", "../other", "other\\candidate"])
def test_candidate_ids_cannot_escape_the_directory_invalidated_on_retry(candidate_id):
    from pydantic import ValidationError

    from batchdx.cyto.config import CandidateConfig

    with pytest.raises(ValidationError, match="directory name"):
        CandidateConfig(id=candidate_id, backend="none")
