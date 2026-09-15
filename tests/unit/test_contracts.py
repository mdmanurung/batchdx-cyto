"""Cross-stage contracts: versioned, strict, and round-trip tested across the R boundary."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto._version import CONTRACT_SCHEMA_VERSION
from batchdx.cyto.io.binary import read_matrix, write_matrix
from batchdx.cyto.schemas import all_contracts, export_json_schemas
from batchdx.cyto.schemas.contracts import (
    CandidateSpecRecord,
    CorrectionResult,
    DecisionResult,
    DesignAudit,
    DetectorResult,
    MetricRecord,
    MetricResult,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import CorrectionStatus, LabelTier, MetricStatus

REQUIRED_METRIC_FIELDS = {
    "run_id", "candidate_id", "metric_namespace", "metric_name", "metric_version",
    "stage", "marker", "population", "comparison_id", "value", "direction", "scale",
    "status", "reason", "n_samples", "n_donors", "n_cells", "coverage",
    "aggregation_weight", "uncertainty_low", "uncertainty_high", "uncertainty_kind",
    "label_tier", "label_confidence", "representation", "provenance",
}


def test_all_seven_contracts_exist_and_are_versioned() -> None:
    contracts = all_contracts()
    assert set(contracts) == {
        "DatasetManifest", "DesignAudit", "SamplingManifest", "DetectorResult",
        "CorrectionResult", "MetricResult", "DecisionResult",
    }
    for name, model in contracts.items():
        assert "schema_version" in model.model_fields, name


def test_contracts_reject_unknown_keys(frozen) -> None:
    _, _, design, _ = frozen
    payload = json.loads(design.to_json())
    payload["a_field_that_does_not_exist"] = 1
    with pytest.raises(Exception, match="(?i)extra"):
        DesignAudit.model_validate(payload)


def test_contracts_round_trip_through_json(frozen, tmp_path: Path) -> None:
    _, dataset, design, sampling = frozen
    for record, model in (
        (dataset.manifest, type(dataset.manifest)),
        (design, DesignAudit),
        (sampling, SamplingManifest),
    ):
        path = tmp_path / f"{model.__name__}.json"
        record.write(path)
        again = model.read(path)
        assert again.model_dump(mode="json") == record.model_dump(mode="json")


def test_metric_record_carries_every_field_the_contract_requires() -> None:
    assert REQUIRED_METRIC_FIELDS <= set(MetricRecord.model_fields)


def test_correction_record_distinguishes_success_ineligible_and_failed() -> None:
    values = {s.value for s in CorrectionStatus}
    assert values == {"success", "ineligible", "failed", "not_available"}

    spec = CandidateSpecRecord(candidate_id="c", backend="none")
    failed = CorrectionResult(candidate=spec, status=CorrectionStatus.FAILED,
                              reason="it failed")
    assert failed.output_path is None
    assert failed.markers_changed == []
    assert "correction_scope" in CandidateSpecRecord.model_fields
    assert "input_hashes" in CandidateSpecRecord.model_fields


def test_json_schemas_are_exported_for_the_r_side(tmp_path: Path) -> None:
    written = export_json_schemas(tmp_path)
    assert set(written) == set(all_contracts())
    for path in written.values():
        schema = json.loads(path.read_text())
        assert schema["type"] == "object"
        assert schema.get("additionalProperties") is False


@pytest.mark.requires_r
def test_r_reads_every_contract_and_the_binary_matrix_format(tmp_path: Path) -> None:
    """The Python/R boundary is round-trip tested, not merely conventional."""

    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not on PATH")

    schema_dir = tmp_path / "schemas"
    export_json_schemas(schema_dir)

    detector = DetectorResult(
        detector="composition", namespace="common", run_id="r", status="ok",
    )
    detector.write(tmp_path / "DetectorResult.json")

    decision = DecisionResult(run_id="r", policy="conservative")
    decision.write(tmp_path / "DecisionResult.json")

    metrics = MetricResult(
        run_id="r", candidate_id="c", metric_group="common",
        records=[MetricRecord(
            run_id="r", candidate_id="c", metric_namespace="common",
            metric_name="emd_within_population", metric_version="1.0.0",
            stage="distribution_alignment", value=0.25, status=MetricStatus.OK,
            label_tier=LabelTier.GATING, label_confidence="medium",
        )],
    )
    metrics.write(tmp_path / "MetricResult.json")

    matrix = np.arange(12, dtype=float).reshape(4, 3) + 0.5
    write_matrix(tmp_path / "m.f64", matrix, ["c1", "c2", "c3", "c4"], ["A", "B", "C"])

    script = tmp_path / "check.R"
    script.write_text(
        'suppressPackageStartupMessages(library(jsonlite))\n'
        'args <- commandArgs(trailingOnly = TRUE)\n'
        'dir <- args[[1]]; lib <- args[[2]]\n'
        'source(file.path(lib, "binmat.R"))\n'
        'stopifnot(length(list.files(file.path(dir, "schemas"))) == 7)\n'
        'det <- fromJSON(file.path(dir, "DetectorResult.json"))\n'
        'stopifnot(det$detector == "composition", det$status == "ok")\n'
        'dec <- fromJSON(file.path(dir, "DecisionResult.json"))\n'
        'stopifnot(dec$policy == "conservative")\n'
        'met <- fromJSON(file.path(dir, "MetricResult.json"))\n'
        'stopifnot(met$records$metric_name[[1]] == "emd_within_population")\n'
        'stopifnot(abs(met$records$value[[1]] - 0.25) < 1e-12)\n'
        'stopifnot(met$records$label_tier[[1]] == "gating")\n'
        'bm <- read_binmat(file.path(dir, "m.f64"))\n'
        'stopifnot(identical(dim(bm$matrix), c(4L, 3L)))\n'
        'stopifnot(identical(colnames(bm$matrix), c("A", "B", "C")))\n'
        'stopifnot(abs(bm$matrix[2, 3] - 5.5) < 1e-12)\n'
        '# write it back and let Python verify the return trip\n'
        'write_binmat(file.path(dir, "back.f64"), bm$matrix * 2, rownames(bm$matrix),\n'
        '             colnames(bm$matrix))\n'
        'cat("ok\\n")\n',
        encoding="utf-8",
    )

    from batchdx.cyto.paths import workflow_root

    proc = subprocess.run(
        ["Rscript", "--vanilla", str(script), str(tmp_path),
         str(workflow_root() / "scripts" / "lib")],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ok" in proc.stdout

    returned, row_ids, col_names, _ = read_matrix(tmp_path / "back.f64")
    assert np.allclose(returned, matrix * 2)
    assert row_ids == ["c1", "c2", "c3", "c4"]
    assert col_names == ["A", "B", "C"]


def test_binary_matrix_rejects_a_shape_mismatch(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="does not match"):
        write_matrix(tmp_path / "bad.f64", np.zeros((2, 2)), ["a"], ["x", "y"])


def test_downstream_never_infers_success_from_a_directory(
    example_project: Path, tmp_path: Path
) -> None:
    """A candidate directory with no result record must not read as a success."""

    from batchdx.cyto.pipeline import (
        RunLayout,
        stage_design,
        stage_partition,
        stage_status,
        stage_validate,
    )

    layout = RunLayout(tmp_path / "run")
    config_path = example_project / "config" / "config.yaml"
    stage_validate(config_path, example_project, layout)
    stage_design(config_path, example_project, layout)
    stage_partition(config_path, example_project, layout)

    layout.candidate_dir("cycombine_default").mkdir(parents=True, exist_ok=True)
    stage_status(config_path, layout)

    payload = json.loads(layout.status.read_text())
    entry = next(e for e in payload["entries"]
                 if e["stage"] == "correct" and e["unit"] == "cycombine_default")
    assert entry["status"] == "failed"
    assert "no correction result file" in entry["reason"]


def test_contract_schema_version_is_stamped_on_every_record(frozen) -> None:
    _, dataset, design, sampling = frozen
    for record in (dataset.manifest, design, sampling):
        assert record.schema_version == CONTRACT_SCHEMA_VERSION
