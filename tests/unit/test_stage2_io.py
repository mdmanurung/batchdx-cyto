"""Stage 2 acceptance gate: ingestion, immutability and manifests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto.config import load_config
from batchdx.cyto.io import DatasetValidationError, ManifestError, load_dataset
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID
from batchdx.cyto.io.export import export_fcs_per_sample, plan_inversion
from batchdx.cyto.io.fcs import PROCESSING_KEYWORD, SCALE_KEYWORD, FCSError, read_fcs, write_fcs
from batchdx.cyto.io.manifests import read_marker_table, read_sample_table


def _config(project: Path):
    return load_config(project / "config" / "config.yaml")


def test_h5ad_roundtrip_preserves_marker_and_channel_order(example_project: Path) -> None:
    config = _config(example_project)
    dataset = load_dataset(config, example_project)
    declared = [m.feature_id for m in read_marker_table(example_project / "config/markers.tsv")]
    assert list(map(str, dataset.adata.var_names)) == declared
    assert [dataset.markers[f].channel for f in declared] == [f"Ch_{f}" for f in declared]


def test_cell_identity_is_stable_and_unique(example_project: Path) -> None:
    config = _config(example_project)
    dataset = load_dataset(config, example_project)
    ids = dataset.obs[CELL_ID]
    assert ids.is_unique
    again = load_dataset(config, example_project)
    assert list(again.obs[CELL_ID]) == list(ids)


def test_baseline_is_immutable_under_read_operations(example_project: Path) -> None:
    config = _config(example_project)
    dataset = load_dataset(config, example_project)
    before = dataset.baseline_hash()

    block = dataset.expression(features=dataset.marker_ids()[:2])
    block += 100.0  # a caller mutating what it was handed must not reach the baseline

    dataset.assert_baseline_unchanged()
    assert dataset.baseline_hash() == before


def test_duplicate_sample_ids_are_rejected(example_project: Path) -> None:
    path = example_project / "config" / "samples.tsv"
    lines = path.read_text().splitlines()
    path.write_text("\n".join(lines + [lines[1]]) + "\n")
    with pytest.raises(ManifestError, match="duplicate sample_id"):
        read_sample_table(path)


def test_ambiguous_channel_mapping_is_rejected(example_project: Path) -> None:
    path = example_project / "config" / "markers.tsv"
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    channel_col = header.index("channel")
    row = lines[2].split("\t")
    row[channel_col] = lines[1].split("\t")[channel_col]
    lines[2] = "\t".join(row)
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ManifestError, match="map to more than one feature_id"):
        read_marker_table(path)


def test_undeclared_scale_is_rejected(example_project: Path) -> None:
    path = example_project / "config" / "processing.json"
    payload = json.loads(path.read_text())
    payload.pop("scale")
    path.write_text(json.dumps(payload))
    with pytest.raises(ManifestError, match="must be declared explicitly"):
        load_dataset(_config(example_project), example_project)


def test_double_transformation_is_rejected(example_project: Path) -> None:
    config_path = example_project / "config" / "config.yaml"
    text = config_path.read_text().replace(
        "preprocessing:\n  enabled: false",
        "preprocessing:\n  enabled: true\n  transformation: asinh\n"
        "  cofactor_by_panel:\n    panelA: 5.0",
    )
    config_path.write_text(text)
    with pytest.raises(DatasetValidationError, match="refusing to transform twice"):
        load_dataset(load_config(config_path), example_project)


def test_already_corrected_input_is_refused(example_project: Path) -> None:
    path = example_project / "config" / "processing.json"
    payload = json.loads(path.read_text())
    payload["corrected"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(DatasetValidationError, match="corrected=true"):
        load_dataset(_config(example_project), example_project)


def test_missing_marker_is_not_treated_as_zero(example_project: Path) -> None:
    path = example_project / "config" / "markers.tsv"
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    row[header.index("feature_id")] = "CD_NOT_MEASURED"
    row[header.index("channel")] = "Ch_CD_NOT_MEASURED"
    lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(DatasetValidationError, match="Declared but absent"):
        load_dataset(_config(example_project), example_project)


def test_nonfinite_values_are_rejected(example_project: Path, tmp_path: Path) -> None:
    import anndata as ad

    path = example_project / "data" / "preprocessed.h5ad"
    adata = ad.read_h5ad(path)
    layer = np.asarray(adata.layers["preprocessed"], dtype=np.float64)
    layer[3, 1] = np.nan
    adata.layers["preprocessed"] = layer
    adata.write_h5ad(path)
    with pytest.raises(DatasetValidationError, match="non-finite"):
        load_dataset(_config(example_project), example_project)


def test_missing_expression_layer_names_what_is_available(example_project: Path) -> None:
    config_path = example_project / "config" / "config.yaml"
    config_path.write_text(
        config_path.read_text().replace("expression_layer: preprocessed",
                                        "expression_layer: not_a_layer")
    )
    with pytest.raises(DatasetValidationError, match="available layers"):
        load_dataset(load_config(config_path), example_project)


def test_multiple_panels_must_be_run_separately(example_project: Path) -> None:
    path = example_project / "config" / "markers.tsv"
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    panel_col = header.index("panel_id")
    row = lines[2].split("\t")
    row[panel_col] = "panelB"
    lines[2] = "\t".join(row)
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(DatasetValidationError, match="more than one panel_id"):
        load_dataset(_config(example_project), example_project)


def test_qc_history_and_transformation_provenance_are_carried(example_project: Path) -> None:
    dataset = load_dataset(_config(example_project), example_project)
    assert dataset.manifest.processing.transformation == "asinh"
    assert dataset.manifest.processing.transformation_parameters["cofactor"] == 5.0
    assert dataset.manifest.processing.qc_history
    assert dataset.manifest.provenance[0].content_hash


def test_fcs_write_read_roundtrip_and_keyword_policy(example_project: Path, tmp_path: Path) -> None:
    dataset = load_dataset(_config(example_project), example_project)
    written = export_fcs_per_sample(
        dataset, dataset.expression(), tmp_path / "fcs", candidate_id="uncorrected"
    )
    assert written

    sample_id = sorted(written)[0]
    contents = read_fcs(written[sample_id])
    assert contents.channel_names == [
        dataset.markers[f].channel for f in map(str, dataset.adata.var_names)
    ]
    assert PROCESSING_KEYWORD in {k.upper() for k in contents.keywords} or any(
        k.upper().endswith(PROCESSING_KEYWORD) for k in contents.keywords
    )

    plan = plan_inversion(dataset.processing)
    assert plan.invertible and plan.cofactor == 5.0

    mask = (dataset.obs[SAMPLE_ID] == sample_id).to_numpy()
    expected = np.sinh(dataset.expression()[mask, :]) * 5.0
    assert np.allclose(contents.events, expected, rtol=1e-5, atol=1e-5)

    scale_values = [v for k, v in contents.keywords.items() if k.upper().endswith(SCALE_KEYWORD)]
    assert scale_values and "inverted" in scale_values[0]


def test_fcs_export_refuses_to_overwrite_a_source_file(tmp_path: Path) -> None:
    path = tmp_path / "already.fcs"
    write_fcs(path, np.zeros((4, 2)), ["A", "B"], scale="test", processing_note="test")
    with pytest.raises(FCSError, match="refusing to overwrite"):
        write_fcs(path, np.zeros((4, 2)), ["A", "B"], scale="test", processing_note="test")


def test_uninvertible_transformation_keeps_transformed_values_and_says_so(
    example_project: Path,
) -> None:
    dataset = load_dataset(_config(example_project), example_project)
    dataset.processing.transformation = "logicle"
    plan = plan_inversion(dataset.processing)
    assert not plan.invertible
    assert "NOT instrument measurements" in plan.note
    values = np.array([[1.0, 2.0]])
    assert np.array_equal(plan.apply(values), values)
