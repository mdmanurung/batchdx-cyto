"""Stage 1 acceptance gate: installed package, example elsewhere, Snakemake end to end."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.slow


def _run(args: list[str], cwd: Path, env: dict | None = None) -> subprocess.CompletedProcess:
    full = dict(os.environ)
    full.update(env or {})
    return subprocess.run(args, cwd=cwd, env=full, capture_output=True, text=True, check=False)


def _cli(*args: str) -> list[str]:
    return [sys.executable, "-m", "batchdx.cyto.cli", *args]


def test_cli_module_entry_point_is_importable_and_reports_what_it_ships(tmp_path: Path) -> None:
    proc = _run(_cli("info", "--json"), cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert set(payload["backends"]) >= {"none", "cycombine", "cytonorm"}
    assert "acquisition_drift" in payload["detectors"]
    assert Path(payload["workflow"]).name == "Snakefile"
    # capability claims are explicit, not implied
    assert payload["backends"]["cycombine"]["fit_apply_verified"] is False
    assert "method" in payload["backends"]["cycombine"]["rejected_parameters"]


def test_example_initialises_every_file_it_references(tmp_path: Path) -> None:
    project = tmp_path / "unrelated" / "place"
    proc = _run(_cli("init", str(project), "--cells-per-sample", "300"), cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr + proc.stdout

    import yaml

    config = yaml.safe_load((project / "config" / "config.yaml").read_text())
    for key in ("path", "samples", "markers", "processing_manifest"):
        assert (project / config["input"][key]).exists(), key
    assert (project / config["evaluation"]["gating_thresholds"]).exists()
    assert (project / "README.md").exists()

    validated = _run(_cli("validate", "--config", "config/config.yaml"), cwd=project)
    assert validated.returncode == 0, validated.stderr + validated.stdout
    assert "input OK" in validated.stdout


def test_snakemake_runs_the_slice_and_a_second_invocation_does_no_analytical_work(
    tmp_path: Path,
) -> None:
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not on PATH")

    project = tmp_path / "project"
    assert _run(_cli("init", str(project), "--cells-per-sample", "300"),
                cwd=tmp_path).returncode == 0

    from batchdx.cyto.paths import snakefile

    outdir = project / "results" / "run1"
    command = [
        "snakemake", "--snakefile", str(snakefile()), "--directory", str(project),
        "--cores", "1", "--config", "config=config/config.yaml",
        f"outdir={outdir}", f"base_dir={project}", "mode=benchmark",
    ]
    first = _run(command, cwd=project)
    assert first.returncode == 0, first.stdout + first.stderr

    assert (outdir / "manifests" / "dataset_manifest.json").exists()
    assert (outdir / "design" / "design_audit.json").exists()
    assert (outdir / "decisions" / "decision_result.json").exists()
    assert (outdir / "report" / "report.html").exists()
    assert (outdir / "status" / "run_status.json").exists()

    # the identity candidate is a genuine candidate, not a placeholder
    identity = json.loads(
        (outdir / "candidates" / "uncorrected" / "correction_result.json").read_text()
    )
    assert identity["status"] == "success"
    assert identity["validation"]["identity_numerically_unchanged"] is True
    assert (outdir / "candidates" / "uncorrected" / "corrected.h5ad").exists()

    second = _run(command, cwd=project)
    assert second.returncode == 0, second.stdout + second.stderr
    assert "Nothing to be done" in second.stderr + second.stdout


def test_a_bad_configuration_propagates_a_useful_error_and_exit_code(tmp_path: Path) -> None:
    project = tmp_path / "project"
    assert _run(_cli("init", str(project), "--cells-per-sample", "200"),
                cwd=tmp_path).returncode == 0

    config_path = project / "config" / "config.yaml"
    config_path.write_text(config_path.read_text() + "\nnot_a_real_key: 1\n")

    proc = _run(_cli("validate", "--config", "config/config.yaml"), cwd=project)
    assert proc.returncode != 0
    assert "Extra inputs are not permitted" in proc.stderr
    assert "not_a_real_key" in proc.stderr


def test_correct_requires_an_explicitly_selected_candidate(tmp_path: Path) -> None:
    project = tmp_path / "project"
    assert _run(_cli("init", str(project), "--cells-per-sample", "200"),
                cwd=tmp_path).returncode == 0

    without = _run(_cli("correct", "--config", "config/config.yaml",
                        "--outdir", "results/x"), cwd=project)
    assert without.returncode != 0
    assert "candidate" in (without.stderr + without.stdout).lower()

    unknown = _run(_cli("correct", "--config", "config/config.yaml",
                        "--outdir", "results/x", "--candidate", "does_not_exist"),
                   cwd=project)
    assert unknown.returncode != 0
    assert "is not in the configuration" in unknown.stderr


def test_packaged_workflow_is_found_from_outside_the_source_checkout(tmp_path: Path) -> None:
    from batchdx.cyto.paths import config_root, examples_root, snakefile, templates_root

    for path in (snakefile(), config_root(), examples_root(), templates_root()):
        assert path.exists(), path
    assert (templates_root() / "report.html.j2").exists()
