"""Stage 11 (partial): Snakemake execution behaviour.

Stage 11 as specified — SLURM profiles, Apptainer images with pinned digests,
checkpoint-driven branches, environment prefetching — is **not** complete in this build
(see ``docs/implementation_status.md``).  What is tested here is the execution
behaviour the Tier A slice depends on: static expansion from the validated
configuration, real artifact dependencies, no-op reruns, invalidation on a changed
configuration, recovery from a deleted output, and a partial report when a candidate
cannot run.

A dry run is not execution validation, so every test below really executes the
workflow.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def project(tmp_path_factory) -> Path:
    from batchdx.cyto.examples import materialise_example

    path = tmp_path_factory.mktemp("workflow")
    materialise_example("paired_reference", path, cells_per_sample=250)
    return path


def _snakemake(project: Path, outdir: Path, mode: str = "benchmark",
               config_entries: list[str] | None = None,
               extra: list[str] | None = None) -> subprocess.CompletedProcess:
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not on PATH")
    from batchdx.cyto.paths import snakefile

    # Snakemake takes one --config group; a second occurrence replaces the first, so
    # every key/value pair goes into the same group.
    command = [
        "snakemake", "--snakefile", str(snakefile()), "--directory", str(project),
        "--cores", "1", "--config", "config=config/config.yaml",
        f"outdir={outdir}", f"base_dir={project}", f"mode={mode}",
    ] + (config_entries or []) + (extra or [])
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
    return subprocess.run(command, cwd=project, capture_output=True, text=True,
                          check=False, env=env)


def test_clean_execution_then_a_no_op_rerun(project: Path) -> None:
    outdir = project / "results" / "clean"
    first = _snakemake(project, outdir)
    assert first.returncode == 0, first.stdout + first.stderr
    assert outdir.joinpath("report", "report.html").exists()

    second = _snakemake(project, outdir)
    assert second.returncode == 0, second.stdout + second.stderr
    assert "Nothing to be done" in second.stdout + second.stderr


def test_deleting_a_real_artifact_forces_recomputation(project: Path) -> None:
    outdir = project / "results" / "recover"
    assert _snakemake(project, outdir).returncode == 0

    target = outdir / "metrics" / "uncorrected" / "metric_result.json"
    decision = outdir / "decisions" / "decision_result.json"
    assert target.exists() and decision.exists()

    # Remove a downstream result and the intermediate it genuinely depends on. Because
    # `decide` declares the metric files as real inputs (rather than depending on a
    # success sentinel), rebuilding the decision has to rebuild the metric result.
    decision.unlink()
    target.unlink()

    rerun = _snakemake(project, outdir)
    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    assert target.exists(), (
        "the metric result is a declared input of the decide rule, so it must be "
        "recomputed; a success sentinel alone would not force this"
    )
    assert decision.exists()


def test_changing_the_configuration_invalidates_downstream_results(project: Path) -> None:
    outdir = project / "results" / "invalidate"
    assert _snakemake(project, outdir).returncode == 0
    fingerprint_before = json.loads(outdir.joinpath("config.resolved.json").read_text())[
        "analytical_fingerprint"
    ]

    config_path = project / "config" / "config.yaml"
    original = config_path.read_text()
    config_path.write_text(original.replace("max_cells_per_sample: 1200",
                                            "max_cells_per_sample: 150"))
    try:
        # the configuration is an input of the validate rule, so Snakemake reruns it
        rerun = _snakemake(project, outdir)
        assert rerun.returncode == 0, rerun.stdout + rerun.stderr
        fingerprint_after = json.loads(outdir.joinpath("config.resolved.json").read_text())[
            "analytical_fingerprint"
        ]
        assert fingerprint_after != fingerprint_before
    finally:
        config_path.write_text(original)


def test_diagnose_mode_stops_before_any_correction(project: Path) -> None:
    outdir = project / "results" / "diagnose"
    proc = _snakemake(project, outdir, mode="diagnose")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert outdir.joinpath("report", "report.html").exists()
    assert not outdir.joinpath("candidates").exists()
    assert not outdir.joinpath("decisions", "decision_result.json").exists()
    assert "diagnosis only" in outdir.joinpath("report", "report.html").read_text()


def test_correct_mode_runs_only_the_selected_candidate(project: Path) -> None:
    outdir = project / "results" / "correct_one"
    proc = _snakemake(project, outdir, mode="correct",
                      config_entries=["candidate=uncorrected"])
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert outdir.joinpath("candidates", "uncorrected", "correction_result.json").exists()
    assert not outdir.joinpath("candidates", "cycombine_default").exists()


def test_correct_mode_without_a_candidate_fails_at_parse_time(project: Path) -> None:
    proc = _snakemake(project, project / "results" / "nope", mode="correct")
    assert proc.returncode != 0
    assert "explicitly selected candidate" in proc.stdout + proc.stderr


def test_an_unavailable_method_produces_a_partial_report_not_a_failed_run(
    project: Path,
) -> None:
    outdir = project / "results" / "partial"
    proc = _snakemake(project, outdir)
    assert proc.returncode == 0, "a scientific skip is not a workflow failure"

    status = json.loads(outdir.joinpath("status", "run_status.json").read_text())
    blocked = [e for e in status["entries"] if e["status"] in ("blocked", "ineligible")]
    assert blocked, "the roster must record the candidates that could not run"
    assert all(e["reason"] for e in blocked)

    html = outdir.joinpath("report", "report.html").read_text()
    for entry in blocked:
        assert entry["unit"] in html


def test_no_rule_hides_a_failure_behind_true_or_blanket_suppression(repo_root: Path) -> None:
    snakefile = (repo_root / "workflow" / "Snakefile").read_text()
    assert "|| true" not in snakefile
    assert "except Exception" not in snakefile
    assert "try:" not in snakefile


def test_rules_declare_logs_threads_resources_and_benchmarks(repo_root: Path) -> None:
    text = (repo_root / "workflow" / "Snakefile").read_text()
    for rule in ("validate", "design_audit", "partition", "diagnose", "correct",
                 "metrics", "decide", "report"):
        assert f"rule {rule}:" in text
    assert text.count("log:") >= 8
    assert text.count("threads:") >= 8
    assert text.count("resources:") >= 8
    assert text.count("benchmark:") >= 7


def test_the_candidate_roster_comes_from_the_configuration_not_the_results_tree(
    project: Path,
) -> None:
    outdir = project / "results" / "roster"
    assert _snakemake(project, outdir).returncode == 0

    stray = outdir / "candidates" / "not_in_the_config"
    stray.mkdir(parents=True, exist_ok=True)
    (stray / "correction_result.json").write_text("{}")

    rerun = _snakemake(project, outdir)
    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    status = json.loads(outdir.joinpath("status", "run_status.json").read_text())
    units = {e["unit"] for e in status["entries"] if e["stage"] == "correct"}
    assert "not_in_the_config" not in units


def test_slurm_profile_is_present_but_labelled_untested(repo_root: Path) -> None:
    profile = repo_root / "workflow" / "profiles" / "slurm" / "config.yaml"
    assert profile.exists()
    text = profile.read_text()
    assert "UNTESTED" in text, (
        "the SLURM profile has not been executed on a cluster in this build and must "
        "say so rather than implying it has"
    )
