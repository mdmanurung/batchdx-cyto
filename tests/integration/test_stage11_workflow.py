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
import re
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

    profile = Path(snakefile()).parent / "profiles" / "local"

    # Snakemake takes one --config group; a second occurrence replaces the first, so
    # every key/value pair goes into the same group.
    command = [
        "snakemake", "--snakefile", str(snakefile()), "--directory", str(project),
        # workflow/profiles/local/config.yaml says the integration suite runs this
        # profile. It did not until this flag was added, so the profile header was the
        # thing that was wrong before; passing it makes the claim checkable rather than
        # aspirational. `--cores 1` below still wins over the profile's cores: 4.
        "--profile", str(profile),
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


def test_missing_backend_output_keeps_the_workflow_and_report_usable(tmp_path: Path) -> None:
    """Exercise the real adapter/stage/DAG boundary without requiring R packages."""
    import yaml

    from batchdx.cyto.examples import materialise_example

    materialise_example("paired_reference", tmp_path, cells_per_sample=100)
    config_path = tmp_path / "config" / "config.yaml"
    raw = yaml.safe_load(config_path.read_text())
    # The interpreter fixture answers version probes but omits the correction matrix.
    # It does not perform R execution or establish backend/scientific validity.
    interpreter = tmp_path / "fake-rscript"
    interpreter.write_text(
        f"#!{sys.executable}\nimport pathlib, sys\n"
        "if '-e' in sys.argv:\n    print('1.0.0')\n"
        "else:\n    print('fixture: correction output intentionally missing')\n"
        "    pathlib.Path(sys.argv[-1]).write_text('{\"status\":\"ok\"}')\n"
    )
    interpreter.chmod(0o755)
    raw["upstream"]["rscript"] = str(interpreter)
    raw["correction"]["candidates"] = [
        {"id": "uncorrected", "backend": "none"},
        {"id": "broken_output", "backend": "cycombine"},
    ]
    config_path.write_text(yaml.safe_dump(raw))
    outdir = tmp_path / "results" / "missing-output"
    proc = _snakemake(tmp_path, outdir)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    candidate_dir = outdir / "candidates" / "broken_output"
    correction = json.loads((candidate_dir / "correction_result.json").read_text())
    assert correction["status"] == "failed"
    assert not (candidate_dir / "corrected.h5ad").exists()
    metrics = json.loads((outdir / "metrics" / "broken_output" / "metric_result.json").read_text())
    assert metrics["records"] == [] and metrics["status"] == "ineligible"
    assert outdir.joinpath("candidates", "uncorrected", "corrected.h5ad").exists()
    html = outdir.joinpath("report", "report.html").read_text()
    assert "broken_output" in html and "failed" in html
    assert "FileNotFoundError" in outdir.joinpath("logs", "correct_broken_output.log").read_text()


def test_no_rule_hides_a_failure_behind_true_or_blanket_suppression(repo_root: Path) -> None:
    snakefile = (repo_root / "workflow" / "Snakefile").read_text()
    assert "|| true" not in snakefile
    assert "except Exception" not in snakefile
    assert "try:" not in snakefile


def test_rules_declare_logs_threads_resources_and_benchmarks(repo_root: Path) -> None:
    text = (repo_root / "workflow" / "Snakefile").read_text()
    for rule in ("validate", "design_audit", "partition", "diagnose", "correct",
                 "metrics", "decide", "status", "report"):
        assert f"rule {rule}:" in text
    assert text.count("log:") >= 9
    assert text.count("threads:") >= 9
    assert text.count("resources:") >= 9
    # Nine rules besides `all`, each with its own benchmark directive, so the threshold
    # has no slack left: dropping one directive anywhere now fails this test. It was
    # `>=7` while `rule status` had none, which is exactly enough slack to hide it.
    assert text.count("benchmark:") >= 9


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


def test_the_slurm_profile_parses_and_leaves_only_site_specific_values_unset(
    repo_root: Path,
) -> None:
    """The profile is never submitted here, so it is checked the only honest way.

    Parsing it and asserting that the *only* unset values are the two that genuinely
    cannot be known off-cluster is a stronger statement than "it exists". A placeholder
    that is site-specific is defensible in a file labelled UNTESTED; a placeholder in a
    memory limit is not, and this fails if one is introduced.
    """

    import yaml

    profile = repo_root / "workflow" / "profiles" / "slurm" / "config.yaml"
    settings = yaml.safe_load(profile.read_text())

    assert settings["executor"] == "slurm"
    # The keys are the ones snakemake maps to command line flags; an invented one would
    # be silently ignored by a profile parser and is therefore worse than an error.
    assert set(settings) <= {
        "executor", "jobs", "keep-going", "rerun-incomplete", "latency-wait",
        "default-resources", "set-resources", "cluster-prefix", "cluster-config",
    }
    # This is the setting the resume behaviour in
    # tests/integration/test_workflow_resilience.py depends on, so it is asserted here
    # rather than only in a comment in the profile.
    assert settings["rerun-incomplete"] is True

    assert settings["default-resources"]["slurm_partition"] == "REPLACE_ME"
    assert settings["default-resources"]["slurm_account"] == "REPLACE_ME"
    for key, value in settings["default-resources"].items():
        if key not in ("slurm_partition", "slurm_account"):
            assert isinstance(value, int), f"{key} is left unset; only the two cluster-"\
                                           "specific strings may be placeholders"

    # Every rule the Snakefile declares is either covered by set-resources or left to
    # default-resources. A per-rule block naming a rule that does not exist would be a
    # typo that a profile parser cannot see.
    rules = set(re.findall(r"^rule (\w+):", (repo_root / "workflow" / "Snakefile").read_text(),
                           re.M))
    assert set(settings["set-resources"]) <= rules


def test_the_local_profile_claim_is_true_because_the_harness_passes_it(
    repo_root: Path,
) -> None:
    """`local/config.yaml` says the integration suite runs that profile.  This asserts it.

    The header used to be an unchecked claim: no test passed `--profile`. Either the claim
    or the harness was wrong, and the claim was the easier one to fix by adding the flag,
    so the flag was added rather than the claim withdrawn.
    """

    profile = repo_root / "workflow" / "profiles" / "local" / "config.yaml"
    assert "this profile is what the integration suite runs" in profile.read_text()

    harness = Path(__file__).read_text()
    assert '"--profile"' in harness, (
        "local/config.yaml claims the integration suite runs it; no test passes it"
    )


def test_the_environment_files_are_human_recipes_and_the_claim_is_checked(
    repo_root: Path,
) -> None:
    """Pin `workflow/envs/README.md`'s "nothing consumes these" claim to a test.

    The README says the two environment files are consumed by a person and not by the
    workflow. A README that says so can rot silently the moment someone adds a `conda:`
    directive, so the statement is asserted instead of trusted.
    """

    snakefile = (repo_root / "workflow" / "Snakefile").read_text()
    rules = re.findall(r"^rule (\w+):", snakefile, re.M)
    assert "all" in rules
    for rule in rules:
        if rule == "all":
            continue
        body = snakefile.split(f"rule {rule}:", 1)[1]
        body = re.split(r"^rule \w+:", body, maxsplit=1, flags=re.M)[0]
        assert "conda:" not in body, (
            f"rule {rule} declares a conda environment; workflow/envs/README.md says no "
            f"rule does, and one of the two is now wrong"
        )

    assert "--use-conda" not in (repo_root / "src" / "batchdx" / "cyto" / "cli.py").read_text()
    # This module names the flag in order to assert its absence, so it is excluded from
    # the sweep; the claim is about the harness as a whole, which is tests/.
    others = [
        p.read_text()
        for p in (repo_root / "tests").rglob("*.py")
        if p.resolve() != Path(__file__).resolve()
    ]
    assert not [t for t in others if "--use-conda" in t]

    workflow = (repo_root / ".github" / "workflows" / "batchdx-cyto-tests.yaml").read_text()
    assert "envs/python.yaml" not in workflow and "envs/openproblems.yaml" not in workflow

    readme = (repo_root / "workflow" / "envs" / "README.md").read_text()
    assert "neither is consumed by the workflow" in readme
