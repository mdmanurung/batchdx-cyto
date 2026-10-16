"""Interrupt a running workflow, then check that resuming it produces the same artifacts.

``rerun-incomplete: true`` is set in both shipped profiles. Nothing in the suite ever
stopped a run, so that setting was a configuration statement with no test behind it, and
the nearest thing that existed -- ``test_stage11_workflow.py`` deleting an output *after* a
completed run -- exercises re-execution, not recovery from a run that stopped while a rule
was executing.

Two signals are tested, because they are not the same event:

``SIGTERM``
    What a scheduler sends. Snakemake handles it gracefully: it stops scheduling, lets the
    job in flight finish, and exits non-zero. No output is left half-written, so this arm
    proves the run is *resumable* but does not exercise ``rerun-incomplete``.
``SIGKILL``
    What a node failure or an OOM kill does. The job in flight dies with no chance to clean
    up, which is the only way to produce a partially written output and therefore the only
    way to test what ``rerun-incomplete`` is for.

How outputs are compared, and why ``benchmarks/`` is excluded, is measured rather than
assumed; the measurement is written down in :func:`_digests`.

The cost is real: the module runs three full benchmark executions and is marked ``slow``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

pytestmark = pytest.mark.slow

#: How long to wait for the interrupted run to get far enough to be worth interrupting.
#: Generous on purpose: exceeding it is a failure with a readable message, whereas a tight
#: bound would be a flake on a loaded machine.
INTERRUPT_DEADLINE_S = 900.0
#: How long to wait for the process to go away after a signal.
EXIT_DEADLINE_S = 600.0


@pytest.fixture(scope="module")
def project(tmp_path_factory) -> Path:
    from batchdx.cyto.examples import materialise_example

    path = tmp_path_factory.mktemp("resilience")
    materialise_example("paired_reference", path, cells_per_sample=250)
    return path


def _command(project: Path, outdir: Path) -> list[str]:
    from batchdx.cyto.paths import snakefile

    # The same profile the rest of the integration suite runs, so `rerun-incomplete` is
    # actually in force here; without it this file would test a setting nobody ships.
    return [
        "snakemake", "--snakefile", str(snakefile()), "--directory", str(project),
        "--profile", str(Path(snakefile()).parent / "profiles" / "local"),
        "--cores", "1",
        "--config", "config=config/config.yaml",
        f"outdir={outdir}", f"base_dir={project}", "mode=benchmark",
    ]


def _env() -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
    return env


def _run(project: Path, outdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(_command(project, outdir), cwd=project, capture_output=True,
                          text=True, check=False, env=_env())


#: JSON keys that record how long a step took.  They differ between any two runs by
#: definition -- the same reason ``benchmarks/`` is excluded -- and asserting on them
#: would be asserting that the machine is idle.  Measured on the CytoNorm fixture:
#: two trainings with the same seed on the same inputs wrote identical corrected FCS
#: bytes and ``r_elapsed_seconds`` values of 11.74 vs 11.69.
_VOLATILE_JSON_KEYS = {"elapsed_seconds", "r_elapsed_seconds"}

#: FlowSOM model files under ``work/cytonorm_train/``.  Measured on the CytoNorm
#: fixture: retraining with the same seed on the same inputs reproduces the corrected
#: FCS bytes exactly but not the model's RDS serialization, so the model file carries
#: no reproducibility signal the corrected outputs do not already carry.  The
#: corrected FCS files themselves stay byte-compared.
_VOLATILE_SUFFIXES = ("CytoNorm_FlowSOM.RDS",)


def _normalize_json(path: Path) -> bytes:
    """JSON bytes with the measured-volatile timing keys removed, recursively."""

    def strip(node: object) -> object:
        if isinstance(node, dict):
            return {k: strip(v) for k, v in node.items() if k not in _VOLATILE_JSON_KEYS}
        if isinstance(node, list):
            return [strip(v) for v in node]
        return node

    normalized = strip(json.loads(path.read_text(encoding="utf-8")))
    return json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _digests(outdir: Path) -> dict[str, str]:
    """Content hashes of every file under ``outdir``, minus the measured-volatile files.

    ``benchmarks/*.tsv`` records how long each rule took, so it differs between any two
    runs by definition and asserting on it would be asserting that the machine is idle.
    Everything else is expected to be byte-identical, and that expectation is measured
    rather than assumed: two uninterrupted runs of this configuration into *different* output
    directories disagree only on the files that record an absolute path
    (``candidates/*/correction_result.json``, ``candidates/*/**/params.json``,
    ``candidates/*/**/result.json``, ``report/report.html``) -- which cannot arise here,
    because the comparison is always against the same output directory -- and on
    ``benchmarks/``.  Three measured exceptions are normalized rather than compared
    raw: the timing keys in JSON artifacts (``elapsed_seconds`` in the R result
    envelopes, ``resources.r_elapsed_seconds`` in the correction results) and the
    FlowSOM model RDS, whose serialization differs between trainings that produce
    byte-identical corrected output.  Nothing else in the tree carries a timestamp,
    hostname or run id.
    """

    digests: dict[str, str] = {}
    for path in sorted(outdir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(outdir).as_posix()
        if rel.startswith("benchmarks/"):
            continue
        if rel.endswith(_VOLATILE_SUFFIXES):
            continue
        if rel.endswith(".json"):
            digests[rel] = hashlib.sha256(_normalize_json(path)).hexdigest()
        else:
            digests[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


def _interrupt_mid_run(project: Path, outdir: Path, sig: int) -> tuple[int, list[str]]:
    """Run, wait until a ``correct`` rule is demonstrably executing, then signal it.

    Returns ``(returncode, the rule log names present at the moment of the signal)``.

    The wait target is ``logs/correct_*.log`` rather than any log: by the time the first
    correction starts, validation, the design audit, partitioning and every detector have
    finished and roughly half the DAG is still to come. Waiting for the first log of any
    kind would signal during ``validate``, before anything interesting had happened.
    """

    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not on PATH")

    proc = subprocess.Popen(_command(project, outdir), cwd=project, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, env=_env(),
                            start_new_session=True)
    logs = outdir / "logs"
    seen: list[str] = []
    try:
        deadline = time.monotonic() + INTERRUPT_DEADLINE_S
        signalled = False
        while time.monotonic() < deadline:
            if proc.poll() is not None:
                break
            if logs.is_dir() and any(logs.glob("correct_*.log")):
                # Let the rule get past Snakemake's job bookkeeping so the signal lands on
                # a running process rather than one still being scheduled.
                time.sleep(0.25)
                if proc.poll() is None:
                    seen = sorted(p.name for p in logs.glob("*.log"))
                    os.killpg(os.getpgid(proc.pid), sig)
                    signalled = True
                    break
            time.sleep(0.05)
        if not signalled:
            pytest.fail(
                f"no logs/correct_*.log appeared within {INTERRUPT_DEADLINE_S:.0f}s (or the "
                f"run finished first), so no rule was demonstrably in flight and the "
                f"interruption this test exists to perform never happened"
            )

        try:
            proc.communicate(timeout=EXIT_DEADLINE_S)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.communicate()
            pytest.fail(f"the workflow ignored signal {sig} and had to be SIGKILLed")
    finally:
        if proc.poll() is None:  # pragma: no cover - only on an unexpected raise
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.communicate()
    return proc.returncode, seen


def _remove_the_downstream_half(outdir: Path) -> None:
    """Delete everything from the correction stage onwards.

    Without this a second run has nothing to do and cannot be interrupted. It is the same
    move ``test_deleting_a_real_artifact_forces_recomputation`` makes, but here it is what
    makes the window for a signal wide enough to land on a rule that is executing.
    """

    for relative in ("candidates", "metrics", "decisions", "report", "status"):
        shutil.rmtree(outdir / relative, ignore_errors=True)


def _assert_resumed_to_the_same_bytes(reference: dict[str, str], outdir: Path) -> None:
    after = _digests(outdir)
    missing = sorted(set(reference) - set(after))
    assert not missing, f"the resumed run did not restore: {missing}"
    changed = sorted(k for k in reference if reference[k] != after[k])
    assert not changed, (
        "the resumed run produced different bytes for artifacts an uninterrupted run had "
        f"already written: {changed}"
    )


@pytest.fixture(scope="module")
def completed(project: Path) -> tuple[Path, dict[str, str]]:
    """One complete benchmark run, and the hashes every test here compares against.

    Module-scoped because it is the expensive part: a second reference run would double the
    cost of this file for no additional coverage.
    """

    outdir = project / "results" / "reference"
    proc = _run(project, outdir)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    reference = _digests(outdir)
    assert reference, "a completed run produced no files to compare against"
    return outdir, reference


def test_a_sigterm_mid_run_leaves_the_run_resumable_and_the_bytes_unchanged(
    project: Path, completed: tuple[Path, dict[str, str]],
) -> None:
    outdir, reference = completed

    _remove_the_downstream_half(outdir)
    returncode, logs_at_signal = _interrupt_mid_run(project, outdir, signal.SIGTERM)

    assert returncode != 0, (
        "a workflow stopped by SIGTERM must not report success; a zero here would mean the "
        "signal never reached it"
    )
    assert logs_at_signal, "no rule log existed at the moment of the signal"
    assert not outdir.joinpath("report", "report.html").exists(), (
        "the run reached its final target despite the signal, so nothing was interrupted"
    )

    rerun = _run(project, outdir)
    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    _assert_resumed_to_the_same_bytes(reference, outdir)


def test_a_sigkill_mid_rule_does_not_leave_an_artifact_the_rerun_trusts(
    project: Path, completed: tuple[Path, dict[str, str]],
) -> None:
    """The arm that exercises ``rerun-incomplete``.

    A SIGTERM run leaves no half-written output, so it can say nothing about whether a
    partially written one is detected. SIGKILL can produce exactly that, and this is the
    only place in the suite where one exists.
    """

    outdir, reference = completed

    _remove_the_downstream_half(outdir)
    returncode, logs_at_signal = _interrupt_mid_run(project, outdir, signal.SIGKILL)

    assert returncode == -signal.SIGKILL, (
        f"expected the workflow to die from SIGKILL, got returncode {returncode}"
    )
    assert logs_at_signal, "no rule log existed at the moment of the signal"

    rerun = _run(project, outdir)
    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    _assert_resumed_to_the_same_bytes(reference, outdir)