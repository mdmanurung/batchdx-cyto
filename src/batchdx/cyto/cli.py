"""Command-line interface.

The analytical commands delegate execution to Snakemake.  There is no second workflow
engine here: ``diagnose``, ``benchmark`` and ``correct`` build a Snakemake invocation
and hand over, so everything a user can do from the CLI they can also do by calling
``snakemake`` directly with the same ``--config`` values.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import typer

from batchdx.cyto._version import __version__
from batchdx.cyto.config import ConfigError, load_config
from batchdx.cyto.paths import config_root, examples_root, snakefile

app = typer.Typer(
    add_completion=False,
    help="Batch-effect diagnosis, correction and evaluation for cytometry.",
    no_args_is_help=True,
)


def _fail(message: str, code: int = 2) -> None:
    typer.secho(message, fg=typer.colors.RED, err=True)
    raise typer.Exit(code)


def _run_snakemake(
    config_path: Path,
    outdir: Path,
    mode: str,
    cores: int,
    base_dir: Path,
    candidate: str | None = None,
    extra: list[str] | None = None,
    dry_run: bool = False,
) -> int:
    if shutil.which("snakemake") is None:
        _fail(
            "snakemake is not on PATH. Install the workflow extra:\n"
            "  pip install 'batchdx-cyto[workflow]'"
        )
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile()),
        "--directory", str(outdir.parent if outdir.parent != Path("") else Path.cwd()),
        "--cores", str(cores),
        "--config",
        f"config={config_path}",
        f"outdir={outdir}",
        f"base_dir={base_dir}",
        f"mode={mode}",
    ]
    if candidate:
        cmd.append(f"candidate={candidate}")
    if dry_run:
        cmd.append("--dry-run")
    cmd.extend(extra or [])

    typer.secho(" ".join(cmd), fg=typer.colors.BRIGHT_BLACK, err=True)
    proc = subprocess.run(cmd, check=False)
    return proc.returncode


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(__version__)


@app.command()
def init(
    target: Path = typer.Argument(..., help="Directory to create the example project in."),
    example: str = typer.Option(
        "paired_reference",
        help="Which shipped example to materialise: paired_reference, no_reference, confounded.",
    ),
    cells_per_sample: int = typer.Option(1500, help="Synthetic cells per sample."),
) -> None:
    """Create a runnable example project, including its data and every referenced file."""

    from batchdx.cyto.examples import materialise_example

    try:
        paths = materialise_example(example, target, cells_per_sample=cells_per_sample)
    except KeyError as exc:
        _fail(str(exc))
    typer.echo(f"created example {example!r} in {target}")
    for name, path in sorted(paths.items()):
        typer.echo(f"  {name}: {path.relative_to(target)}")
    typer.echo("\nnext:")
    typer.echo(f"  batchdx-cyto diagnose --config {target}/config/config.yaml "
               f"--outdir {target}/results/run1")


@app.command()
def validate(
    config: Path = typer.Option(..., "--config", "-c", exists=True, help="Configuration file."),
    base_dir: Path = typer.Option(Path("."), help="Directory paths in the config resolve against."),
    check_input: bool = typer.Option(True, help="Also load and validate the input data."),
) -> None:
    """Strictly validate a configuration, and optionally the input it points at."""

    try:
        cfg = load_config(config)
    except ConfigError as exc:
        _fail(str(exc))
    typer.secho(f"configuration OK (schema_version {cfg.schema_version})", fg=typer.colors.GREEN)

    if check_input:
        from batchdx.cyto.io import DatasetValidationError, load_dataset
        from batchdx.cyto.io.manifests import ManifestError

        try:
            dataset = load_dataset(cfg, base_dir)
        except (DatasetValidationError, ManifestError) as exc:
            _fail(f"input is not usable:\n  {exc}", code=3)
        typer.secho(
            f"input OK: {dataset.n_cells} cells, {len(dataset.markers)} features, "
            f"{len(dataset.samples)} samples, scale {dataset.processing.scale!r}",
            fg=typer.colors.GREEN,
        )


@app.command()
def diagnose(
    config: Path = typer.Option(..., "--config", "-c", exists=True),
    outdir: Path = typer.Option(..., "--outdir", "-o"),
    base_dir: Path = typer.Option(Path(".")),
    cores: int = typer.Option(1),
    dry_run: bool = typer.Option(False, "--dry-run"),
) -> None:
    """Diagnose batch structure and render a diagnosis-only report (no correction)."""

    code = _run_snakemake(config.resolve(), outdir.resolve(), "diagnose", cores,
                          base_dir.resolve(), dry_run=dry_run)
    raise typer.Exit(code)


@app.command()
def benchmark(
    config: Path = typer.Option(..., "--config", "-c", exists=True),
    outdir: Path = typer.Option(..., "--outdir", "-o"),
    base_dir: Path = typer.Option(Path(".")),
    cores: int = typer.Option(1),
    dry_run: bool = typer.Option(False, "--dry-run"),
) -> None:
    """Run every configured candidate, score them, and assess them under the policy."""

    code = _run_snakemake(config.resolve(), outdir.resolve(), "benchmark", cores,
                          base_dir.resolve(), dry_run=dry_run)
    raise typer.Exit(code)


@app.command()
def correct(
    config: Path = typer.Option(..., "--config", "-c", exists=True),
    outdir: Path = typer.Option(..., "--outdir", "-o"),
    candidate: str = typer.Option(
        ..., "--candidate",
        help="Candidate id to apply. Required: a correction is never chosen for you.",
    ),
    base_dir: Path = typer.Option(Path(".")),
    cores: int = typer.Option(1),
    dry_run: bool = typer.Option(False, "--dry-run"),
) -> None:
    """Apply one explicitly selected candidate configuration."""

    try:
        cfg = load_config(config)
    except ConfigError as exc:
        _fail(str(exc))
    known = [c.id for c in cfg.correction.candidates]
    if candidate not in known:
        _fail(f"candidate {candidate!r} is not in the configuration; available: {known}")

    code = _run_snakemake(config.resolve(), outdir.resolve(), "correct", cores,
                          base_dir.resolve(), candidate=candidate, dry_run=dry_run)
    raise typer.Exit(code)


@app.command()
def report(
    config: Path = typer.Option(..., "--config", "-c", exists=True),
    outdir: Path = typer.Option(..., "--outdir", "-o"),
    base_dir: Path = typer.Option(Path(".")),
    cores: int = typer.Option(1),
) -> None:
    """Regenerate the HTML report from existing results, without rerunning correction."""

    code = _run_snakemake(config.resolve(), outdir.resolve(), "report", cores,
                          base_dir.resolve())
    raise typer.Exit(code)


@app.command()
def info(
    as_json: bool = typer.Option(False, "--json", help="Emit machine-readable output."),
) -> None:
    """Show which backends, detectors and metrics this build actually ships."""

    from batchdx.cyto.adapters import available_backends
    from batchdx.cyto.diagnostics import available_detectors
    from batchdx.cyto.metrics import DEFERRED_METRICS, registered_metrics

    payload = {
        "version": __version__,
        "workflow": str(snakefile()),
        "examples": str(examples_root()),
        "shipped_config": str(config_root()),
        "backends": {
            name: {
                "supports_fit_apply": caps.supports_fit_apply,
                "fit_apply_verified": caps.fit_apply_verified,
                "reference_layouts": list(caps.reference_layouts),
                "accepted_parameters": sorted(caps.accepted_parameters),
                "rejected_parameters": caps.rejected_parameters,
                "notes": caps.notes,
            }
            for name, caps in sorted(available_backends().items())
        },
        "detectors": sorted(available_detectors()),
        "metrics": sorted(f"{ns}.{name}" for ns, name in registered_metrics()),
        "deferred_metrics": DEFERRED_METRICS,
    }
    if as_json:
        typer.echo(json.dumps(payload, indent=2))
        return
    typer.echo(f"batchdx-cyto {payload['version']}")
    typer.echo(f"  workflow: {payload['workflow']}")
    typer.echo("  backends: " + ", ".join(payload["backends"]))
    typer.echo("  detectors: " + ", ".join(payload["detectors"]))
    typer.echo("  metrics: " + ", ".join(payload["metrics"]))


@app.command()
def audit(
    out_dir: Path = typer.Option(Path("docs/audit"), help="Where to write the audit artifacts."),
    upstream_dir: Path = typer.Option(
        Path(".upstream"), help="Directory holding the pinned upstream checkouts."
    ),
    rscript: str = typer.Option("Rscript"),
) -> None:
    """Run the Stage 0 upstream audit scripts and store their captured output."""

    # Resolve the audit scripts from the installed package first (the wheel ships
    # them under workflow_assets via hatch force-include), then from the repo
    # layout so an editable checkout keeps working.
    packaged = Path(__file__).resolve().parent / "workflow_assets" / "scripts"
    repo = Path(__file__).resolve().parents[3] / "scripts"
    scripts_dir = packaged if (packaged / "audit_upstream.py").exists() else repo
    py_script = scripts_dir / "audit_upstream.py"
    r_script = scripts_dir / "audit_upstream.R"
    out_dir.mkdir(parents=True, exist_ok=True)

    failures = []
    if py_script.exists():
        proc = subprocess.run(
            [sys.executable, str(py_script), "--out", str(out_dir / "audit_python.json"),
             "--upstream-dir", str(upstream_dir)],
            check=False,
        )
        if proc.returncode != 0:
            failures.append("audit_upstream.py")
    else:
        failures.append(f"{py_script} not found")

    if shutil.which(rscript) and r_script.exists():
        proc = subprocess.run(
            [rscript, "--vanilla", str(r_script), str(out_dir / "audit_r.json"),
             str(upstream_dir)],
            check=False,
        )
        if proc.returncode != 0:
            failures.append("audit_upstream.R")
    else:
        typer.secho(
            f"skipping the R audit: {rscript!r} not on PATH or {r_script} missing",
            fg=typer.colors.YELLOW, err=True,
        )

    if failures:
        _fail("audit failed for: " + ", ".join(failures), code=4)
    typer.secho(f"audit artifacts written to {out_dir}", fg=typer.colors.GREEN)


def main() -> None:
    app()


if __name__ == "__main__":  # pragma: no cover
    main()
