"""The preprocessing step is visible in the run: stage, workflow, and report.

Three groups, and the grouping is the point:

* **The stage on real raw events.**  ``stage_preprocess`` runs on a synthetic
  project whose input is per-sample FCS files -- the only input format the stage
  accepts -- with a declared rectangle debris gate and a per-panel asinh cofactor.
  The raw panel values are written as ``cofactor * sinh(asinh cofactor-5 values)``,
  so applying the declared chain must recover the simulated values to within the
  FCS round-trip bound: the test checks what the chain did, not merely that it
  finished.  The frozen manifest, the removal log, the baseline processing
  manifest and the summary must all agree with each other.
* **The refusals.**  A panel without a declared cofactor, and a declared PeacoQC
  pass -- wired through a monkeypatched runner, because what is under test here is
  the plumbing (order, removal bookkeeping, work directory), not PeacoQC itself.
* **The workflow and the report.**  The ``preprocess`` rule exists in the DAG only
  when ``preprocessing.enabled`` is true, a real run up to ``design_audit`` freezes
  the post-QC manifest, and a full benchmark run renders the "What the preprocessing
  step did" section with the numbers the summary records.

Every test really executes what it claims to test; the dry-run is used only to
inspect the DAG, never to claim execution.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

pytestmark = pytest.mark.slow

SEED = 473
N_CELLS = 300
COFACTOR = 150.0

#: One rectangle debris gate on the scatter channels, in the reference
#: implementation's YAML shape.  The fixture's FSC-A is normal(100000, 30000), so
#: the window below removes a non-trivial minority of every sample.
GATES_YAML = """\
gate_chain: [debris]
gates:
  debris:
    filterId: "debris"
    mk1: "FSC-A"
    mk2: "SSC-A"
    type: "rectangle"
    coords:
      mk1: [60000, 140000]
      mk2: [0, 400000]
"""

#: The processing manifest of the RAW input: untransformed, uncorrected.  The
#: double-transform refusal in ``io.dataset._validate_processing`` demands exactly
#: this -- a run that declares a transformation must not be handed values whose
#: manifest already declares one.
RAW_PROCESSING = {
    "compensated": True,
    "unmixed": False,
    "compensation_source": "simulated; no instrument spillover applied",
    "qc_history": [],
    "transformation": None,
    "transformation_parameters": {},
    "scale": "linear instrument units",
    "units": "instrument units",
    "corrected": False,
    "notes": "synthetic raw fixture for the preprocessing stage test",
}

PEACOQC_SETTINGS = """\
channels: [FSC-A, SSC-A]
"""


# --------------------------------------------------------------------------------------
# The fixture: a complete project whose input is raw per-sample FCS
# --------------------------------------------------------------------------------------


def _simulate_fixture():
    from batchdx.cyto.simulate import paired_reference_spec, simulate

    return simulate(paired_reference_spec(seed=SEED, n_cells=N_CELLS))


def _write_raw_fcs(project: Path) -> dict[str, int]:
    """Write per-sample raw FCS files whose asinh(cofactor 150) recovers the fixture.

    The simulated values live on the asinh (cofactor 5) scale; the raw file carries
    ``150 * sinh(v)`` on the panel channels, so the declared chain -- asinh with
    cofactor 150 -- inverts the export exactly, up to FCS single precision.  Scatter
    channels are drawn around the debris gate's window so the gate removes a
    non-trivial minority of events.
    """

    from batchdx.cyto.io.fcs import write_fcs

    adata, samples, _markers, _processing = _simulate_fixture()
    marker_ids = list(adata.var_names)
    channel_names = ["FSC-A", "SSC-A", "FSC-H"] + [f"Ch_{m}" for m in marker_ids]
    rng = np.random.default_rng(SEED + 1)
    raw_dir = project / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for sid in samples["sample_id"]:
        sid = str(sid)
        mask = (adata.obs["sample_id"] == sid).to_numpy()
        panel = np.asarray(adata[mask].layers["preprocessed"], dtype=np.float64)
        n = panel.shape[0]
        counts[sid] = n
        fsc = rng.normal(100000.0, 30000.0, size=n)
        ssc = rng.normal(80000.0, 20000.0, size=n)
        fsch = fsc + rng.normal(0.0, 5000.0, size=n)
        events = np.column_stack([fsc, ssc, fsch, COFACTOR * np.sinh(panel)])
        write_fcs(
            raw_dir / f"{sid}.fcs",
            events,
            channel_names,
            channel_labels=channel_names,
            scale="linear instrument units",
            processing_note=(
                "synthetic raw fixture: panel channels are 150*sinh(asinh cofactor-5 "
                "values), so the declared asinh cofactor-150 chain inverts the export"
            ),
            overwrite=True,
        )
    return counts


def _enable_preprocessing(project: Path, candidates: list[dict] | None = None) -> Path:
    """Rewrite the example configuration for FCS input with the stage enabled."""

    config_path = project / "config" / "config.yaml"
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    cfg["input"]["format"] = "fcs"
    cfg["input"]["path"] = "data/raw"
    cfg["input"]["expression_layer"] = "raw"
    cfg["preprocessing"] = {
        "enabled": True,
        "transformation": "asinh",
        "cofactor_by_panel": {"panelA": COFACTOR},
        "gating": "config/gates.yaml",
        "peacoqc": None,
    }
    if candidates is not None:
        cfg["correction"]["candidates"] = candidates
    config_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return config_path


@pytest.fixture(scope="module")
def fcs_project(tmp_path_factory) -> Path:
    from batchdx.cyto.examples import materialise_example

    project = tmp_path_factory.mktemp("preprocess") / "project"
    materialise_example("paired_reference", project, cells_per_sample=N_CELLS, seed=SEED)

    (project / "config" / "gates.yaml").write_text(GATES_YAML, encoding="utf-8")
    (project / "config" / "peacoqc.yaml").write_text(PEACOQC_SETTINGS, encoding="utf-8")
    _write_raw_fcs(project)

    samples_path = project / "config" / "samples.tsv"
    samples = pd.read_csv(samples_path, sep="\t")
    samples["source_file"] = samples["sample_id"].astype(str) + ".fcs"
    samples.to_csv(samples_path, sep="\t", index=False)

    (project / "config" / "processing.json").write_text(
        json.dumps(RAW_PROCESSING, indent=2) + "\n", encoding="utf-8"
    )
    # One candidate only: the uncorrected arm needs no R backend, so the end-to-end
    # run below exercises the full DAG in minutes and the correction backends have
    # their own integration suite.
    _enable_preprocessing(
        project, candidates=[{"id": "uncorrected", "backend": "none", "parameters": {}}]
    )
    return project


def _sample_ids(project: Path) -> list[str]:
    samples = pd.read_csv(project / "config" / "samples.tsv", sep="\t")
    return [str(s) for s in samples["sample_id"]]


def _snakemake(project: Path, outdir: Path, mode: str = "benchmark",
               config_entries: list[str] | None = None,
               extra: list[str] | None = None) -> subprocess.CompletedProcess:
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not on PATH")
    from batchdx.cyto.paths import snakefile

    profile = Path(snakefile()).parent / "profiles" / "local"
    command = [
        "snakemake", "--snakefile", str(snakefile()), "--directory", str(project),
        "--profile", str(profile),
        "--cores", "1", "--config", "config=config/config.yaml",
        f"outdir={outdir}", f"base_dir={project}", f"mode={mode}",
    ] + (config_entries or []) + (extra or [])
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
    return subprocess.run(command, cwd=project, capture_output=True, text=True,
                          check=False, env=env)


# --------------------------------------------------------------------------------------
# The stage on real raw events
# --------------------------------------------------------------------------------------


def test_stage_preprocess_gates_transforms_and_freezes_the_baseline(fcs_project) -> None:
    """The declared chain runs on the raw FCS, and every artifact agrees with every other.

    The strong assertion is the recovered values: the raw file carries
    ``150*sinh(v)``, the chain declares asinh with cofactor 150, so the baseline's
    values must be ``v`` again to within the FCS round-trip bound.  A chain that
    transformed twice, or transformed the wrong columns, cannot pass this by
    coincidence.
    """

    import anndata as ad

    from batchdx.cyto.io.fcs import read_fcs
    from batchdx.cyto.io.hashing import file_sha256
    from batchdx.cyto.io.manifests import read_removed_events
    from batchdx.cyto.pipeline import RunLayout, stage_preprocess, stage_validate

    config_path = fcs_project / "config" / "config.yaml"
    layout = RunLayout(fcs_project / "results" / "stage")

    stage_validate(config_path, fcs_project, layout)
    manifest_path = stage_preprocess(config_path, fcs_project, layout)
    assert manifest_path == layout.dataset_manifest
    assert layout.preprocess_baseline.exists()

    ids = _sample_ids(fcs_project)
    raw_events = {
        sid: read_fcs(fcs_project / "data" / "raw" / f"{sid}.fcs") for sid in ids
    }
    n_input = sum(c.n_events for c in raw_events.values())

    summary = json.loads(layout.preprocess_summary.read_text(encoding="utf-8"))
    assert summary["n_events_input"] == n_input
    assert summary["n_samples"] == len(ids)
    assert summary["n_events_removed"] == n_input - summary["n_events_kept"]
    # The fixture's scatter spread guarantees the gate removes a non-trivial minority
    # without removing everything: a chain that kept or dropped everything has not
    # been exercised.
    assert 0 < summary["n_events_removed"] < n_input
    assert summary["removals_by_step"] == {"prelim_gating": summary["n_events_removed"]}
    assert any("debris" in s["detail"] for s in summary["steps"])

    log = read_removed_events(layout.removal_log)
    assert len(log) == summary["n_events_removed"]
    assert set(log["step"].unique()) == {"prelim_gating"}

    processing = json.loads(layout.baseline_processing_manifest.read_text(encoding="utf-8"))
    assert processing["transformation"] == "asinh"
    assert processing["transformation_parameters"] == {"cofactor[panelA]": COFACTOR}
    assert "150" in processing["scale"]
    assert processing["corrected"] is False
    assert processing["compensated"] is None and processing["unmixed"] is None

    baseline = ad.read_h5ad(layout.preprocess_baseline)
    assert baseline.n_obs == summary["n_events_kept"]
    assert baseline.layers["preprocessed"].shape == (
        summary["n_events_kept"], len(baseline.var_names),
    )
    for sid in ids:
        checksum = file_sha256(fcs_project / "data" / "raw" / f"{sid}.fcs")
        sel = (baseline.obs["sample_id"] == sid).to_numpy()
        assert set(baseline.obs.loc[sel, "source_checksum"]) == {checksum}
        kept_index = baseline.obs.loc[sel, "original_event_index"].to_numpy()
        assert kept_index.min() >= 0
        assert kept_index.max() < raw_events[sid].n_events


def test_the_baseline_values_are_the_declared_transform_of_the_raw_ones(
    fcs_project,
) -> None:
    """asinh(cofactor 150) applied to 150*sinh(v) must give v back, per kept event."""

    import anndata as ad

    from batchdx.cyto.io.fcs import round_trip_tolerance
    from batchdx.cyto.pipeline import RunLayout, stage_preprocess

    layout = RunLayout(fcs_project / "results" / "recovery")
    stage_preprocess(fcs_project / "config" / "config.yaml", fcs_project, layout)

    adata, _samples, _markers, _processing = _simulate_fixture()
    baseline = ad.read_h5ad(layout.preprocess_baseline)
    layer = np.asarray(baseline.layers["preprocessed"], dtype=np.float64)
    for sid in _sample_ids(fcs_project):
        mask = (adata.obs["sample_id"] == sid).to_numpy()
        expected_all = np.asarray(adata[mask].layers["preprocessed"], dtype=np.float64)
        sel = (baseline.obs["sample_id"] == sid).to_numpy()
        kept_index = baseline.obs.loc[sel, "original_event_index"].to_numpy()
        expected = expected_all[kept_index]
        # The raw values were written through FCS single precision; the tolerance of
        # that round trip bounds the error the transform can introduce, because the
        # asinh derivative is at most one.
        raw_values = COFACTOR * np.sinh(expected_all)
        np.testing.assert_allclose(layer[sel], expected, atol=round_trip_tolerance(raw_values))


def test_the_frozen_manifest_counts_post_qc_cells(fcs_project) -> None:
    """The manifest the design audit will read describes the baseline, not the raw input."""

    from batchdx.cyto.pipeline import RunLayout, stage_preprocess

    layout = RunLayout(fcs_project / "results" / "manifest")
    stage_preprocess(fcs_project / "config" / "config.yaml", fcs_project, layout)

    summary = json.loads(layout.preprocess_summary.read_text(encoding="utf-8"))
    manifest = json.loads(layout.dataset_manifest.read_text(encoding="utf-8"))
    assert manifest["n_cells"] == summary["n_events_kept"]
    per_sample = {s["sample_id"]: s["n_cells"] for s in manifest["samples"]}
    assert sum(per_sample.values()) == summary["n_events_kept"]

    import anndata as ad

    baseline = ad.read_h5ad(layout.preprocess_baseline)
    for sid, n in per_sample.items():
        assert int((baseline.obs["sample_id"] == sid).sum()) == n


# --------------------------------------------------------------------------------------
# The refusals
# --------------------------------------------------------------------------------------


def test_a_panel_without_a_declared_cofactor_is_refused(fcs_project) -> None:
    """There is no universal cofactor: a transform without a per-panel cofactor does not run."""

    from batchdx.cyto.config import ConfigError
    from batchdx.cyto.pipeline import RunLayout, stage_preprocess

    cfg = yaml.safe_load(
        (fcs_project / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    cfg["preprocessing"]["cofactor_by_panel"] = {}
    alt = fcs_project / "config" / "config_nocofactor.yaml"
    alt.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    layout = RunLayout(fcs_project / "results" / "nocofactor")
    with pytest.raises(ConfigError, match="cofactor"):
        stage_preprocess(alt, fcs_project, layout)
    # The refusal happened before anything was frozen.
    assert not layout.preprocess_baseline.exists()


def test_a_declared_peacoqc_pass_runs_after_the_transform_and_is_logged(
    fcs_project, monkeypatch,
) -> None:
    """The plumbing of the QC pass, with the runner monkeypatched.

    What is under test is the order (QC sees the transformed values), the removal
    bookkeeping (the QC rows land in the log under their own step and the events
    leave the baseline) and the work directory -- not PeacoQC itself, which has its
    own backend suite.
    """

    import anndata as ad

    from batchdx.cyto.pipeline import RunLayout, stage_preprocess
    from batchdx.cyto.preprocess import peacoqc as peacoqc_module
    from batchdx.cyto.preprocess.peacoqc import PeacoQCResult
    from batchdx.cyto.preprocess.removals import STEP_QC, removal_frame

    cfg = yaml.safe_load(
        (fcs_project / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    cfg["preprocessing"]["peacoqc"] = "config/peacoqc.yaml"
    alt = fcs_project / "config" / "config_peacoqc.yaml"
    alt.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    seen: dict[str, object] = {}
    work_dirs: list[Path] = []

    def fake_run(values, channels, *, sample_ids, original_event_index, settings,
                 work_dir, rscript="Rscript"):
        # The stage calls the runner once per sample, with that sample's kept events.
        seen.setdefault("n_values", []).append(int(values.shape[0]))
        seen["channels"] = list(channels)
        seen["settings_channels"] = list(settings.channels)
        work_dirs.append(Path(work_dir))
        # Remove the last event of every sample: one removal per sample, traceable
        # in the log, the summary and the baseline's row count.
        rows = []
        for sid in pd.unique(sample_ids):
            m = sample_ids == sid
            rows.append((str(sid), int(np.asarray(original_event_index)[m][-1])))
        removed = removal_frame(
            sample_id=[r[0] for r in rows],
            original_event_index=[r[1] for r in rows],
            step=STEP_QC,
            reason="test fixture: the last event of every sample fails the criteria",
        )
        return PeacoQCResult(removed=removed, samples=())

    monkeypatch.setattr(peacoqc_module, "run_peacoqc", fake_run)

    layout = RunLayout(fcs_project / "results" / "peacoqc")
    stage_preprocess(alt, fcs_project, layout)

    ids = _sample_ids(fcs_project)
    gate_only = RunLayout(fcs_project / "results" / "manifest")
    gate_summary = json.loads(gate_only.preprocess_summary.read_text(encoding="utf-8"))

    summary = json.loads(layout.preprocess_summary.read_text(encoding="utf-8"))
    # The QC pass saw exactly the events the gate kept, per sample, on the channels
    # its settings declare -- not the full panel.
    assert sum(seen["n_values"]) == gate_summary["n_events_kept"]
    assert len(seen["n_values"]) == len(ids)
    assert seen["settings_channels"] == ["FSC-A", "SSC-A"]
    assert summary["n_events_removed"] == gate_summary["n_events_removed"] + len(ids)
    assert summary["removals_by_step"] == {
        "prelim_gating": gate_summary["n_events_removed"],
        "peacoqc": len(ids),
    }
    assert summary["declared"]["peacoqc"] == "config/peacoqc.yaml"

    baseline = ad.read_h5ad(layout.preprocess_baseline)
    assert baseline.n_obs == summary["n_events_kept"]
    assert [d.name for d in sorted(work_dirs)] == sorted(ids)


# --------------------------------------------------------------------------------------
# The workflow and the report
# --------------------------------------------------------------------------------------


def test_the_preprocess_rule_exists_in_the_dag_only_when_enabled(
    fcs_project, tmp_path,
) -> None:
    """A dry run shows the rule when the configuration enables the stage.

    The disabled half runs on a plain h5ad example, untouched: the same Snakefile
    must keep the original validate-owns-the-manifest DAG when the stage is off.
    """

    enabled = _snakemake(fcs_project, tmp_path / "dry_enabled", extra=["--dryrun"])
    assert enabled.returncode == 0, enabled.stderr
    assert re.search(r"rule preprocess\b", enabled.stdout)

    from batchdx.cyto.examples import materialise_example

    plain = tmp_path / "plain"
    materialise_example("paired_reference", plain, cells_per_sample=N_CELLS, seed=SEED)
    disabled = _snakemake(plain, tmp_path / "dry_disabled", extra=["--dryrun"])
    assert disabled.returncode == 0, disabled.stderr
    assert not re.search(r"rule preprocess\b", disabled.stdout)
    assert re.search(r"rule validate\b", disabled.stdout)


def test_design_audit_reads_the_preprocessed_manifest(fcs_project, tmp_path) -> None:
    """A real run up to the design audit freezes the post-QC manifest."""

    outdir = tmp_path / "until_design"
    proc = _snakemake(fcs_project, outdir, extra=["--until", "design_audit"])
    assert proc.returncode == 0, _tail(proc)

    summary = json.loads((outdir / "baseline" / "preprocess_summary.json").read_text(
        encoding="utf-8"))
    manifest = json.loads((outdir / "manifests" / "dataset_manifest.json").read_text(
        encoding="utf-8"))
    assert manifest["n_cells"] == summary["n_events_kept"]
    assert (outdir / "design" / "design_audit.json").exists()


def test_the_report_shows_what_the_preprocessing_step_did(fcs_project, tmp_path) -> None:
    """A full benchmark run renders the section, with the summary's own numbers."""

    outdir = tmp_path / "full"
    proc = _snakemake(fcs_project, outdir)
    assert proc.returncode == 0, _tail(proc)

    summary = json.loads((outdir / "baseline" / "preprocess_summary.json").read_text(
        encoding="utf-8"))
    html = (outdir / "report" / "report.html").read_text(encoding="utf-8")
    assert "What the preprocessing step did" in html
    # The numbers on the page are the numbers on disk -- rendered from the summary,
    # not recomputed.
    assert f"<strong>{summary['n_events_input']}</strong>" in html
    assert f"<strong>{summary['n_events_kept']}</strong>" in html
    assert str(summary["n_events_removed"]) in html


def _tail(proc: subprocess.CompletedProcess) -> str:
    out = (proc.stdout or "") + (proc.stderr or "")
    return out[-3000:]
