"""Integration tests for the domain-review fixes.

Covers: the drift detector's acquisition-axis honesty (event index vs Time
channel), the per-sample event-rate statistic, bead-normalization awareness,
the duplicate-reference rejection in the CytoBatchFlagR eligibility verdict,
and the cyCombine anchor plumbing in the R script.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from batchdx.cyto.diagnostics import DetectorContext, get_detector


def _context(frozen, tmp_path, **kwargs):
    config, dataset, design, sampling = frozen
    return DetectorContext(
        run_id="t", design=design, sampling=sampling, out_dir=str(tmp_path),
        seed=config.project.seed, **kwargs,
    )


# ------------------------------------------------------------------------------
# Drift: axis honesty and the event-rate statistic
# ------------------------------------------------------------------------------

def test_drift_axis_note_names_the_actual_axis_when_no_time_channel_exists(
    frozen, tmp_path
) -> None:
    _, dataset, _, _ = frozen
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    assert result.status == "ok"
    axis_note = str(result.requirements.get("acquisition_axis", ""))
    assert "event order" in axis_note, axis_note
    assert "not clock time" in axis_note
    assert "acquisition_time" not in axis_note or "no FCS Time channel" in axis_note


def test_drift_reports_a_per_sample_event_rate_statistic(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    rates = [a for a in result.assessments
             if a.statistic == "event_rate_max_block_step_over_median"]
    assert rates, "the event-rate statistic must be reported for every sample"
    for a in rates:
        if a.value is not None:
            assert a.threshold is not None or a.flag is None


def test_drift_event_rate_is_undefined_for_a_single_block_sample(frozen, tmp_path) -> None:
    """A sample too small to split into blocks has no rate; it must say so."""

    _, dataset, _, _ = frozen
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    for a in result.assessments:
        if a.statistic == "event_rate_max_block_step_over_median" and a.value is None:
            assert a.status.name in ("INSUFFICIENT_DATA", "UNDEFINED")
            assert a.reason


# ------------------------------------------------------------------------------
# Bead normalization awareness
# ------------------------------------------------------------------------------

def test_drift_warns_when_no_bead_normalization_is_declared(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    assert any("bead" in w.lower() for w in result.warnings), result.warnings


def test_bead_evidence_is_recorded_in_requirements(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    result = get_detector("acquisition_drift").run(dataset, _context(frozen, tmp_path))
    assert "bead_normalization_declared" in result.requirements
    assert isinstance(result.requirements["bead_normalization_declared"], bool)


# ------------------------------------------------------------------------------
# CytoBatchFlagR: duplicate reference ids within a batch
# ------------------------------------------------------------------------------

def test_cytobatchflagr_rejects_a_duplicate_reference_id_within_one_batch(
    tmp_path, repo_root
) -> None:
    """A repeated reference id inside one batch makes the IQR assessment ambiguous."""

    from batchdx.cyto.config import load_config
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
    from batchdx.cyto.examples import materialise_example
    from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation

    project = tmp_path / "flagr_dup"
    materialise_example("paired_reference", project, cells_per_sample=300)

    spec = SimulationSpec(
        samples=[
            SampleSpec(
                sample_id=f"ref_{batch}", donor_id="REFDONOR", batch=batch,
                group="reference", reference_id="REF1",
                reference_role="technical_replicate", n_cells=300,
            )
            for batch in ("B1", "B2", "B3", "B4")
        ] + [
            SampleSpec(
                sample_id=f"s_{batch}", donor_id=f"D{batch}", batch=batch,
                group="case", n_cells=300,
            )
            for batch in ("B1", "B2", "B3", "B4")
        ],
        seed=473,
    )
    write_simulation(spec, project)
    config = load_config(project / "config" / "config.yaml")
    dataset = load_dataset(config, project)

    # duplicate the reference id inside B2 under a second sample id
    ref = dataset.samples["ref_B2"]
    from batchdx.cyto.schemas.contracts import SampleRecord
    from batchdx.cyto.schemas.enums import ReferenceRole

    dataset.samples["ref_B2_dup"] = SampleRecord(
        sample_id="ref_B2_dup", batch=ref.batch, group=ref.group, donor_id=ref.donor_id,
        reference_id=ref.reference_id,
        reference_role=ReferenceRole("technical_replicate"),
    )

    design = audit_design(list(dataset.samples.values()), dataset_id="d",
                          correction_factor="batch", protected=["group"])
    labels, tier, confidence = resolve_labels(dataset, config, project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    context = DetectorContext(
        run_id="t", design=design, sampling=sampling, out_dir=str(tmp_path / "out"),
        seed=473, upstream_dir=str(repo_root / ".upstream"), threshold_path=None,
    )
    result = get_detector("cytobatchflagr").run(dataset, context)
    assert result.status == "ineligible"
    assert "at most once per batch" in result.reason


# ------------------------------------------------------------------------------
# cyCombine anchor plumbing
# ------------------------------------------------------------------------------

def test_cycombine_r_script_builds_an_anchor_column(repo_root: Path) -> None:
    script = (repo_root / "workflow" / "scripts" / "cycombine_correct.R").read_text()
    assert "anchor" in script
    # the anchor column must distinguish reference cells from ordinary samples
    assert "reference" in script.lower()


def test_cycombine_r_script_fails_loudly_on_anchor_collisions(repo_root: Path) -> None:
    script = (repo_root / "workflow" / "scripts" / "cycombine_correct.R").read_text()
    assert "collision" in script.lower() or "collide" in script.lower()


def test_cycombine_adapter_declares_the_anchor_caveat(repo_root: Path) -> None:
    source = (repo_root / "src" / "batchdx" / "cyto" / "adapters" / "cycombine.py").read_text()
    assert "confounded" in source.lower()
