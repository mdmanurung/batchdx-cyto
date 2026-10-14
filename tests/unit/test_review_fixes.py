"""Regression tests for the domain-review fixes.

Covers: the CytoNorm 2.0 control-free eligibility tier, config refusals for
unimplemented QC flags, the asinh-cofactor cross-check, the FCS channel-name
fallback, and Time-channel metadata exposure.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from batchdx.cyto.design import audit_design
from batchdx.cyto.schemas.contracts import SampleRecord
from batchdx.cyto.schemas.enums import ReferenceRole


def _sample(sample_id, batch, group=None, donor=None, reference=None, role="none", **extra):
    return SampleRecord(
        sample_id=sample_id, batch=batch, group=group, donor_id=donor,
        reference_id=reference, reference_role=ReferenceRole(role),
        technical_factors=extra,
    )


def _audit(samples, **kwargs):
    defaults = dict(dataset_id="test", correction_factor="batch", protected=["group"],
                    donor_factor="donor_id", minimum_batches=2,
                    require_shared_reference=False)
    defaults.update(kwargs)
    return audit_design(samples, **defaults)


def _eligible(audit, name):
    return next(e for e in audit.eligibility if e.name == name)


# ------------------------------------------------------------------------------
# CytoNorm 2.0 control-free tier
# ------------------------------------------------------------------------------

def test_cytonorm_control_free_tier_needs_only_two_batches() -> None:
    samples = []
    for batch in ("B1", "B2"):
        samples.append(_sample(f"s_{batch}", batch, group="case", donor=f"D{batch}"))
    entry = _eligible(_audit(samples), "cytonorm_control_free")
    assert entry.eligible


def test_cytonorm_control_free_tier_is_ineligible_with_one_batch() -> None:
    samples = [_sample("s1", "B1", group="case", donor="D1")]
    entry = _eligible(_audit(samples), "cytonorm_control_free")
    assert not entry.eligible


def test_control_free_tier_is_weaker_than_the_reference_tier() -> None:
    """The control-free tier must be eligible where the reference tier is not."""

    samples = []
    for batch in ("B1", "B2"):
        samples.append(_sample(f"s_{batch}", batch, group="case", donor=f"D{batch}"))
    audit = _audit(samples)
    assert _eligible(audit, "cytonorm_control_free").eligible
    assert not _eligible(audit, "cytonorm").eligible


def test_reference_tier_still_dominates_when_a_reference_exists() -> None:
    samples = []
    for batch in ("B1", "B2"):
        samples.append(_sample(f"ref_{batch}", batch, group="reference", reference="REF1",
                               role="technical_replicate"))
        samples.append(_sample(f"s_{batch}", batch, group="case", donor=f"D{batch}"))
    audit = _audit(samples)
    assert _eligible(audit, "cytonorm").eligible
    assert _eligible(audit, "cytonorm_control_free").eligible


# ------------------------------------------------------------------------------
# Config refusals
# ------------------------------------------------------------------------------

def _base_config_kwargs(**over):
    base = dict(
        enabled=True, transformation="asinh",
        cofactor_by_panel={"default_panel": 5},
    )
    base.update(over)
    return base


def test_config_refuses_peacoqc_which_is_not_implemented() -> None:
    from batchdx.cyto.config import PreprocessingConfig

    with pytest.raises(ValueError, match="not implemented"):
        PreprocessingConfig(**_base_config_kwargs(peacoqc=True))


def test_config_refuses_a_gating_setting_which_is_not_implemented() -> None:
    from batchdx.cyto.config import PreprocessingConfig

    with pytest.raises(ValueError, match="not implemented"):
        PreprocessingConfig(**_base_config_kwargs(gating="spiscy"))


def test_config_accepts_unset_qc_flags() -> None:
    from batchdx.cyto.config import PreprocessingConfig

    cfg = PreprocessingConfig(**_base_config_kwargs())
    assert not cfg.peacoqc


# ------------------------------------------------------------------------------
# asinh cofactor cross-check (processing manifest vs declared scale)
# ------------------------------------------------------------------------------

def _processing(scale, transformation="asinh", parameters=None):
    from batchdx.cyto.schemas.contracts import ProcessingManifest

    return ProcessingManifest(
        scale=scale, transformation=transformation,
        transformation_parameters=parameters or {},
    )


def test_scale_cofactor_disagreement_is_rejected(example_project: Path) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.io.dataset import _validate_processing

    config = load_config(example_project / "config" / "config.yaml")
    manifest = _processing("asinh_cofactor5", parameters={"cofactor": 150})
    with pytest.raises(Exception, match="must agree"):
        _validate_processing(manifest, config)


def test_scale_cofactor_agreement_is_accepted(example_project: Path) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.io.dataset import _validate_processing

    config = load_config(example_project / "config" / "config.yaml")
    manifest = _processing("asinh_cofactor5", parameters={"cofactor": 5})
    _validate_processing(manifest, config)  # must not raise


# ------------------------------------------------------------------------------
# FCS channel resolution and Time channel
# ------------------------------------------------------------------------------

def test_fcs_read_reports_channel_key_and_no_time_for_a_synthetic_file(
    example_project: Path, tmp_path: Path
) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.io.export import export_fcs_per_sample
    from batchdx.cyto.io.fcs import read_fcs

    config = load_config(example_project / "config" / "config.yaml")
    dataset = load_dataset(config, example_project)
    out_dir = tmp_path / "fcs"
    paths = list(export_fcs_per_sample(
        dataset, dataset.expression(), out_dir,
        candidate_id="review-fix-test", invert=False,
    ).values())
    raw = read_fcs(paths[0])
    assert raw.channel_key in ("$PnN", "$PnS")
    # the synthetic export carries no Time channel
    assert raw.time_channel is None
    assert raw.time_values is None


def test_time_channel_pattern_matches_the_declared_convention() -> None:
    from batchdx.cyto.io.fcs import TIME_CHANNEL_PATTERN

    assert TIME_CHANNEL_PATTERN.search("Time")
    assert TIME_CHANNEL_PATTERN.search("time")
    assert not TIME_CHANNEL_PATTERN.search("CD45")
    assert not TIME_CHANNEL_PATTERN.search("Runtime_length")