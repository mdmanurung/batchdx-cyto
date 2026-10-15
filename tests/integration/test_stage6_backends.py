"""Stage 6 acceptance gate: backends, output validation, and honest failure.

The ``cycombine`` backend cannot be executed in the environment this build was
developed in, because cyCombine and its Bioconductor dependencies could not be
installed (``docs/blocked.md``).  Its real smoke test is marked
``requires_cycombine`` and skips visibly.  Everything that does not need the R package
— capability declarations, argument rejection, the count-model refusal, the identity
path, output validation and the failure path — is tested here unconditionally.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from batchdx.cyto.adapters import (
    BackendContext,
    BackendError,
    BackendNotAvailable,
    available_backends,
    get_backend,
    resolve_parameters,
)
from batchdx.cyto.adapters.r_runner import r_package_version
from batchdx.cyto.adapters.validate import OutputValidationError, validate_output
from batchdx.cyto.schemas.enums import CorrectionStatus


def _context(frozen, tmp_path, **kwargs):
    config, _, design, sampling = frozen
    return BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                          seed=config.project.seed, **kwargs)


def test_identity_backend_is_numerically_unchanged(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    matrix, result = get_backend("none").run(
        dataset, "uncorrected", {}, _context(frozen, tmp_path)
    )
    assert result.status is CorrectionStatus.SUCCESS
    assert np.array_equal(matrix, dataset.expression())
    assert result.validation["identity_numerically_unchanged"] is True
    assert result.markers_changed == []
    assert result.samples_changed == []
    assert result.output_scale == dataset.processing.scale


def test_every_backend_declares_its_capabilities_explicitly() -> None:
    capabilities = available_backends()
    assert {"none", "cycombine", "cytonorm"} <= set(capabilities)
    shipped = {n: c for n, c in capabilities.items() if not n.startswith("fixture_")}
    for name, caps in shipped.items():
        assert caps.input_scale and caps.output_scale
        assert caps.resource_profile
        if caps.supports_fit_apply and not caps.fit_apply_verified:
            assert "not" in caps.notes.lower() or "verified" in caps.notes.lower(), (
                f"{name} claims fit/apply without verification and does not say so"
            )


def test_cycombine_refuses_the_count_data_route(frozen, tmp_path) -> None:
    backend = get_backend("cycombine")
    with pytest.raises(BackendError) as excinfo:
        resolve_parameters(backend.capabilities, {"method": "ComBat_seq"})
    assert "ComBat_seq" in str(excinfo.value)
    assert "count" in str(excinfo.value)

    # even asking for the right method by name is refused: the route is pinned
    with pytest.raises(BackendError):
        resolve_parameters(backend.capabilities, {"method": "ComBat"})


def test_the_r_script_would_refuse_a_non_combat_method(repo_root: Path) -> None:
    script = (repo_root / "workflow" / "scripts" / "cycombine_correct.R").read_text()
    assert 'identical(as.character(p$method), "ComBat")' in script
    assert "refusing method" in script


def test_unsupported_arguments_are_rejected_not_ignored() -> None:
    backend = get_backend("cycombine")
    with pytest.raises(BackendError, match="does not accept parameter"):
        resolve_parameters(backend.capabilities, {"not_a_real_parameter": 1})

    resolved = resolve_parameters(backend.capabilities, {"xdim": 12})
    assert resolved["xdim"] == 12
    # omitted parameters resolve to the pinned documented defaults and are serialised
    assert resolved["norm_method"] == "scale"
    assert resolved["ydim"] == 8


def test_a_requested_protected_covariate_is_never_quietly_dropped(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    # Point the context at an interpreter that does not exist, so this covers the
    # unavailable path deterministically.  Gating on whether cyCombine happens to be
    # installed made this test a no-op in exactly the environment where the real
    # backend works, which is how an untested path stays untested.
    context = _context(frozen, tmp_path, rscript="Rscript-does-not-exist")
    backend = get_backend("cycombine")
    with pytest.raises(BackendNotAvailable) as excinfo:
        backend.run(dataset, "c", {"covar": "not_a_column"}, context)
    # The absence is reported explicitly and no data comes back.
    assert "not installed" in str(excinfo.value)


def test_a_protected_covariate_is_refused_when_the_backend_is_available(frozen, tmp_path) -> None:
    """The covariate contract holds on the real path, not only when cyCombine is absent.

    A requested protected covariate must never be silently dropped, because dropping it
    changes what the correction is adjusting for.  This is asserted against the
    installed backend; if cyCombine is genuinely unavailable the test says so rather
    than passing vacuously.
    """

    if not r_package_version("cyCombine").available:
        pytest.skip("cyCombine is not installed; see docs/blocked.md")

    _, dataset, _, _ = frozen
    backend = get_backend("cycombine")
    with pytest.raises(BackendError) as excinfo:
        backend.run(dataset, "c", {"covar": "not_a_column"}, _context(frozen, tmp_path))
    message = str(excinfo.value)
    assert "covariate" in message or "column" in message, message


def test_absent_backend_reports_not_available_with_an_actionable_reason(
    frozen, tmp_path
) -> None:
    _, dataset, _, _ = frozen
    # Same reasoning as above: the absent-interpreter case is reproduced rather than
    # waited for, so this path is covered in every environment.
    context = _context(frozen, tmp_path, rscript="Rscript-does-not-exist")
    with pytest.raises(BackendNotAvailable) as excinfo:
        get_backend("cycombine").run(dataset, "c", {}, context)
    message = str(excinfo.value)
    assert "cyCombine" in message
    assert "never substitutes another algorithm" in message


def test_cytonorm_is_registered_as_a_fit_apply_backend(frozen, tmp_path) -> None:
    """CytoNorm is implemented, so this asserts the contract it now has to honour.

    This test previously asserted the pre-implementation state (`fit_apply_verified is
    False` and `run()` raising "not implemented").  CytoNorm 2.0.12 at the pinned commit
    now installs and executes, so that state is gone and asserting it would pin the
    package to being permanently unable to use the one backend whose interface is
    genuinely fit-and-apply.  The behaviour that made this test worth having — that an
    artificial aggregate control is refused — is kept in
    `test_cytonorm_refuses_an_artificial_aggregate_control` below, and the executed
    round trip, event-identity verification and split-fit refusal are covered by
    tests/integration/test_cytonorm_backend.py.
    """

    backend = get_backend("cytonorm")
    assert backend.capabilities.supports_fit_apply is True
    assert "measured_reference_in_every_batch" in backend.capabilities.reference_layouts
    # A transform list is refused as a parameter: the baseline is already transformed,
    # so declaring one would mean applying a second transform.
    assert "transformList" in backend.capabilities.rejected_parameters
    # The verified flag is earned by an executed test, never asserted: the dedicated
    # suite fails the build if it is True while the package is absent.
    if not r_package_version("CytoNorm").available:
        assert backend.capabilities.fit_apply_verified is False
    assert backend.capabilities.fit_apply_verified in (True, False)


def test_cytonorm_refuses_an_artificial_aggregate_control() -> None:
    backend = get_backend("cytonorm")
    with pytest.raises(BackendError, match="aggregate control"):
        resolve_parameters(backend.capabilities, {"aggregate_control": True})


def test_output_validation_catches_a_changed_excluded_marker(frozen) -> None:
    _, dataset, _, _ = frozen
    markers = dataset.marker_ids("correction")
    corrected = dataset.expression()
    excluded = [m for m in dataset.marker_ids() if m not in markers[:2]]
    corrected[:, dataset.feature_index([excluded[0]])[0]] += 0.5

    with pytest.raises(OutputValidationError, match="outside the correction scope"):
        validate_output(dataset, corrected, markers[:2], candidate_id="c")


def test_output_validation_catches_shape_and_finiteness_problems(frozen) -> None:
    _, dataset, _, _ = frozen
    markers = dataset.marker_ids()
    with pytest.raises(OutputValidationError, match="shape"):
        validate_output(dataset, dataset.expression()[:10, :], markers, candidate_id="c")

    broken = dataset.expression()
    broken[0, 0] = np.inf
    with pytest.raises(OutputValidationError, match="non-finite"):
        validate_output(dataset, broken, markers, candidate_id="c")


def test_a_runtime_failure_never_returns_unchanged_data_as_success(
    frozen, tmp_path, example_project: Path
) -> None:
    from batchdx.cyto.pipeline import RunLayout, stage_design, stage_partition, stage_validate
    from batchdx.cyto.schemas.contracts import CorrectionResult

    layout = RunLayout(tmp_path / "run")
    config_path = example_project / "config" / "config.yaml"
    config_path.write_text(
        config_path.read_text().replace(
            "  - id: cycombine_default\n    backend: cycombine\n    parameters: {}\n",
            "  - id: boom\n    backend: fixture_failure\n    parameters: {}\n",
        )
    )
    stage_validate(config_path, example_project, layout)
    stage_design(config_path, example_project, layout)
    stage_partition(config_path, example_project, layout)

    from batchdx.cyto.pipeline import stage_correct, stage_metrics

    stage_correct("boom", config_path, example_project, layout)
    result = CorrectionResult.read(layout.correction_result("boom"))
    assert result.status is CorrectionStatus.FAILED
    assert result.output_path is None
    assert result.validation["produced_output"] is False
    # The `correct` rule declares corrected.h5ad as an output, so a failed candidate
    # writes the unmodified baseline to keep the declared output existing. It must
    # never be mistaken for a correction: the result says FAILED, output_path is
    # None, and metrics read no matrix from a non-SUCCESS candidate.
    assert layout.corrected_h5ad("boom").exists()
    import anndata as ad

    stored = ad.read_h5ad(layout.corrected_h5ad("boom"))
    assert stored.uns.get("candidate_id") == "boom"

    stage_metrics("boom", config_path, example_project, layout)
    from batchdx.cyto.schemas.contracts import MetricResult

    metrics = MetricResult.read(layout.metric_result("boom"))
    assert metrics.status == "ineligible"
    assert metrics.records == []
    assert "not a zero score" in metrics.reason


def test_correction_magnitude_is_reported_by_marker_sample_and_population(
    frozen, tmp_path
) -> None:
    _, dataset, _, _ = frozen
    _, result = get_backend("fixture_offset_removal").run(
        dataset, "offset", {}, _context(frozen, tmp_path)
    )
    assert result.magnitude_by_marker
    assert result.magnitude_by_sample
    assert result.magnitude_by_population
    assert result.markers_changed
    assert result.samples_changed


def test_candidate_identity_is_stable_under_parameter_reordering(frozen, tmp_path) -> None:
    _, dataset, _, _ = frozen
    context = _context(frozen, tmp_path)
    _, first = get_backend("none").run(dataset, "uncorrected", {}, context)
    _, second = get_backend("none").run(dataset, "uncorrected", {}, context)
    assert first.candidate.canonical_signature == second.candidate.canonical_signature
    assert first.output_hash == second.output_hash


@pytest.mark.requires_cycombine
def test_cycombine_smoke_test(frozen, tmp_path) -> None:
    if not r_package_version("cyCombine").available:
        pytest.skip(
            "cyCombine is not installed in this environment; see docs/blocked.md. "
            "The Stage 6 gate for the cycombine backend is recorded as BLOCKED, not passed."
        )
    _, dataset, _, _ = frozen
    matrix, result = get_backend("cycombine").run(
        dataset, "cycombine_default", {}, _context(frozen, tmp_path)
    )
    assert result.status is CorrectionStatus.SUCCESS
    assert matrix.shape == dataset.expression().shape
    assert result.validation["excluded_markers_unchanged"]
    assert not np.array_equal(matrix, dataset.expression())
