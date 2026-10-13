"""Real subprocess boundary failures reach correction, detection and metric records."""

import sys

import pytest

from batchdx.cyto.adapters import BackendContext, BackendError, get_backend, r_runner
from batchdx.cyto.diagnostics import DetectorContext, get_detector
from batchdx.cyto.metrics.native import compute_native_metrics
from batchdx.cyto.schemas.enums import MetricStatus


@pytest.mark.parametrize("mode", ["missing", "malformed", "timeout"])
@pytest.mark.parametrize("consumer", ["correction", "detector", "metrics"])
def test_invalid_r_results_cannot_be_reported_as_success(frozen, tmp_path, monkeypatch,
                                                        mode, consumer):
    interpreter = tmp_path / "fake-rscript"
    interpreter.write_text(
        f"#!{sys.executable}\nimport pathlib, sys, time\n"
        "if '-e' in sys.argv:\n    print('1.0.0')\n    sys.exit(0)\n"
        "print('partial R fixture output', flush=True)\n"
        f"mode = {mode!r}\n"
        "if mode == 'timeout':\n    time.sleep(2)\n"
        "if mode == 'malformed':\n    pathlib.Path(sys.argv[-1]).write_text('{')\n"
    )
    interpreter.chmod(0o755)
    monkeypatch.setattr(r_runner, "DEFAULT_R_TIMEOUT", 0.5)
    config, dataset, design, sampling = frozen
    context = BackendContext(design=design, sampling=sampling, work_dir=str(tmp_path),
                             seed=config.project.seed, rscript=str(interpreter))
    reason = "timed out" if mode == "timeout" else "invalid result.json"
    if consumer == "correction":
        with pytest.raises(BackendError, match=reason) as caught:
            get_backend("cycombine").run(dataset, "broken", {}, context)
        assert "partial R fixture output" in str(caught.value)
    elif consumer == "detector":
        result = get_detector("cycombine_express").run(
            dataset, DetectorContext(run_id="t", design=design, sampling=sampling,
                                     out_dir=str(tmp_path), seed=config.project.seed,
                                     rscript=str(interpreter)),
        )
        assert result.status == "failed"
        assert not result.assessments
        assert reason in result.reason
        assert "partial R fixture output" in "\n".join(result.warnings)
    else:
        baseline = dataset.expression()
        records = compute_native_metrics(dataset, baseline, baseline, sampling, design,
                                         run_id="t", candidate_id="broken", context=context,
                                         work_dir=str(tmp_path / "metrics"))
        assert len(records) == 2
        for record in records:
            assert record.status is MetricStatus.FAILED
            assert record.value is None
            assert reason in record.reason
            assert "partial R fixture output" in record.reason
    logs = list(tmp_path.rglob("*.log"))
    assert len(logs) == 1
    assert reason in logs[0].read_text()


def test_correction_refuses_unavailability_before_reading_any_old_matrix(frozen, tmp_path,
                                                                        monkeypatch):
    from batchdx.cyto.adapters import cycombine
    from batchdx.cyto.adapters.r_runner import RResult
    from batchdx.cyto.schemas.contracts import SoftwareVersion

    monkeypatch.setattr(cycombine, "r_version",
                        lambda *_: SoftwareVersion(name="R", available=True))
    monkeypatch.setattr(cycombine, "r_package_version",
                        lambda name, *_: SoftwareVersion(name=name, available=True))
    monkeypatch.setattr(cycombine, "run_r_script",
                        lambda *a, **kw: RResult(0, "", "", {"status": "not_available"}))

    def unexpected_read(*args):
        pytest.fail("a non-success result must not read an old correction matrix")

    monkeypatch.setattr(cycombine, "read_matrix", unexpected_read)
    config, dataset, design, sampling = frozen
    with pytest.raises(BackendError, match="did not report success"):
        get_backend("cycombine").run(
            dataset, "broken", {}, BackendContext(design=design, sampling=sampling,
                                                   work_dir=str(tmp_path), seed=config.project.seed),
        )
