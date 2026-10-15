"""Checked R subprocess contracts, without requiring correction packages."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from batchdx.cyto.adapters import r_runner


@pytest.fixture
def interpreter(tmp_path):
    path = tmp_path / "fake-rscript"
    path.write_text(
        f"#!{sys.executable}\n"
        "import json, pathlib, sys, time\n"
        "params = json.loads(pathlib.Path(sys.argv[-2]).read_text())\n"
        "print('partial stdout', flush=True)\n"
        "print('partial stderr', file=sys.stderr, flush=True)\n"
        "time.sleep(params.get('sleep', 0))\n"
        "if 'result' in params:\n"
        "    pathlib.Path(sys.argv[-1]).write_text(params['result'])\n"
        "sys.exit(params.get('exit', 0))\n"
    )
    path.chmod(0o755)
    return str(path)


@pytest.mark.parametrize("result", [
    None, "{", "[]", "null", "{}", '{"status":"failed"}',
    '{"status":"ok","value":NaN}', '{"status":"ok","value":Infinity}',
])
def test_successful_exit_requires_a_valid_result(interpreter, tmp_path, result):
    params = {} if result is None else {"result": result}
    log_path = tmp_path / "logs" / "r.log"
    with pytest.raises(r_runner.RScriptFailed) as caught:
        r_runner.run_r_script("cycombine_correct.R", params, tmp_path / "work",
                              rscript=interpreter, log_path=log_path)
    assert caught.value.returncode == 0
    assert "result.json" in str(caught.value)
    assert "partial stdout" in caught.value.log
    assert "partial stderr" in caught.value.log
    assert log_path.read_text() == caught.value.log


@pytest.mark.parametrize("field", ["warnings", "confounded_nodes"])
@pytest.mark.parametrize("values", [[], ["node 12"], ["node 12", "node 34"]])
def test_list_cardinality_is_preserved(interpreter, tmp_path, field, values):
    payload = {"status": "ok", field: values}
    result = r_runner.run_r_script("cycombine_correct.R", {"result": json.dumps(payload)},
                                   tmp_path, rscript=interpreter)
    assert result.payload == payload


@pytest.mark.parametrize("value", ["one warning", None, [1], {}])
def test_invalid_warning_arrays_are_refused(interpreter, tmp_path, value):
    with pytest.raises(r_runner.RScriptFailed, match="warnings"):
        r_runner.run_r_script("cycombine_correct.R",
                              {"result": json.dumps({"status": "ok", "warnings": value})},
                              tmp_path, rscript=interpreter)


def test_missing_payload_cannot_reuse_a_previous_success(interpreter, tmp_path):
    (tmp_path / "result.json").write_text('{"status":"ok"}')
    with pytest.raises(r_runner.RScriptFailed):
        r_runner.run_r_script("cycombine_correct.R", {}, tmp_path, rscript=interpreter)
    assert not (tmp_path / "result.json").exists()


def test_nonzero_exit_preserves_logs(interpreter, tmp_path):
    with pytest.raises(r_runner.RScriptFailed) as caught:
        r_runner.run_r_script("cycombine_correct.R", {"exit": 7}, tmp_path,
                              rscript=interpreter, log_path=tmp_path / "r.log")
    assert caught.value.returncode == 7
    assert "status 7" in str(caught.value)
    assert "partial stdout" in caught.value.log
    assert (tmp_path / "r.log").read_text() == caught.value.log


def test_timeout_retains_partial_logs_without_inventing_an_exit_code(interpreter, tmp_path):
    with pytest.raises(r_runner.RScriptFailed, match="timed out") as caught:
        r_runner.run_r_script("cycombine_correct.R", {"sleep": 2}, tmp_path,
                              rscript=interpreter, log_path=tmp_path / "r.log", timeout=0.5)
    assert caught.value.returncode is None
    assert "partial stdout" in caught.value.log
    assert "partial stderr" in caught.value.log
    assert "0.5" in caught.value.log
    assert (tmp_path / "r.log").read_text() == caught.value.log


@pytest.mark.parametrize("timeout", [None, 12])
def test_execution_has_a_bounded_default(monkeypatch, tmp_path, timeout):
    monkeypatch.setattr(r_runner, "r_available", lambda _: True)
    seen = []

    def run(command, **kwargs):
        seen.append(kwargs["timeout"])
        Path(command[-1]).write_text('{"status":"ok"}')
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(r_runner.subprocess, "run", run)
    r_runner.run_r_script("cycombine_correct.R", {}, tmp_path, timeout=timeout)
    assert seen == [14400 if timeout is None else timeout]


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan")])
def test_timeout_must_be_positive_and_finite(interpreter, tmp_path, timeout):
    with pytest.raises(ValueError, match="timeout"):
        r_runner.run_r_script("cycombine_correct.R", {}, tmp_path,
                              rscript=interpreter, timeout=timeout)


@pytest.mark.parametrize("probe", [r_runner.r_version, lambda: r_runner.r_package_version("foo")])
def test_availability_probes_are_bounded_and_report_timeout(monkeypatch, probe):
    monkeypatch.setattr(r_runner, "r_available", lambda _: True)

    def run(command, **kwargs):
        assert kwargs["timeout"] == 30
        raise subprocess.TimeoutExpired(command, kwargs["timeout"], output=b"partial")

    monkeypatch.setattr(r_runner.subprocess, "run", run)
    version = probe()
    assert not version.available
    assert "timed out" in version.detail


def test_explicit_unavailability_is_preserved(interpreter, tmp_path):
    payload = {"status": "not_available", "missing_packages": ["foo"]}
    result = r_runner.run_r_script("cytobatchflagr_detect.R",
                                   {"result": json.dumps(payload)}, tmp_path,
                                   rscript=interpreter)
    assert result.payload == payload
