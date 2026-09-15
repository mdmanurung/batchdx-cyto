"""Stage 0 acceptance gate: the audit is executable, and adapters cite only real functions.

The gate is deliberately strict in one direction: a function name referenced by any
adapter that does not appear in the captured audit output fails the build.  It is
lenient in the other: when the pinned checkouts have not been fetched, the tests skip
with a visible reason rather than passing vacuously.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

LOCK = "config/upstream_lock.yaml"
CHECKOUT_DIRS = {
    "cycombine": "cyCombine",
    "cytonorm": "CytoNorm",
    "cytobatchflagr": "CytoBatchFlagR",
    "openproblems_cyto_batch_integration": "task_cyto_batch_integration",
    "spiscy": "spiscy",
}


@pytest.fixture(scope="module")
def lock(repo_root: Path) -> dict:
    return yaml.safe_load((repo_root / LOCK).read_text(encoding="utf-8"))


def _require_checkouts(upstream_dir: Path) -> None:
    missing = [d for d in CHECKOUT_DIRS.values() if not (upstream_dir / d).exists()]
    if missing:
        pytest.skip(
            f"pinned upstream checkouts are absent ({missing}); run scripts/fetch_upstream.sh"
        )


def test_lock_lists_a_commit_licence_and_strategy_for_every_resource(lock: dict) -> None:
    assert lock["schema_version"] == 1
    for key, spec in lock["resources"].items():
        assert spec.get("commit"), f"{key} has no pinned commit"
        assert len(str(spec["commit"])) == 40, f"{key} is not pinned to a full commit sha"
        assert spec.get("licence"), f"{key} has no recorded licence"
        assert spec.get("adapter_strategy"), f"{key} has no recorded adapter strategy"


def test_fetch_script_pins_the_same_commits_as_the_lock(repo_root: Path, lock: dict) -> None:
    script = (repo_root / "scripts" / "fetch_upstream.sh").read_text(encoding="utf-8")
    for key, spec in lock["resources"].items():
        assert str(spec["commit"]) in script, f"{key}'s pinned commit is not in fetch_upstream.sh"


def test_python_audit_runs_and_reconciles_commits_and_metrics(
    repo_root: Path, upstream_dir: Path, tmp_path: Path
) -> None:
    _require_checkouts(upstream_dir)
    out = tmp_path / "audit_python.json"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "audit_upstream.py"),
         "--upstream-dir", str(upstream_dir), "--out", str(out),
         "--lock", str(repo_root / LOCK)],
        capture_output=True, text=True, check=False, cwd=repo_root,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["failures"] == []
    for key in CHECKOUT_DIRS:
        entry = payload["resources"][key]
        assert entry["commit_matches"], f"{key}: {entry['observed_commit']}"

    op = payload["resources"]["openproblems_cyto_batch_integration"]["metrics"]
    assert op["present"]
    assert len(op["metric_names"]) == 9
    assert "cms" not in op["metric_names"], (
        "the plan named a CellMixS 'cms' metric; it is absent at the pinned commit and "
        "must stay recorded as unavailable"
    )


def test_python_audit_fails_on_a_commit_mismatch(
    repo_root: Path, upstream_dir: Path, tmp_path: Path
) -> None:
    """The gate has to be able to fail, or it is not a gate."""

    _require_checkouts(upstream_dir)
    lock_text = (repo_root / LOCK).read_text(encoding="utf-8")
    tampered = tmp_path / "tampered_lock.yaml"
    tampered.write_text(
        lock_text.replace("fcd42dbbe92bfa6bffa6094bf9b4de536e71c249", "0" * 40),
        encoding="utf-8",
    )
    proc = subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "audit_upstream.py"),
         "--upstream-dir", str(upstream_dir), "--out", str(tmp_path / "x.json"),
         "--lock", str(tampered)],
        capture_output=True, text=True, check=False, cwd=repo_root,
    )
    assert proc.returncode != 0
    assert "COMMIT MISMATCH" in proc.stdout


@pytest.mark.requires_r
def test_r_audit_runs_and_finds_every_adapter_referenced_function(
    repo_root: Path, upstream_dir: Path, tmp_path: Path, lock: dict
) -> None:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not on PATH")
    _require_checkouts(upstream_dir)

    out = tmp_path / "audit_r.json"
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(repo_root / "scripts" / "audit_upstream.R"),
         str(out), str(upstream_dir), str(repo_root / LOCK)],
        capture_output=True, text=True, check=False, cwd=repo_root,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    payload = json.loads(out.read_text(encoding="utf-8"))

    for key in ("cycombine", "cytonorm", "cytobatchflagr"):
        entry = payload["resources"][key]
        declared = lock["resources"][key]["adapter_calls"]
        assert entry["mode"] in ("installed", "source")
        assert not entry.get("missing"), (
            f"{key}: adapter-referenced function(s) {entry['missing']} are absent from "
            "the audited sources"
        )
        for name in declared:
            assert name in entry["signatures"], f"{key}.{name} has no captured signature"
            assert entry["signatures"][name]["signature"].startswith(name)


@pytest.mark.requires_r
def test_audit_output_documents_the_facts_the_adapters_rely_on(
    repo_root: Path, upstream_dir: Path, tmp_path: Path
) -> None:
    """Every upstream premise an adapter encodes must be visible in the audit artifact."""

    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not on PATH")
    _require_checkouts(upstream_dir)

    out = tmp_path / "audit_r.json"
    subprocess.run(
        ["Rscript", "--vanilla", str(repo_root / "scripts" / "audit_upstream.R"),
         str(out), str(upstream_dir), str(repo_root / LOCK)],
        capture_output=True, text=True, check=True, cwd=repo_root,
    )
    payload = json.loads(out.read_text(encoding="utf-8"))
    signatures = payload["resources"]["cycombine"]["signatures"]

    # The adapter pins method="ComBat" because the count route is reachable here.
    assert "ComBat_seq" in signatures["correct_data"]["signature"]
    # The cluster detector's batch-count rule is a property of this function.
    detect = (upstream_dir / "cyCombine" / "R" / "detect_batch_effect.R").read_text()
    assert "length(levels(df$batch)) >= 3 & length(levels(df$batch)) <= 30" in detect

    # CytoNorm genuinely separates fitting from application.
    cytonorm = payload["resources"]["cytonorm"]["signatures"]
    assert "files" in cytonorm["CytoNorm.train"]["arguments"]
    assert "model" in cytonorm["CytoNorm.normalize"]["arguments"]

    # CytoBatchFlagR's 1.5 IQR fence, read from the source rather than assumed.
    functions = (upstream_dir / "CytoBatchFlagR" / "R" / "CytoBatchFlagR_functions.R").read_text()
    assert "iqr <- 1.5" in functions


def test_cytobatchflagr_licence_forbids_vendoring_and_the_adapter_says_so(
    repo_root: Path, lock: dict
) -> None:
    spec = lock["resources"]["cytobatchflagr"]
    assert spec["licence"] == "GPL-3.0"
    assert spec["licence_allows_reuse"] is False
    assert "UNMODIFIED" in spec["adapter_strategy"]

    decision = repo_root / "docs" / "decisions" / "0003-cytobatchflagr-integration.md"
    assert decision.exists(), "the licence decision must be recorded before the adapter"
    text = decision.read_text(encoding="utf-8")
    assert "GPL-3" in text

    adapter = (repo_root / "src" / "batchdx" / "cyto" / "diagnostics"
               / "cytobatchflagr.py").read_text(encoding="utf-8")
    assert "arm's length" in adapter or "arm's-length" in adapter


def test_every_adapter_referenced_r_function_appears_in_the_lock(repo_root: Path, lock: dict) -> None:
    """A guard against an adapter quietly growing a call the audit never checks."""

    sources = {
        "cycombine": [
            repo_root / "src/batchdx/cyto/adapters/cycombine.py",
            repo_root / "workflow/scripts/cycombine_correct.R",
            repo_root / "workflow/scripts/cycombine_detect.R",
        ],
        "cytobatchflagr": [
            repo_root / "workflow/scripts/cytobatchflagr_detect.R",
        ],
    }
    import re

    for key, paths in sources.items():
        declared = set(lock["resources"][key]["adapter_calls"])
        for path in paths:
            text = path.read_text(encoding="utf-8")
            if key == "cycombine":
                used = {
                    name.rstrip(".")
                    for name in re.findall(r"cyCombine::([A-Za-z_.][A-Za-z0-9_.]*)", text)
                }
            else:
                used = {
                    name for name in
                    re.findall(r"\b([a-z_][a-z0-9_]*)\s*\(", text)
                    if name in {"neg_mfi_df_all", "pos_mfi_df_all", "percent_pos_df_all",
                                "iqr_dataframe", "iqr_metric_function", "detect_outlier"}
                }
            unlisted = used - declared
            assert not unlisted, (
                f"{path.name} calls {sorted(unlisted)} which the lock does not declare, so "
                "the Stage 0 audit would never check them"
            )
