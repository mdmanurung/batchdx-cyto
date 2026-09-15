#!/usr/bin/env python3
"""Stage 0 audit, Python side.

Reports the resolved Python environment, checks every pinned checkout against the
commit recorded in ``config/upstream_lock.yaml``, reads each resource's licence from
the checkout itself, and enumerates the OpenProblems metric components that actually
exist at the pinned commit.

The script exits non-zero when a pinned commit does not match, when a metric this
package claims to reconcile against is absent upstream, or when a metric registered in
the ``compatibility`` namespace has not been verified against an upstream
implementation.  Nothing here is inferred from an upstream summary: every fact comes
from a file in the checkout.

Usage:
    python scripts/audit_upstream.py --upstream-dir .upstream --out docs/audit/audit_python.json
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
from importlib import metadata
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

PYTHON_DEPENDENCIES = [
    "numpy", "pandas", "scipy", "anndata", "h5py", "pydantic", "typer", "jinja2",
    "pyyaml", "flowio", "pyarrow", "matplotlib", "snakemake",
]

LICENCE_PATTERNS = [
    (re.compile(r"GNU GENERAL PUBLIC LICENSE\s+Version 3", re.I), "GPL-3.0"),
    (re.compile(r"GNU GENERAL PUBLIC LICENSE\s+Version 2", re.I), "GPL-2.0"),
    (re.compile(r"\bMIT License\b", re.I), "MIT"),
    (re.compile(r"Apache License", re.I), "Apache-2.0"),
    (re.compile(r"BSD 3-Clause", re.I), "BSD-3-Clause"),
]


def git_commit(path: Path) -> str | None:
    if not (path / ".git").exists():
        return None
    try:
        out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=False)
    except OSError:
        return None
    return out.stdout.strip() or None


def detect_licence(path: Path) -> dict:
    # Several R packages ship a short LICENSE stub next to the full text in LICENSE.md,
    # so every candidate is scanned and the first recognised licence wins.
    seen: dict | None = None
    for name in ("LICENSE", "LICENSE.md", "LICENCE", "LICENCE.md", "COPYING"):
        candidate = path / name
        if not candidate.exists():
            continue
        text = candidate.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        for pattern, spdx in LICENCE_PATTERNS:
            if pattern.search(text):
                return {"file": name, "detected": spdx, "first_line": lines[0] if lines else ""}
        if seen is None:
            seen = {"file": name, "detected": "unrecognised",
                    "first_line": lines[0] if lines else ""}
    if seen is not None:
        return seen
    description = path / "DESCRIPTION"
    if description.exists():
        for line in description.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("License:"):
                return {"file": "DESCRIPTION", "detected": line.split(":", 1)[1].strip(),
                        "first_line": line}
    return {"file": None, "detected": "not found", "first_line": ""}


def read_lock(path: Path) -> dict:
    import yaml

    return yaml.safe_load(path.read_text(encoding="utf-8"))


def openproblems_metrics(checkout: Path) -> dict:
    """Enumerate metric components from the Viash configs actually present."""

    import yaml

    metrics_dir = checkout / "src" / "metrics"
    if not metrics_dir.is_dir():
        return {"present": False, "components": {}, "metric_names": []}
    components: dict[str, list[dict]] = {}
    for config in sorted(metrics_dir.glob("*/config.vsh.yaml")):
        try:
            raw = yaml.safe_load(config.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:  # pragma: no cover - defensive
            components[config.parent.name] = [{"name": "<unparsable>", "error": str(exc)}]
            continue
        info = (raw or {}).get("info", {}) or {}
        entries = info.get("metrics", []) or []
        components[config.parent.name] = [{
            "name": e.get("name"),
            "label": e.get("label"),
            "maximize": e.get("maximize"),
            "min": e.get("min"),
            "max": e.get("max"),
            "implementation": sorted(p.name for p in config.parent.glob("script.*")),
        } for e in entries]
    names = [e["name"] for entries in components.values() for e in entries if e.get("name")]
    return {"present": True, "components": components, "metric_names": sorted(names)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-dir", default=".upstream")
    parser.add_argument("--lock", default=str(REPO_ROOT / "config" / "upstream_lock.yaml"))
    parser.add_argument("--out", default=str(REPO_ROOT / "docs" / "audit" / "audit_python.json"))
    args = parser.parse_args()

    upstream = Path(args.upstream_dir)
    lock = read_lock(Path(args.lock))

    environment = {"python": sys.version.split()[0], "platform": platform.platform()}
    dependencies = {}
    for name in PYTHON_DEPENDENCIES:
        try:
            dependencies[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            dependencies[name] = None

    directories = {
        "cycombine": "cyCombine",
        "cytonorm": "CytoNorm",
        "cytobatchflagr": "CytoBatchFlagR",
        "openproblems_cyto_batch_integration": "task_cyto_batch_integration",
        "spiscy": "spiscy",
    }

    failures: list[str] = []
    resources: dict[str, dict] = {}
    for key, spec in lock["resources"].items():
        path = upstream / directories.get(key, key)
        commit = git_commit(path)
        pinned = str(spec.get("commit") or "")
        entry = {
            "repository": spec.get("repository"),
            "checkout": str(path),
            "checkout_present": path.exists(),
            "pinned_commit": pinned,
            "observed_commit": commit,
            "commit_matches": bool(commit) and commit == pinned,
            "declared_licence": spec.get("licence"),
            "observed_licence": detect_licence(path) if path.exists() else None,
            "adapter_strategy": (spec.get("adapter_strategy") or "").strip(),
        }
        if entry["checkout_present"] and not entry["commit_matches"]:
            failures.append(
                f"{key}: checkout is at {commit} but the lock pins {pinned}"
            )
        resources[key] = entry

    op_path = upstream / directories["openproblems_cyto_batch_integration"]
    op = openproblems_metrics(op_path) if op_path.exists() else {
        "present": False, "components": {}, "metric_names": []
    }
    resources["openproblems_cyto_batch_integration"]["metrics"] = op

    declared = lock["resources"]["openproblems_cyto_batch_integration"].get(
        "metrics_present_at_commit", []
    )
    if op["present"]:
        missing = sorted(set(declared) - set(op["metric_names"]))
        extra = sorted(set(op["metric_names"]) - set(declared))
        resources["openproblems_cyto_batch_integration"]["metric_reconciliation"] = {
            "declared_but_absent": missing, "present_but_undeclared": extra,
        }
        if missing:
            failures.append(
                "openproblems: the lock declares metric(s) that are absent at the pinned "
                f"commit: {missing}"
            )

    # Our own registry: nothing may claim to be an upstream compatibility component
    # without having been verified against one.
    from batchdx.cyto.metrics.registry import registered_metrics

    unverified = sorted(
        f"{ns}.{name}" for (ns, name), spec in registered_metrics().items()
        if ns == "compatibility" and not spec.upstream_verified
    )
    if unverified:
        failures.append(
            "metrics registered in the 'compatibility' namespace without upstream "
            f"verification: {unverified}"
        )

    payload = {
        "generated_by": "scripts/audit_upstream.py",
        "environment": environment,
        "python_dependencies": dependencies,
        "upstream_dir": str(upstream),
        "lock": str(args.lock),
        "resources": resources,
        "registered_metrics": sorted(f"{ns}.{name}" for ns, name in registered_metrics()),
        "failures": failures,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("== batchdx-cyto upstream audit (Python) ==")
    print(f"python {environment['python']} on {environment['platform']}")
    for name, version in dependencies.items():
        print(f"   {name}: {version or 'NOT INSTALLED'}")
    print()
    for key, entry in resources.items():
        status = "ok" if entry["commit_matches"] else (
            "MISSING CHECKOUT" if not entry["checkout_present"] else "COMMIT MISMATCH"
        )
        licence = (entry["observed_licence"] or {}).get("detected", "unknown")
        print(f"-- {key} [{status}] licence={licence} "
              f"declared={entry['declared_licence']}")
        print(f"   pinned   {entry['pinned_commit']}")
        print(f"   observed {entry['observed_commit']}")
    if op["present"]:
        print()
        print("-- openproblems metric components at the pinned commit:")
        for component, entries in op["components"].items():
            names = ", ".join(str(e.get("name")) for e in entries)
            print(f"   {component}: {names}")
        print(f"   total metric names: {len(op['metric_names'])}")
    print()
    print(f"audit written to {out}")

    if failures:
        print("\nFAIL:")
        for f in failures:
            print(f"  - {f}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
