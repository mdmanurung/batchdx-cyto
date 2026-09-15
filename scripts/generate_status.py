#!/usr/bin/env python3
"""Derive ``docs/implementation_status.md`` from an executed test report.

Status is machine-derived on purpose.  Every acceptance gate in ``scripts/gates.yaml``
names the test identifiers that decide it; this script reads a ``pytest --json-report``
artifact and writes the status document from what actually ran.  A gate line without a
corresponding executed test cannot be produced, and a gate carrying a
``blocked_reason`` is never reported as PASSED however many of its tests pass.

Usage:
    pytest --json-report --json-report-file=build/report.json
    python scripts/generate_status.py --report build/report.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

PASSED = "PASSED"
FAILED = "FAILED"
BLOCKED = "BLOCKED"
PARTIAL = "PARTIAL"
DEFERRED = "DEFERRED"
NOT_RUN = "NOT RUN"

BADGE = {
    PASSED: "✅", FAILED: "❌", BLOCKED: "⛔", PARTIAL: "🟡",
    DEFERRED: "⬜", NOT_RUN: "❔",
}


def load_outcomes(report_path: Path) -> dict[str, str]:
    """Map ``tests/<file>::<test>`` to the outcome pytest recorded for it."""

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    outcomes: dict[str, str] = {}
    for test in payload.get("tests", []):
        node = test["nodeid"]
        outcome = test.get("outcome", "unknown")
        # a test that errored during setup counts as failed, not as absent
        if test.get("setup", {}).get("outcome") == "failed":
            outcome = "failed"
        outcomes[node] = outcome
    return outcomes


def classify(gate: dict, outcomes: dict[str, str]) -> tuple[str, dict[str, int], list[str]]:
    tests = gate.get("tests") or []
    counts = {"passed": 0, "failed": 0, "skipped": 0, "not_run": 0}
    details: list[str] = []

    for node in tests:
        outcome = outcomes.get(node)
        if outcome is None:
            counts["not_run"] += 1
            details.append(f"{node} — not executed in this report")
        elif outcome == "passed":
            counts["passed"] += 1
        elif outcome == "skipped":
            counts["skipped"] += 1
            details.append(f"{node} — skipped")
        else:
            counts["failed"] += 1
            details.append(f"{node} — {outcome}")

    if not tests:
        status = DEFERRED
    elif counts["failed"]:
        status = FAILED
    elif counts["not_run"]:
        status = NOT_RUN
    elif gate.get("blocked_reason"):
        status = BLOCKED
    elif gate.get("partial_reason"):
        status = PARTIAL
    elif counts["skipped"]:
        status = PARTIAL
    else:
        status = PASSED

    # A blocked gate never reads as passed, whatever the counts say.
    if gate.get("blocked_reason") and status in (PASSED, PARTIAL):
        status = BLOCKED
    return status, counts, details


def render(gates: list[dict], outcomes: dict[str, str], report_path: Path,
           summary: dict) -> str:
    generated = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines: list[str] = [
        "# Implementation status",
        "",
        "**This file is generated. Do not edit it by hand.**",
        "",
        f"Derived by `scripts/generate_status.py` from `{report_path.name}` on {generated}.",
        "Every line below corresponds to test identifiers that were executed; a claim",
        "without an executed test cannot appear here.",
        "",
        "Regenerate with:",
        "",
        "```",
        "pytest --json-report --json-report-file=build/report.json",
        "python scripts/generate_status.py --report build/report.json",
        "```",
        "",
        "## Legend",
        "",
        "| Badge | Meaning |",
        "|---|---|",
        f"| {BADGE[PASSED]} PASSED | every test for the gate executed and passed |",
        f"| {BADGE[PARTIAL]} PARTIAL | the gate is met in part; the omission is stated |",
        f"| {BADGE[BLOCKED]} BLOCKED | an external obstacle prevents part of the gate from being tested at all |",
        f"| {BADGE[FAILED]} FAILED | at least one test for the gate failed |",
        f"| {BADGE[DEFERRED]} DEFERRED | not implemented in this build |",
        f"| {BADGE[NOT_RUN]} NOT RUN | a named test was not present in this report |",
        "",
        "## Summary by tier",
        "",
    ]

    by_tier: dict[str, list[tuple[dict, str]]] = {}
    statuses: dict[str, tuple[str, dict, list[str]]] = {}
    for gate in gates:
        status, counts, details = classify(gate, outcomes)
        statuses[gate["stage"]] = (status, counts, details)
        by_tier.setdefault(gate["tier"], []).append((gate, status))

    lines += ["| Tier | Stage | Title | Status | Tests passed |", "|---|---|---|---|---|"]
    for tier in sorted(by_tier):
        for gate, status in by_tier[tier]:
            _, counts, _ = statuses[gate["stage"]]
            total = len(gate.get("tests") or [])
            lines.append(
                f"| {tier} | {gate['stage']} | {gate['title']} | "
                f"{BADGE[status]} {status} | {counts['passed']}/{total} |"
            )

    lines += ["", "## Gate detail", ""]
    for gate in gates:
        status, counts, details = statuses[gate["stage"]]
        lines += [
            f"### Stage {gate['stage']} (Tier {gate['tier']}) — {gate['title']}",
            "",
            f"**{BADGE[status]} {status}** — {counts['passed']} passed, "
            f"{counts['failed']} failed, {counts['skipped']} skipped, "
            f"{counts['not_run']} not executed.",
            "",
        ]
        if gate.get("blocked_reason"):
            lines += ["> **Blocked:** " + " ".join(gate["blocked_reason"].split()), ""]
        if gate.get("partial_reason"):
            lines += ["> **Partial:** " + " ".join(gate["partial_reason"].split()), ""]
        if details:
            lines += ["Tests needing attention:", ""]
            lines += [f"* `{d}`" for d in details]
            lines.append("")
        if gate.get("tests"):
            lines += ["<details><summary>All test identifiers for this gate</summary>", ""]
            for node in gate["tests"]:
                outcome = outcomes.get(node, "not executed")
                lines.append(f"* `{node}` — {outcome}")
            lines += ["", "</details>", ""]

    lines += [
        "## Report totals",
        "",
        f"* collected: {summary.get('total', 'unknown')}",
        f"* passed: {summary.get('passed', 0)}",
        f"* failed: {summary.get('failed', 0)}",
        f"* skipped: {summary.get('skipped', 0)}",
        f"* errors: {summary.get('error', 0)}",
        "",
        "Tests outside the gate lists (helpers and additional checks) are counted in the",
        "totals but do not decide a gate.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, help="pytest --json-report output")
    parser.add_argument("--gates", default=str(REPO_ROOT / "scripts" / "gates.yaml"))
    parser.add_argument("--out", default=str(REPO_ROOT / "docs" / "implementation_status.md"))
    args = parser.parse_args()

    report_path = Path(args.report)
    if not report_path.exists():
        print(f"test report not found: {report_path}")
        print("run: pytest --json-report --json-report-file=build/report.json")
        return 2

    gates = yaml.safe_load(Path(args.gates).read_text(encoding="utf-8"))["gates"]
    outcomes = load_outcomes(report_path)
    summary = json.loads(report_path.read_text(encoding="utf-8")).get("summary", {})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(gates, outcomes, report_path, summary), encoding="utf-8")
    print(f"wrote {out}")

    failed = [g["stage"] for g in gates if classify(g, outcomes)[0] == FAILED]
    if failed:
        print(f"gates with failing tests: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
