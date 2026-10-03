"""The report must not silently drop, cap without saying so, or leak anything.

One defect class is covered here, from several angles: **a thing that was measured, or
that was recorded, and then does not appear on the page**.  A reader of a report cannot
distinguish "did not run" from "ran and said nothing", so every omission this file tests
is a false statement rather than an omission of redundancy:

* a detector that ran and succeeded appears nowhere (it is visible only through the
  failures table, so a *successful* native detector is the one case that is invisible);
* an assessment with ``flag=None`` — every descriptive statistic reported without a
  threshold — is dropped by a ``selectattr('flag')``;
* a table is cut at 250 or 200 rows with no statement that it was cut;
* the marker explorer takes a prefix of the panel in ``markers.tsv`` row order and never
  says which markers it left out;
* a sample identifier is redacted in the samples table and printed in the clear in the
  diagnostics table, and an absolute path in the resolved configuration is printed as it
  is;
* the interval that Stage 8 computes for every record reaches neither the table nor the
  downloadable CSV;
* a candidate run under ``allow_exploratory_override`` leaves no trace anywhere.

The heavy fixture is module-scoped on purpose.  The Stage 10 report file re-runs the
whole pipeline per test and takes twelve minutes; here one run serves every assertion,
and each test that needs to mutate the run works on its own copy of it, so no test can
be affected by another test's edit.
"""

from __future__ import annotations

import base64
import inspect
import json
import re
import shutil
from pathlib import Path

import pytest
import yaml

from batchdx.cyto.examples import materialise_example
from batchdx.cyto.pipeline import (
    RunLayout,
    stage_correct,
    stage_decide,
    stage_design,
    stage_diagnose,
    stage_metrics,
    stage_partition,
    stage_report,
    stage_validate,
)
from batchdx.cyto.reporting import html
from batchdx.cyto.schemas.contracts import (
    DetectorAssessment,
    DetectorResult,
    MethodEligibility,
    MetricResult,
    SoftwareVersion,
)
from batchdx.cyto.schemas.enums import MetricStatus

#: A native detector that the report used to drop on the floor.  Nothing in the
#: ``acquisition_drift``/``composition`` gate mentions it, so before the fix a
#: successful result for it rendered nowhere.
NATIVE_DETECTOR = "cycombine_express"

#: The two statistics a ``selectattr('flag')`` filter discards: one per detector, both
#: emitted without a ``flag`` argument on purpose.
UNFLAGGED_STATISTICS = ("max_block_step_over_mad", "population_frequency")

CSV_HEADER = re.compile(r"data:text/csv;base64,([A-Za-z0-9+/=]+)")


def _csv_bytes(html_text: str) -> bytes:
    match = CSV_HEADER.search(html_text)
    assert match, "the report carries no downloadable metric CSV"
    return base64.b64decode(match.group(1))


def _csv_text(html_text: str) -> str:
    return _csv_bytes(html_text).decode("utf-8")


#: The resolved configuration, as JSON, read back out of the rendered page rather than
#: string-matched: the template emits it through ``tojson``, which escapes ``<`` as
#: ``<``, so a raw substring search for ``<redacted>/…`` would never match.
CONFIG_BLOCK = re.compile(
    r"Resolved configuration</h3>\s*<pre><code>(.*?)</code></pre>", re.S
)


def _printed_config(html_text: str) -> dict:
    match = CONFIG_BLOCK.search(html_text)
    assert match, "the report does not print its resolved configuration"
    return json.loads(match.group(1))


def _all_strings(payload) -> list[str]:
    if isinstance(payload, str):
        return [payload]
    if isinstance(payload, dict):
        return [s for value in payload.values() for s in _all_strings(value)]
    if isinstance(payload, list):
        return [s for value in payload for s in _all_strings(value)]
    return []


def _no_external_assets(html_text: str) -> None:
    """The report is one offline file; nothing here may add a fetch."""

    assert not re.search(r'(?:src|href)\s*=\s*["\'](?:https?:)?//', html_text, re.I)


def _detector_html(html_text: str, name: str) -> str:
    """The rendered block for one detector, from its anchor to the next detector's.

    Anchored on ``id="detector-<name>"`` rather than on the heading text, because the
    detector's own sub-headings share the ``<h3>`` level and would end the block early.
    """

    start = html_text.index(f'<h3 id="detector-{name}">')
    rest = html_text[start + 1 :]
    nxt = rest.find('<h3 id="detector-')
    return html_text[start: start + 1 + nxt] if nxt != -1 else html_text[start:]


# --------------------------------------------------------------------------------------
# One run, shared
# --------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def project(tmp_path_factory) -> Path:
    target = tmp_path_factory.mktemp("completeness")
    materialise_example("paired_reference", target, cells_per_sample=400)
    return target


@pytest.fixture(scope="module")
def run(project: Path) -> dict:
    """Validate, design, partition, both common detectors, the identity candidate.

    Only the identity candidate is corrected, so the fixture needs no R: this file is
    about the report's completeness, not about the backends, and a fixture that depends
    on which R packages happen to be installed cannot assert anything about a detector
    or a metric.
    """

    layout = RunLayout(project / "run")
    config_path = project / "config" / "config.yaml"
    stage_validate(config_path, project, layout)
    stage_design(config_path, project, layout)
    stage_partition(config_path, project, layout)
    for detector in ("acquisition_drift", "composition"):
        stage_diagnose(detector, config_path, project, layout)
    stage_correct("uncorrected", config_path, project, layout)
    stage_metrics("uncorrected", config_path, project, layout)
    stage_decide(config_path, project, layout)
    stage_report(config_path, project, layout)
    return {"project": project, "config": config_path, "layout": layout}


@pytest.fixture
def own_run(run, tmp_path: Path) -> dict:
    """A private copy of the run, for tests that write to it."""

    root = tmp_path / "run"
    shutil.copytree(run["layout"].root, root)
    return {**run, "layout": RunLayout(root)}


@pytest.fixture
def restore_project_tables(project: Path):
    """Put back the two project files a test rewrites, so the next test sees the fixture.

    The example project is module-scoped, so a test that reverses ``markers.tsv`` or
    points ``upstream.checkout_dir`` somewhere private would otherwise decide what every
    later test in the file sees.
    """

    paths = [project / "config" / name for name in ("config.yaml", "markers.tsv")]
    originals = [(p, p.read_text()) for p in paths]
    yield
    for path, text in originals:
        path.write_text(text)


def _rerender(state: dict) -> str:
    stage_report(state["config"], state["project"], state["layout"])
    return state["layout"].report.read_text(encoding="utf-8")


#: The markers the explorer says it chose, in the order it prints them.
SHOWN_MARKERS = re.compile(r"<code>markers\.tsv</code>: <code>(.*?)</code>", re.S)


def _row_cells(fragment: str) -> list[str]:
    """Every cell of one table row, empties preserved so the positions stay honest."""

    parts = re.sub(r"<[^>]+>", "\x01", fragment).split("\x01")
    return [p.strip() for p in parts[1:-1]]


def _metric_table_columns(html_text: str) -> list[str]:
    body = html_text.split('<table id="metricTable">', 1)[1].split("</table>", 1)[0]
    return _row_cells(re.search(r"<tr>(.*?)</tr>", body, re.S).group(1))


def _metric_table_cells(html_text: str) -> list[list[str]]:
    """One list of cell strings per rendered metric row, empties preserved."""

    body = html_text.split('<table id="metricTable">', 1)[1].split("</table>", 1)[0]
    width = len(_metric_table_columns(html_text))
    return [cells for cells in (_row_cells(row) for row in re.findall(r"<tr>(.*?)</tr>", body, re.S))
            if len(cells) == width]



# --------------------------------------------------------------------------------------
# 1.1  every detector, whatever its namespace and status
# --------------------------------------------------------------------------------------


def test_every_detector_result_on_disk_is_rendered_with_its_status(run) -> None:
    layout = run["layout"]
    html_text = layout.report.read_text(encoding="utf-8")

    on_disk = sorted(
        p for p in (layout.root / "diagnostics").glob("*/detector_result.json")
    )
    assert on_disk, "the fixture ran no detector"
    for path in on_disk:
        result = DetectorResult.read(path)
        block = _detector_html(html_text, result.detector)
        assert f'<span class="pill ok">{result.status}</span>' in block or (
            f'{result.status}</span>' in block
        ), f"{result.detector} is not rendered with its status {result.status}"
        assert f"{result.namespace} namespace" in block


def test_a_successful_native_detector_is_rendered_not_only_its_failures(own_run) -> None:
    """The defect, stated as a test: a *successful* native detector used to vanish.

    A non-``ok`` detector still reached the page through the failures table, which is
    why the gate on two detector names went unnoticed: only the success case was
    invisible.  Every field the contract carries is checked, because a detector rendered
    as a bare name and a status is the other half of the same defect.
    """

    layout = own_run["layout"]
    version = "2.0.12-native"
    result = DetectorResult(
        detector=NATIVE_DETECTOR,
        detector_version=version,
        namespace="native",
        run_id="paired_reference_example",
        status="ok",
        reason="per-marker EMD across batch pairs",
        eligibility=MethodEligibility(
            name=NATIVE_DETECTOR, kind="detector", eligible=True,
            reason_code="eligible", reason="all declared requirements met",
            requirements_met={"n_batches >= 2": True},
        ),
        assessments=[
            DetectorAssessment(
                assessment_id="express:B1:CD4", detector=NATIVE_DETECTOR, namespace="native",
                scope="marker_pair", batch="B1", marker="CD4", statistic="emd",
                value=0.83, flag=True, threshold=0.5, status=MetricStatus.OK,
                evidence_class="descriptive", reason="above the declared EMD threshold",
            ),
            DetectorAssessment(
                assessment_id="express:B2:CD3", detector=NATIVE_DETECTOR, namespace="native",
                scope="marker_pair", batch="B2", marker="CD3", statistic="emd",
                value=0.02, status=MetricStatus.OK, evidence_class="descriptive",
                reason="reported without a threshold",
            ),
        ],
        tables={"emd_by_marker": "cycombine_express/emd_by_marker.csv"},
        plots={"mds": "cycombine_express/mds.png"},
        requirements={"n_batches >= 2": True, "marker_panel": "evaluation"},
        parameters={"n_markers": 2, "permutation_seed": 473},
        warnings=["a warning the page must carry"],
        software=[SoftwareVersion(name="cyCombine", version="0.2.5", source="R library")],
    )
    result.write(layout.detector_result(NATIVE_DETECTOR))

    html_text = _rerender(own_run)
    block = _detector_html(html_text, NATIVE_DETECTOR)

    assert f"{NATIVE_DETECTOR}\n" in block or NATIVE_DETECTOR in block
    assert f'{result.status}</span>' in block
    assert "native namespace" in block
    assert version in block, "detector_version is fetched and must be rendered"
    assert "all declared requirements met" in block, "eligibility is fetched and dropped"
    assert "n_batches &gt;= 2" in block or "n_batches >= 2" in block, (
        "requirements are fetched and dropped"
    )
    assert "permutation_seed" in block, "parameters are fetched and dropped"
    assert "cyCombine" in block and "0.2.5" in block, "software is fetched and dropped"
    assert "emd_by_marker.csv" in block, "tables are fetched and dropped"
    assert "mds.png" in block, "plots are fetched and dropped"
    assert "a warning the page must carry" in block
    assert "per-marker EMD across batch pairs" in block


def test_a_native_detector_that_really_ran_is_reported_wherever_it_landed(run) -> None:
    """Not a fixture detector: the real adapter, whatever state it reaches.

    ``stage_diagnose`` writes whatever the adapter produced — ``ok``, ``ineligible`` or
    ``not_available`` — and every one of those is a report the reader needs.  The
    assertion is therefore on the status that actually happened, not on a hoped-for one.
    """

    layout = run["layout"]
    stage_diagnose(NATIVE_DETECTOR, run["config"], run["project"], layout)
    produced = DetectorResult.read(layout.detector_result(NATIVE_DETECTOR))
    assert produced.detector == NATIVE_DETECTOR
    assert produced.namespace == "native"

    html_text = _rerender(run)
    block = _detector_html(html_text, NATIVE_DETECTOR)
    assert f"{produced.status}</span>" in block
    if produced.reason:
        assert produced.reason.split(".")[0] in block
    _no_external_assets(html_text)


# --------------------------------------------------------------------------------------
# 1.2  unflagged assessments, and stated truncations
# --------------------------------------------------------------------------------------


def test_assessments_reported_without_a_flag_are_rendered_in_their_own_table(run) -> None:
    html_text = run["layout"].report.read_text(encoding="utf-8")

    for name in ("acquisition_drift", "composition"):
        block = _detector_html(html_text, name)
        for statistic in UNFLAGGED_STATISTICS:
            if statistic in block:
                assert "Assessments reported without a flag" in block, (
                    f"{name}: {statistic} is rendered somewhere that does not say it "
                    "carries no flag"
                )
    drift = _detector_html(html_text, "acquisition_drift")
    assert "max_block_step_over_mad" in drift, (
        "the drift detector's own headline discontinuity statistic is dropped"
    )
    composition = _detector_html(html_text, "composition")
    assert "population_frequency" in composition, (
        "every per-sample population frequency row is dropped"
    )


def test_unflagged_assessments_are_all_present_not_a_sample_of_them(run) -> None:
    """The count on the page is the count in the artifact, flagged and unflagged."""

    layout = run["layout"]
    html_text = layout.report.read_text(encoding="utf-8")
    for path in sorted((layout.root / "diagnostics").glob("*/detector_result.json")):
        result = DetectorResult.read(path)
        block = _detector_html(html_text, result.detector)
        unflagged = [a for a in result.assessments if not a.flag]
        flagged = [a for a in result.assessments if a.flag]
        assert f"<strong>{len(result.assessments)}</strong> assessment(s)" in block
        assert f"<strong>{len(flagged)}</strong> flagged" in block
        assert f"<strong>{len(unflagged)}</strong> not" in block
        # every distinct statistic in the artifact appears in the rendered block
        for statistic in {a.statistic for a in result.assessments}:
            assert statistic in block, f"{result.detector}: {statistic} is not rendered"


def test_a_truncated_assessment_table_states_how_many_rows_it_dropped(
    own_run, monkeypatch
) -> None:
    monkeypatch.setattr(html, "ASSESSMENT_ROW_LIMIT", 5)
    html_text = _rerender(own_run)
    composition = _detector_html(html_text, "composition")
    assert "Showing 5 of 96 rows" in composition
    assert "the rest are in the artifact on disk" in composition
    # and nothing pretends the table is complete
    assert "Showing all 96 rows" not in composition


def test_a_truncated_joint_structure_table_states_how_many_rows_it_dropped(
    own_run, monkeypatch
) -> None:
    total = sum(
        len(json.loads(path.read_text())["records"])
        for path in own_run["layout"].root.glob("metrics/*/metric_result.json")
    )
    joint = sum(
        1
        for path in own_run["layout"].root.glob("metrics/*/metric_result.json")
        for record in json.loads(path.read_text())["records"]
        if record["metric_name"] == "joint_structure_distance"
    )
    assert joint > 5
    assert total > joint

    monkeypatch.setattr(html, "JOINT_ROW_LIMIT", 5)
    html_text = _rerender(own_run)
    assert f"Showing 5 of {joint} rows" in html_text
    # the full count of joint records exceeds the cap, so the cap really bit
    assert joint > 5


# --------------------------------------------------------------------------------------
# 1.3  the marker caps are declared, and the choice is deterministic
# --------------------------------------------------------------------------------------


def test_the_marker_explorer_declares_what_it_omits(
    own_run, monkeypatch, restore_project_tables
) -> None:
    monkeypatch.setattr(html, "EXPLORER_MARKER_LIMIT", 3)
    monkeypatch.setattr(html, "COMPARISON_MARKER_LIMIT", 2)
    html_text = _rerender(own_run)

    sampling = json.loads((own_run["layout"].sampling_manifest).read_text())
    declared = sorted(sampling["evaluation_markers"])
    assert len(declared) > 3

    shown = SHOWN_MARKERS.search(html_text)
    assert shown is not None
    assert shown.group(1) == ", ".join(declared[:3])
    assert f"Showing 3 of {len(declared)} evaluation" in html_text

    omitted = declared[3:]
    assert f"Not shown here ({len(omitted)}):" in html_text
    for marker in omitted:
        assert marker in html_text
    # the cap is stated as a cap, not as "this is the panel"
    assert "No evaluation marker is omitted" not in html_text
    assert "Showing 2 of 8 evaluation" in html_text


def test_the_markers_shown_do_not_depend_on_the_row_order_of_markers_tsv(
    own_run, monkeypatch, restore_project_tables
) -> None:
    """The source list is in manifest declaration order, so the cap must not inherit it.

    Reversing ``markers.tsv`` is exactly the change that used to decide which markers a
    reader sees.  With the cap declared in sorted order the page is unchanged, and a
    reader of two reports of the same panel sees the same markers.
    """

    monkeypatch.setattr(html, "EXPLORER_MARKER_LIMIT", 3)
    before = _rerender(own_run)
    before_markers = json.loads(
        own_run["layout"].sampling_manifest.read_text()
    )["evaluation_markers"]

    markers_path = own_run["project"] / "config" / "markers.tsv"
    lines = markers_path.read_text().splitlines()
    markers_path.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")

    # The frozen manifest is rebuilt, so the ordering really did change underneath, and
    # it is not the sorted order: if it were, this test could not fail.
    stage_validate(own_run["config"], own_run["project"], own_run["layout"])
    stage_design(own_run["config"], own_run["project"], own_run["layout"])
    stage_partition(own_run["config"], own_run["project"], own_run["layout"])
    reordered = json.loads(own_run["layout"].sampling_manifest.read_text())
    declared = reordered["evaluation_markers"]
    assert declared != sorted(declared), "the fixture did not actually reorder the panel"
    assert declared == before_markers[::-1] or set(declared) == set(before_markers)

    after = _rerender(own_run)
    shown_before, shown_after = SHOWN_MARKERS.search(before), SHOWN_MARKERS.search(after)
    assert shown_before is not None and shown_after is not None
    assert shown_before.group(1) == shown_after.group(1)
    assert shown_after.group(1) == ", ".join(
        sorted(reordered["evaluation_markers"])[:3]
    )


# --------------------------------------------------------------------------------------
# 1.4 / 1.5  redaction: the parameter is gone, the real leaks are closed
# --------------------------------------------------------------------------------------


def test_the_metric_row_projection_takes_no_redactor_parameter() -> None:
    """The parameter lied about its effect, so it is gone rather than kept decorative.

    Kept as a signature assertion rather than a comment, because a parameter that claims
    to redact and does not is exactly the kind of thing a reader of the signature trusts.
    """

    assert "redactor" not in inspect.signature(html._metric_rows).parameters
    assert "redactor" in html._metric_rows.__doc__


def test_a_sample_identifier_is_redacted_in_the_diagnostics_table_too(own_run) -> None:
    """The leak ``_metric_rows``'s dead parameter was standing in for.

    The samples table redacts ``sample_id`` through the redactor, and a detector
    assessment carries the same field, so printing it raw put every sample id of the run
    into a file whose header says the identifiers are redacted.
    """

    layout = own_run["layout"]
    composition = DetectorResult.read(layout.detector_result("composition"))
    sample_id = next(
        a.sample_id for a in composition.assessments
        if a.sample_id and a.statistic == "population_frequency"
    )
    assert sample_id, "the composition detector reported no per-sample frequency row"

    html_text = _rerender(own_run)
    assert sample_id not in html_text, f"raw sample id {sample_id!r} reached the report"
    assert "population_frequency" in html_text, "the row was dropped rather than redacted"


def test_an_absolute_path_in_the_resolved_configuration_is_rewritten(
    own_run, restore_project_tables
) -> None:
    config_path = own_run["config"]
    raw = yaml.safe_load(config_path.read_text())
    outside = "/srv/private/flow/cycombine-checkout"
    raw["upstream"]["checkout_dir"] = outside
    config_path.write_text(yaml.safe_dump(raw, sort_keys=False))

    html_text = _rerender(own_run)
    assert outside not in html_text, "an absolute path reached the report verbatim"

    printed = _printed_config(html_text)
    assert printed["upstream"]["checkout_dir"] == "<redacted>/cycombine-checkout"
    assert printed["input"]["path"] == "<run>/data/preprocessed.h5ad"
    assert printed["input"]["samples"] == "<run>/config/samples.tsv"
    assert printed["evaluation"]["gating_thresholds"] == "<run>/config/thresholds.tsv"
    for value in _all_strings(printed):
        assert not value.startswith("/"), f"{value!r} is an absolute path on the page"
    _no_external_assets(html_text)


def test_a_detectors_own_per_cell_payload_is_pseudonymised_and_declared(own_run) -> None:
    """Found while rendering ``parameters`` for item 1.1: the leak is not only in config.

    ``cycombine_detect.py:100-111`` puts a per-cell ``sample`` column and the absolute
    paths it wrote into ``DetectorResult.parameters``, because the pinned upstream API
    groups by a column literally named ``sample``.  Rendering that field — which the
    report previously did not do at all — would print every sample identifier of the run
    once per evaluation cell, and the machine's directory layout besides.
    """

    layout = own_run["layout"]
    dataset_ids = [
        a.sample_id
        for a in DetectorResult.read(layout.detector_result("composition")).assessments
        if a.sample_id
    ]
    per_cell = [dataset_ids[0]] * 12
    absolute = "/home/somebody/private/diagnostics/cycombine_express"
    DetectorResult(
        detector=NATIVE_DETECTOR, namespace="native", run_id="r", status="ok",
        parameters={
            "mode": "express",
            "out_dir": absolute,
            "sample": per_cell,
            "batch": ["B1"] * 12,
            "seed": 473,
        },
    ).write(layout.detector_result(NATIVE_DETECTOR))

    html_text = _rerender(own_run)
    block = _detector_html(html_text, NATIVE_DETECTOR)

    for sample_id in set(dataset_ids):
        assert sample_id not in html_text, f"raw sample id {sample_id!r} reached the report"
    assert absolute not in html_text
    assert "out_dir" in block and "mode" in block and "express" in block
    assert "elided_values" in block, "a 12-element list was printed rather than declared"
    assert "12" in block
    assert "S-" in block, "the sample column was dropped instead of pseudonymised"
    _no_external_assets(html_text)


def test_a_metric_reason_that_names_a_sample_is_pseudonymised(own_run) -> None:
    """``MetricRecord.reason`` is free text, and the table prints it."""

    from batchdx.cyto.metrics.records import metric_record
    from batchdx.cyto.metrics.registry import POPULATION_FREQUENCY_SHIFT
    from batchdx.cyto.schemas.contracts import SamplingManifest

    layout = own_run["layout"]
    sampling = SamplingManifest.read(layout.sampling_manifest)
    sample_id = next(
        a.sample_id
        for a in DetectorResult.read(layout.detector_result("composition")).assessments
        if a.sample_id
    )
    payload = json.loads(layout.metric_result("uncorrected").read_text())
    payload["records"].append(
        metric_record(
            POPULATION_FREQUENCY_SHIFT, sampling,
            run_id="r", candidate_id="uncorrected", population="B",
            value=0.1, status=MetricStatus.OK, n_samples=2, n_donors=2, n_cells=60,
            coverage=1.0, aggregation_weight=1.0,
            reason=f"sample {sample_id} is the only one in its batch",
        ).model_dump(mode="json")
    )
    layout.metric_result("uncorrected").write_text(json.dumps(payload, indent=2))

    html_text = _rerender(own_run)
    assert sample_id not in html_text
    assert "is the only one in its batch" in html_text
    _csv_text(html_text)
    assert sample_id not in _csv_text(html_text), "the CSV carries the raw sample id too"


# --------------------------------------------------------------------------------------
# 1.6  n_donors on every family
# --------------------------------------------------------------------------------------


DONOR_FAMILIES = (
    "emd_within_population",
    "mad_signed_change",
    "population_frequency_shift",
    "joint_structure_distance",
)


def test_n_donors_is_reported_on_every_metric_family_that_has_donors(run) -> None:
    layout = run["layout"]
    records = [
        record
        for path in layout.root.glob("metrics/*/metric_result.json")
        for record in json.loads(path.read_text())["records"]
    ]
    families = {r["metric_name"] for r in records}
    for family in DONOR_FAMILIES:
        assert family in families, f"{family} was not produced by the fixture at all"
    # A group no sample contributed to has no donor to count; that is a missing value,
    # not a family that forgot to report one, so the check is on the groups that exist.
    scored = [
        r for r in records
        if r["metric_name"] in DONOR_FAMILIES and r["n_samples"]
    ]
    assert scored, "no record in any family carries a sample count"
    for record in scored:
        assert record["n_donors"], (
            f"{record['metric_name']} ({record['population']}, {record['comparison_id']}) "
            f"reports no donor count while {record['n_samples']} sample(s) each declare one"
        )
        assert record["n_donors"] != 0, "a missing donor count is not a zero"
        assert record["n_donors"] <= record["n_samples"]

    # and the column the reader sees is filled, not blank, for those rows
    html_text = layout.report.read_text(encoding="utf-8")
    columns = _metric_table_columns(html_text)
    assert "n_donors" in columns
    index = columns.index("n_donors")
    metric_index = columns.index("metric")
    rows = _metric_table_cells(html_text)
    for family in DONOR_FAMILIES:
        family_rows = [r for r in rows if r[metric_index] == family and r[index]]
        assert family_rows, f"the report's n_donors column is blank on every {family} row"


def test_a_family_with_no_declared_donor_reports_missing_not_zero() -> None:
    from batchdx.cyto.metrics.common import _donors

    donors_of = {"s1": "D1", "s2": None}
    assert _donors(["s1", "s2"], donors_of) == 1
    assert _donors(["s2"], donors_of) is None, (
        "no sample declares a donor: the record must carry None so the page shows an "
        "empty cell, not a zero that reads as a measurement"
    )
    assert _donors([], donors_of) is None


# --------------------------------------------------------------------------------------
# 1.7  the exploratory override reaches the audit trail
# --------------------------------------------------------------------------------------


#: CytoNorm declares that it needs a measured reference in every batch, and the
#: no-reference example has none — an ineligibility whose ``override_allowed`` is true,
#: so it is exactly the branch that used to be a bare ``pass``.
OVERRIDE_CANDIDATE = "cytonorm_exploratory"


def _no_reference_project(target: Path) -> Path:
    materialise_example("no_reference", target, cells_per_sample=200)
    return target


def _override_config(project: Path, allow: bool) -> None:
    config_path = project / "config" / "config.yaml"
    raw = yaml.safe_load(config_path.read_text())
    raw["design"]["allow_exploratory_override"] = allow
    raw["correction"]["candidates"].append(
        {"id": OVERRIDE_CANDIDATE, "backend": "cytonorm", "parameters": {}}
    )
    # An R interpreter that cannot exist keeps the branch fast and environment
    # independent: the override is recorded before the backend is called, so the record
    # under test is the same whether the package is installed or not.
    raw["upstream"]["rscript"] = "Rscript-does-not-exist"
    config_path.write_text(yaml.safe_dump(raw, sort_keys=False))


def _override_run(project: Path, layout: RunLayout, candidates: list[str]) -> None:
    config_path = project / "config" / "config.yaml"
    stage_validate(config_path, project, layout)
    stage_design(config_path, project, layout)
    stage_partition(config_path, project, layout)
    for candidate in candidates:
        stage_correct(candidate, config_path, project, layout)
        stage_metrics(candidate, config_path, project, layout)
    stage_decide(config_path, project, layout)


def test_an_ineligible_backend_is_refused_and_records_no_override(tmp_path: Path) -> None:
    project = _no_reference_project(tmp_path / "no_override")
    _override_config(project, allow=False)
    layout = RunLayout(tmp_path / "run")
    _override_run(project, layout, [OVERRIDE_CANDIDATE])

    correction = json.loads(layout.correction_result(OVERRIDE_CANDIDATE).read_text())
    assert correction["status"] == "ineligible"
    assert not layout.exploratory_overrides.exists()
    decision = json.loads(layout.decision.read_text())
    assert decision["overrides"] == []


def test_an_exploratory_override_is_recorded_in_the_audit_trail_and_the_report(
    tmp_path: Path,
) -> None:
    project = _no_reference_project(tmp_path / "override")
    _override_config(project, allow=True)
    layout = RunLayout(tmp_path / "run")
    _override_run(project, layout, [OVERRIDE_CANDIDATE])

    assert layout.exploratory_overrides.exists()
    recorded = json.loads(layout.exploratory_overrides.read_text())["overrides"]
    assert [entry["candidate_id"] for entry in recorded] == [OVERRIDE_CANDIDATE]
    entry = recorded[0]
    assert entry["backend"] == "cytonorm"
    assert "missing_reference_in_some_batch" in entry["overrode"]
    assert "allow_exploratory_override" in entry["declared_by"]

    # the correction result carries it too, so a reader of that artifact alone sees it
    correction = json.loads(layout.correction_result(OVERRIDE_CANDIDATE).read_text())
    assert correction["status"] != "ineligible", "the override did not let the backend run"
    assert any("EXPLORATORY OVERRIDE" in w for w in correction["warnings"]), correction

    # the decision result is the one the report reads
    decision = json.loads(layout.decision.read_text())
    assert len(decision["overrides"]) == 1
    assert decision["overrides"][0]["candidate_id"] == OVERRIDE_CANDIDATE
    # every value is a string, as the contract declares
    assert all(isinstance(v, str) for v in decision["overrides"][0].values())

    stage_report(project / "config" / "config.yaml", project, layout)
    html_text = layout.report.read_text(encoding="utf-8")
    assert "Exploratory overrides" in html_text
    assert OVERRIDE_CANDIDATE in html_text
    assert "cannot support a recommendation" in html_text
    _no_external_assets(html_text)


def test_re_running_correct_does_not_duplicate_the_override_record(tmp_path: Path) -> None:
    project = _no_reference_project(tmp_path / "rerun")
    _override_config(project, allow=True)
    layout = RunLayout(tmp_path / "run")
    config_path = project / "config" / "config.yaml"
    _override_run(project, layout, [OVERRIDE_CANDIDATE])
    stage_correct(OVERRIDE_CANDIDATE, config_path, project, layout)
    stage_decide(config_path, project, layout)

    recorded = json.loads(layout.exploratory_overrides.read_text())["overrides"]
    assert [entry["candidate_id"] for entry in recorded] == [OVERRIDE_CANDIDATE]
    decision = json.loads(layout.decision.read_text())
    assert len(decision["overrides"]) == 1


# --------------------------------------------------------------------------------------
# 1.8  the interval reaches the table and the CSV
# --------------------------------------------------------------------------------------


#: A prespecified contrast over the two functional markers, on every frozen population.
#: Same shape as the one tests/scientific/test_stage8_preservation.py runs, so the
#: records below are produced by the real Stage 8 code and not restated by hand.
FUNCTIONAL_CONTRAST = dict(
    contrast_id="case_vs_control",
    factor="group",
    level_a="case",
    level_b="control",
    expected_direction="increase_in_a",
    markers=("CD69", "HLADR"),
    n_resamples=200,
    seed=11,
)


def test_the_interval_and_its_method_reach_the_table_and_the_csv(own_run) -> None:
    from batchdx.cyto.config import load_config
    from batchdx.cyto.design import audit_design
    from batchdx.cyto.io import load_dataset
    from batchdx.cyto.metrics import ContrastSpec, compute_preservation_metrics
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels

    project, layout = own_run["project"], own_run["layout"]
    config = load_config(own_run["config"])
    dataset = load_dataset(config, project)
    design = audit_design(
        list(dataset.samples.values()),
        dataset_id=dataset.manifest.dataset_id,
        correction_factor=config.design.correction_batch,
        nested_within=config.design.nested_within,
        protected=config.design.protect,
        donor_factor=config.design.donor,
        minimum_batches=config.design.minimum_batches,
        minimum_batches_origin=config.design.minimum_batches_origin,
        require_shared_reference=config.design.require_shared_reference,
    )
    labels, tier, confidence = resolve_labels(dataset, config, project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)

    preservation = compute_preservation_metrics(
        dataset, dataset.expression(), sampling, design,
        [ContrastSpec(**FUNCTIONAL_CONTRAST)],
        run_id=config.project.id, candidate_id="uncorrected",
    )
    with_interval = [r for r in preservation.records if r.uncertainty_low is not None]
    assert with_interval, "Stage 8 produced no record carrying an interval"
    kinds = {r.uncertainty_kind for r in with_interval}
    assert any(k.startswith("donor_block_bootstrap_percentile") for k in kinds), kinds
    assert any(k.startswith("restricted_permutation_null") for k in kinds), kinds

    MetricResult(
        run_id=preservation.run_id, candidate_id="uncorrected",
        metric_group=preservation.metric_group, records=preservation.records,
        status=preservation.status, reason=preservation.reason,
    ).write(layout.metric_result("uncorrected"))

    html_text = _rerender(own_run)
    csv_text = _csv_text(html_text)

    for column in ("uncertainty_low", "uncertainty_high", "uncertainty_kind"):
        assert column in html_text, f"{column} is not a column of the metric table"
        assert column in csv_text.splitlines()[0], f"{column} is not in the CSV header"
    for kind in kinds:
        assert kind in html_text, f"interval kind {kind} does not reach the page"
        assert kind in csv_text, f"interval kind {kind} does not reach the CSV"
    record = with_interval[0]
    low = f"{record.uncertainty_low:.6g}"
    high = f"{record.uncertainty_high:.6g}"
    assert low in html_text and high in html_text
    assert low in csv_text and high in csv_text
    _no_external_assets(html_text)


def test_a_record_without_an_interval_reports_an_empty_cell_not_a_zero_width_one(
    own_run,
) -> None:
    html_text = _rerender(own_run)
    rows = re.findall(
        r"<tr>((?:(?!</tr>).)*emd_within_population.*?)</tr>", html_text, re.S
    )
    assert rows, "the metric table rendered no EMD row"
    body = re.sub(r"<[^>]+>", "\x00", rows[0]).split("\x00")
    for value in body:
        assert value not in {"0", "0.0", "0.00000"}, (
            "an absent interval must be an empty cell, not a zero-width interval"
        )


# --------------------------------------------------------------------------------------
# the module as a whole still honours its two invariants
# --------------------------------------------------------------------------------------


def test_the_report_is_still_one_offline_self_contained_file(run) -> None:
    html_text = run["layout"].report.read_text(encoding="utf-8")
    assert html_text.startswith("<!doctype html>")
    _no_external_assets(html_text)
    assert "<svg" in html_text
    assert 'id="metricFilter"' in html_text
    assert not re.search(r"<script[^>]+src=", html_text, re.I), "no external script"


def test_the_rendered_page_never_prints_a_bare_json_decode_error(run) -> None:
    """A guard on the guards: the new sections must not throw at render time."""

    html_text = run["layout"].report.read_text(encoding="utf-8")
    for marker in ("Traceback", "UndefinedError", "TemplateSyntaxError", "jinja2"):
        assert marker not in html_text


def test_decision_candidates_all_appear_on_the_page(run) -> None:
    """The traceability rule, restated for the sections this change touched."""

    decision = json.loads(run["layout"].decision.read_text())
    html_text = run["layout"].report.read_text(encoding="utf-8")
    for assessment in decision["assessments"]:
        assert assessment["candidate_id"] in html_text
        assert assessment["state"] in html_text
    for record in (
        record
        for path in run["layout"].root.glob("metrics/*/metric_result.json")
        for record in json.loads(path.read_text())["records"]
    ):
        assert record["metric_name"] in html_text
