"""Render the unified offline HTML report from validated structured artifacts.

The renderer reads only what is on disk.  It never recomputes a score, so every number
in the report is traceable to a record, and regenerating the report cannot trigger a
correction.
"""

from __future__ import annotations

import base64
import copy
import csv
import io
import json
from pathlib import Path

import numpy as np
from jinja2 import Environment, FileSystemLoader, select_autoescape

from batchdx.cyto._version import __version__
from batchdx.cyto.config import load_config
from batchdx.cyto.diagnostics.composition import frequency_table
from batchdx.cyto.io import load_dataset
from batchdx.cyto.io.dataset import SAMPLE_ID
from batchdx.cyto.metrics.registry import (
    DEFERRED_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    registered_metrics,
)
from batchdx.cyto.paths import templates_root
from batchdx.cyto.reporting import plots
from batchdx.cyto.reporting.redact import Redactor
from batchdx.cyto.schemas.contracts import (
    CorrectionResult,
    DecisionResult,
    DesignAudit,
    DetectorResult,
    MetricRecord,
    MetricResult,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import MetricStatus

#: How many evaluation markers the batch explorer draws, and how many of those the
#: before/after comparison draws.  Both are *declared* caps, not silent ones: the page
#: states how many of the panel's evaluation markers it shows and names the ones it
#: omits, because an unstated cap reads as "this is all of it".
EXPLORER_MARKER_LIMIT = 8
COMPARISON_MARKER_LIMIT = 3

#: Rows shown before a diagnostics or preservation table says so.  The full record set
#: is always on disk; the cap is only about what the page has to scroll.
ASSESSMENT_ROW_LIMIT = 250
JOINT_ROW_LIMIT = 200

#: Longest free-form list printed from a detector artifact.  A detector legitimately
#: records one entry per frozen evaluation cell — ``cycombine_express`` carries the
#: per-cell batch, label and sample as the pinned API requires — and printing that is
#: both unreadable and a per-cell identifier dump, so the list is replaced by its length.
PARAMETER_LIST_LIMIT = 8


def _is_path_text(text: str) -> bool:
    """Whether a whole string is a filesystem path rather than prose that mentions one."""

    return text.startswith(("/", "~", "./", "../", ".\\")) or "://" in text


def _redact_free_form(value, redactor: Redactor, known_ids: frozenset[str]):
    """Redact identifiers and paths inside a detector's own free-form payload.

    ``DetectorResult.parameters``, ``requirements``, ``tables``, ``plots`` and
    ``warnings`` are declared ``dict[str, Any]`` / ``list[str]``, so they can carry
    anything a detector chose to record — and one of them does: the pinned
    ``cycombine_express`` API takes a per-cell ``sample`` column, so its parameters hold
    one raw sample identifier per evaluation cell plus the absolute paths it wrote.  The
    rule applied here is declared rather than guessed:

    * a run identifier occurring anywhere in a string is replaced by that identifier's
      pseudonym, so it matches the samples table;
    * a string that *is* a path is rewritten by :meth:`Redactor.path`;
    * a list longer than :data:`PARAMETER_LIST_LIMIT` is replaced by its length and its
      first entries, so the omission is visible on the page;
    * anything else is printed as the detector wrote it.

    The known limit: a path *inside* prose — a captured R log, say — is not rewritten,
    because rewriting a log would falsify it.  The page says so where it happens.
    """

    if isinstance(value, str):
        text = value
        for identifier in known_ids:
            if identifier in text:
                text = text.replace(identifier, redactor.identifier(identifier, "S"))
        return redactor.path(text) if _is_path_text(text) else text
    if isinstance(value, dict):
        return {k: _redact_free_form(v, redactor, known_ids) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        if len(value) > PARAMETER_LIST_LIMIT:
            return {
                "elided_values": len(value),
                "first_values": [
                    _redact_free_form(v, redactor, known_ids) for v in value[:2]
                ],
                "note": (
                    f"a list of {len(value)} values, of which "
                    f"{len(value) - 2} are not shown; the whole array is in the detector "
                    "result on disk"
                ),
            }
        return [_redact_free_form(v, redactor, known_ids) for v in value]
    return value

#: The resolved configuration keys whose value is a filesystem path.  These are the only
#: values ``config.resolved_dict()`` can carry that leak the producing machine, and they
#: are enumerated rather than pattern-matched so that a newly added key is not silently
#: printed in the clear: adding a path-valued key means adding it here.
CONFIG_PATH_KEYS: tuple[str, ...] = (
    "input.path",
    "input.samples",
    "input.markers",
    "input.processing_manifest",
    "evaluation.labels",
    "evaluation.gating_thresholds",
    "upstream.lock",
    "upstream.checkout_dir",
)


def _redact_config_paths(payload: dict, redactor: Redactor) -> dict:
    """The resolved configuration with every path-valued key rewritten by the redactor.

    The report prints the whole resolved configuration, so an unrewritten
    ``input.path`` puts an absolute path — and with it the directory layout of the
    machine, and often the operator's name — into a file whose whole purpose is to be
    circulated.  The rewrite is a copy: ``config.resolved_dict()`` is what the
    fingerprint was computed from and is not mutated here.
    """

    out = copy.deepcopy(payload)
    for dotted in CONFIG_PATH_KEYS:
        section, _, leaf = dotted.partition(".")
        block = out.get(section)
        if isinstance(block, dict) and isinstance(block.get(leaf), str):
            block[leaf] = redactor.path(block[leaf])
    return out


def _data_uri_csv(rows: list[dict], columns: list[str]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    encoded = base64.b64encode(buffer.getvalue().encode()).decode()
    return f"data:text/csv;base64,{encoded}"


def _metric_rows(records: list[MetricRecord]) -> list[dict]:
    """Flatten metric records into the table rows and the downloadable CSV.

    No redactor is involved, and that is a decision rather than an omission.  A
    ``MetricRecord`` carries no sample, donor or reference identifier: its identifying
    fields are ``candidate_id``, ``marker``, ``population`` and ``comparison_id``, and
    ``comparison_id`` is a derived key over batch and factor-level names that the design
    table and the samples table already publish in the clear.  Redacting it here would
    make this table disagree with those two rather than make it safer.  The rows that
    *do* carry a redacted identifier are projected with the redactor at their own site:
    detector assessment sample ids in :func:`_detector_views`, and the path-valued
    configuration keys in :func:`_redact_config_paths`.
    """

    rows = []
    for r in records:
        rows.append({
            "candidate_id": r.candidate_id,
            "namespace": r.metric_namespace,
            "metric": r.metric_name,
            "version": r.metric_version,
            "stage": r.stage,
            "marker": r.marker or "",
            "population": r.population or "",
            "comparison": r.comparison_id or "",
            "value": "" if r.value is None else f"{r.value:.6g}",
            # The interval is reported next to the value, with the method that produced
            # it, so a reader can tell a bootstrap percentile interval from the central
            # interval of a restricted-permutation null.  An absent interval is an empty
            # cell, never a zero-width interval at the value.
            "uncertainty_low": "" if r.uncertainty_low is None else f"{r.uncertainty_low:.6g}",
            "uncertainty_high": "" if r.uncertainty_high is None else f"{r.uncertainty_high:.6g}",
            "uncertainty_kind": r.uncertainty_kind or "",
            "direction": str(r.direction),
            "scale": r.scale,
            "status": str(r.status),
            "reason": r.reason or "",
            "n_samples": r.n_samples if r.n_samples is not None else "",
            "n_donors": r.n_donors if r.n_donors is not None else "",
            "n_cells": r.n_cells if r.n_cells is not None else "",
            "coverage": "" if r.coverage is None else f"{r.coverage:.3f}",
            "weight": f"{r.aggregation_weight:.3f}",
            "label_tier": str(r.label_tier) if r.label_tier else "",
            "label_confidence": r.label_confidence or "",
            "representation": (
                f"{r.representation.space}/{r.representation.marker_mask}/"
                f"{r.representation.scaling}/{r.representation.distance}/k={r.representation.k}"
                if r.representation else ""
            ),
        })
    return rows


def _shown(total: int, limit: int) -> str:
    """The declared cap as text: how many rows are on the page out of how many exist."""

    if total <= limit:
        return f"Showing all {total} row{'' if total == 1 else 's'}."
    return f"Showing {limit} of {total} rows; the rest are in the artifact on disk."


def _scrub_metric_reasons(
    rows: list[dict], redactor: Redactor, known_ids: frozenset[str]
) -> list[dict]:
    """Put the same redaction on a metric row's free-text reason.

    Kept out of :func:`_metric_rows` on purpose.  That function is a projection of the
    record and nothing else — which is why it takes no redactor — so the one field that
    is not the record's own value is scrubbed here, where the run's identifier set is
    known.  ``MetricRecord.reason`` is free text written by a metric, and a metric that
    names a sample in its reason must not print it in the clear.
    """

    return [
        {**row, "reason": _redact_free_form(row["reason"], redactor, known_ids)}
        if row["reason"] else row
        for row in rows
    ]


def _detector_views(
    detectors: list[DetectorResult], redactor: Redactor, known_ids: frozenset[str]
) -> list[dict]:
    """Every detector that produced a result, with every field it carries.

    Nothing here is filtered on the detector's name or namespace.  A gate that admits
    only some detectors is the same failure this package exists to prevent, just in the
    renderer: a native detector that ran and succeeded would be reported nowhere, and a
    reader could not tell the difference between "did not run" and "ran and said nothing".

    Assessments are split into flagged and unflagged rather than filtered to the flagged
    ones.  ``flag`` is ``None`` on every descriptive statistic a detector reports
    without a threshold — the drift detector's own headline discontinuity statistic and
    every per-sample population frequency are two of them — so a ``selectattr('flag')``
    drops the measurements and keeps only the alarms.
    """

    def free(value):
        return _redact_free_form(value, redactor, known_ids)

    views: list[dict] = []
    for d in detectors:
        flagged, unflagged = [], []
        for a in d.assessments:
            row = {
                "sample_id": redactor.identifier(a.sample_id, "S") if a.sample_id else "",
                "batch": a.batch or "",
                "marker": a.marker or "",
                "population": a.population or "",
                "statistic": a.statistic,
                "value": "" if a.value is None else f"{a.value:.4g}",
                "threshold": "" if a.threshold is None else f"{a.threshold:.4g}",
                "status": str(a.status),
                "evidence_class": a.evidence_class,
                "reason": free(a.reason) if a.reason else "",
            }
            (flagged if a.flag else unflagged).append(row)
        views.append({
            "detector": d.detector,
            "namespace": d.namespace,
            "status": d.status,
            "reason": free(d.reason) if d.reason else "",
            "detector_version": d.detector_version or "not declared",
            "warnings": [free(w) for w in d.warnings],
            "eligibility": d.eligibility,
            "requirements": free(d.requirements),
            "parameters": free(d.parameters),
            "tables": free(d.tables),
            "plots": free(d.plots),
            "software": d.software,
            "n_assessments": len(d.assessments),
            "n_flagged": len(flagged),
            "n_unflagged": len(unflagged),
            "flagged": flagged[:ASSESSMENT_ROW_LIMIT],
            "unflagged": unflagged[:ASSESSMENT_ROW_LIMIT],
            "flagged_note": _shown(len(flagged), ASSESSMENT_ROW_LIMIT),
            "unflagged_note": _shown(len(unflagged), ASSESSMENT_ROW_LIMIT),
        })
    return views


def _correction_rows(
    corrections: dict[str, CorrectionResult], redactor: Redactor, known_ids: frozenset[str]
) -> list[dict]:
    """One rendered row per candidate, with its free-text reason scrubbed.

    ``CorrectionResult.reason`` is written by the adapter, not by this package: a backend
    that reports a failure by quoting its own log puts that log's absolute paths in the
    reason, so the field gets the same treatment as a detector's free-form payload.
    """

    return [
        {
            "candidate_id": candidate_id,
            "backend": result.candidate.backend,
            "status": result.status.value,
            "n_markers_changed": len(result.markers_changed),
            "n_samples_changed": len(result.samples_changed),
            "output_scale": result.output_scale or "",
            "reason": _redact_free_form(result.reason, redactor, known_ids)
            if result.reason else "",
            "warnings": [_redact_free_form(w, redactor, known_ids) for w in result.warnings],
        }
        for candidate_id, result in sorted(corrections.items())
    ]


def _marker_selection(all_markers: list[str], limit: int) -> tuple[list[str], list[str]]:
    """The markers a figure section draws, and the ones it omits.

    The source list follows the row order of ``markers.tsv``, so an unsorted prefix
    would decide which markers a reader sees by the order somebody happened to write a
    spreadsheet in.  Sorting makes the page's choice a property of the panel rather than
    of the file, and the omitted markers are returned so the page can name them.
    """

    ordered = sorted(all_markers)
    return ordered[:limit], ordered[limit:]


def render_report(
    config_path: str | Path,
    base_dir: str | Path,
    layout,
    diagnosis_only: bool = False,
) -> Path:
    """Render ``report.html``.  ``diagnosis_only`` renders before any correction exists."""

    config = load_config(config_path)
    redactor = Redactor(enabled=config.report.redact_identifiers, base_dir=base_dir)
    dataset = load_dataset(config, base_dir)

    design = DesignAudit.read(layout.design_audit)
    sampling = (
        SamplingManifest.read(layout.sampling_manifest)
        if layout.sampling_manifest.exists() else None
    )

    detectors: list[DetectorResult] = []
    for name in sorted({p.parent.name for p in
                        (layout.root / "diagnostics").glob("*/detector_result.json")}):
        detectors.append(DetectorResult.read(layout.detector_result(name)))

    corrections: dict[str, CorrectionResult] = {}
    metrics: dict[str, MetricResult] = {}
    if not diagnosis_only:
        for candidate in config.correction.candidates:
            cpath = layout.correction_result(candidate.id)
            if cpath.exists():
                corrections[candidate.id] = CorrectionResult.read(cpath)
            mpath = layout.metric_result(candidate.id)
            if mpath.exists():
                metrics[candidate.id] = MetricResult.read(mpath)

    decision = (
        DecisionResult.read(layout.decision)
        if (not diagnosis_only and layout.decision.exists()) else None
    )

    all_records: list[MetricRecord] = []
    for result in metrics.values():
        all_records.extend(result.records)
    if layout.reductions.exists():
        payload = json.loads(layout.reductions.read_text(encoding="utf-8"))
        all_records.extend(MetricRecord.model_validate(r) for r in payload["records"])

    metric_rows = _metric_rows(all_records)
    metric_columns = list(metric_rows[0]) if metric_rows else ["metric", "status"]

    # -- figures ---------------------------------------------------------------------
    batch_of = dataset.sample_batches()
    obs = dataset.obs
    evaluation_markers = list(
        sampling.evaluation_markers if sampling else dataset.marker_ids("evaluation")
    )
    markers, markers_omitted = _marker_selection(evaluation_markers, EXPLORER_MARKER_LIMIT)
    comparison_markers, comparison_omitted = _marker_selection(
        evaluation_markers, COMPARISON_MARKER_LIMIT
    )
    baseline = dataset.expression()
    figures: list[dict] = []
    for marker in markers:
        col = dataset.feature_index([marker])[0]
        by_batch: dict[str, np.ndarray] = {}
        for sid in obs[SAMPLE_ID].unique():
            batch = batch_of.get(str(sid))
            if batch is None:
                continue
            mask = (obs[SAMPLE_ID] == sid).to_numpy()
            by_batch.setdefault(str(batch), []).append(baseline[mask, col])
        merged = {b: np.concatenate(v) for b, v in by_batch.items()}
        figures.append({"title": marker, "svg": plots.density_by_batch(
            merged, f"{marker} by batch (uncorrected)", "expression")})

    comparison_figures: list[dict] = []
    if not diagnosis_only:
        import anndata as ad

        for candidate_id, result in sorted(corrections.items()):
            path = layout.corrected_h5ad(candidate_id)
            if result.status.value != "success" or not path.exists():
                continue
            adata = ad.read_h5ad(path)
            layer = str(adata.uns.get("expression_layer", "corrected"))
            corrected = np.asarray(adata.layers[layer], dtype=np.float64)
            for marker in comparison_markers:
                col = dataset.feature_index([marker])[0]
                before, after = {}, {}
                for batch in sorted({b for b in batch_of.values() if b}):
                    mask = obs[SAMPLE_ID].map(batch_of).to_numpy() == batch
                    before[batch] = baseline[mask, col]
                    after[batch] = corrected[mask, col]
                comparison_figures.append({
                    "candidate": candidate_id, "marker": marker,
                    "svg": plots.before_after(before, after, marker, "expression"),
                })

    composition_figure = None
    if sampling is not None:
        freq = frequency_table(dataset, sampling.labels_by_cell)
        composition_figure = plots.frequency_bars(
            {redactor.identifier(k, "S"): v for k, v in freq.items()},
            {redactor.identifier(k, "S"): (batch_of.get(k) or "") for k in freq},
        )

    # -- assembled context -----------------------------------------------------------
    samples = [{
        "sample_id": redactor.identifier(rec.sample_id, "S"),
        "donor_id": redactor.identifier(rec.donor_id, "D") if rec.donor_id else "",
        "batch": rec.batch or "",
        "group": rec.group or "",
        "stimulation": rec.stimulation or "",
        "timepoint": rec.timepoint or "",
        "reference_id": redactor.identifier(rec.reference_id, "R") if rec.reference_id else "",
        "reference_role": str(rec.reference_role),
        "n_cells": rec.n_cells or 0,
    } for rec in dataset.manifest.samples]

    # Every identifier the run knows, so a detector's or a metric's free text can be
    # scrubbed against the same set the samples table was redacted with.  Collected after
    # ``samples`` so every one of them already has a pseudonym.
    known_ids = frozenset(
        str(value)
        for rec in dataset.manifest.samples
        for value in (rec.sample_id, rec.donor_id, rec.reference_id)
        if value
    )
    metric_rows = _scrub_metric_reasons(metric_rows, redactor, known_ids)

    failures = []
    for result in detectors:
        if result.status != "ok":
            failures.append({
                "kind": "detector", "unit": result.detector, "status": result.status,
                "reason": _redact_free_form(result.reason, redactor, known_ids)
                if result.reason else "",
            })
    for candidate_id, result in sorted(corrections.items()):
        if result.status.value != "success":
            failures.append({
                "kind": "candidate", "unit": candidate_id, "status": result.status.value,
                "reason": _redact_free_form(result.reason, redactor, known_ids)
                if result.reason else "",
            })
    for record in all_records:
        if record.status in (MetricStatus.FAILED, MetricStatus.NOT_IMPLEMENTED):
            failures.append({
                "kind": "metric", "unit": record.metric_name, "status": str(record.status),
                "reason": _redact_free_form(record.reason, redactor, known_ids)
                if record.reason else "",
            })

    joint_rows = [
        r for r in metric_rows if r["metric"] == "joint_structure_distance"
    ]

    context = {
        "title": config.report.title,
        "version": __version__,
        "diagnosis_only": diagnosis_only,
        "config": _redact_config_paths(config.resolved_dict(), redactor),
        "fingerprint": config.analytical_fingerprint(),
        "dataset": dataset.manifest,
        "processing": dataset.processing,
        "samples": samples,
        "design": design,
        "sampling": sampling,
        "detectors": _detector_views(detectors, redactor, known_ids),
        "corrections": corrections,
        "correction_rows": _correction_rows(corrections, redactor, known_ids),
        "decision": decision,
        "metric_rows": metric_rows,
        "metric_columns": metric_columns,
        "metric_csv": _data_uri_csv(metric_rows, metric_columns),
        "figures": figures,
        "explorer_markers": markers,
        "explorer_markers_total": len(evaluation_markers),
        "explorer_markers_omitted": markers_omitted,
        "comparison_markers": comparison_markers,
        "comparison_markers_omitted": comparison_omitted,
        "comparison_figures": comparison_figures,
        "composition_figure": composition_figure,
        "joint_rows": joint_rows[:JOINT_ROW_LIMIT],
        "joint_note": _shown(len(joint_rows), JOINT_ROW_LIMIT),
        "failures": failures,
        "registered_metrics": sorted(
            (f"{ns}.{name}", spec) for (ns, name), spec in registered_metrics().items()
        ),
        "deferred_metrics": DEFERRED_METRICS,
        "upstream_not_found": UPSTREAM_METRICS_NOT_FOUND,
        "redacted": config.report.redact_identifiers,
        "parameter_list_limit": PARAMETER_LIST_LIMIT,
    }

    env = Environment(
        loader=FileSystemLoader(str(templates_root())),
        autoescape=select_autoescape(["html"]),
    )
    html = env.get_template("report.html.j2").render(**context)

    out = layout.report
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    redactor.write_map(out.parent / "identifier_map.json")
    return out
