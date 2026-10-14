"""Render the unified offline HTML report from validated structured artifacts.

The renderer reads only what is on disk.  It never recomputes a score, so every number
in the report is traceable to a record, and regenerating the report cannot trigger a
correction.
"""

from __future__ import annotations

import base64
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
from batchdx.cyto.reporting.redact import Redactor, run_salt
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

MAX_EXPLORER_MARKERS = 8


def _data_uri_csv(rows: list[dict], columns: list[str]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    encoded = base64.b64encode(buffer.getvalue().encode()).decode()
    return f"data:text/csv;base64,{encoded}"


def _metric_rows(records: list[MetricRecord], redactor: Redactor) -> list[dict]:
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


def render_report(
    config_path: str | Path,
    base_dir: str | Path,
    layout,
    diagnosis_only: bool = False,
) -> Path:
    """Render ``report.html``.  ``diagnosis_only`` renders before any correction exists."""

    config = load_config(config_path)
    redactor = Redactor(enabled=config.report.redact_identifiers, base_dir=base_dir,
                        salt=run_salt(layout.root))
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

    # The exports.* configuration is consumed here: when FCS export is requested,
    # every successful candidate's corrected matrix is written as one FCS per
    # sample under the run's exports directory. Nothing is exported by default.
    if config.exports.fcs and not diagnosis_only:
        from batchdx.cyto.io.dataset import write_candidate_output  # noqa: F401
        from batchdx.cyto.io.export import export_fcs_per_sample

        for candidate in config.correction.candidates:
            result_path = layout.correction_result(candidate.id)
            corrected = layout.corrected_h5ad(candidate.id)
            if not (result_path.exists() and corrected.exists()):
                continue
            result = CorrectionResult.read(result_path)
            if result.status is not CorrectionStatus.SUCCESS:
                continue
            import anndata as ad

            adata = ad.read_h5ad(corrected)
            layer_name = str(adata.uns.get("expression_layer", "corrected"))
            matrix = np.asarray(adata.layers[layer_name], dtype=np.float64)
            out_dir = layout.root / "exports" / candidate.id / "fcs"
            export_fcs_per_sample(
                dataset, matrix, out_dir, candidate_id=candidate.id, invert=True,
            )

    all_records: list[MetricRecord] = []
    for result in metrics.values():
        all_records.extend(result.records)
    if layout.reductions.exists():
        payload = json.loads(layout.reductions.read_text(encoding="utf-8"))
        all_records.extend(MetricRecord.model_validate(r) for r in payload["records"])

    metric_rows = _metric_rows(all_records, redactor)
    metric_columns = list(metric_rows[0]) if metric_rows else ["metric", "status"]

    # -- figures ---------------------------------------------------------------------
    batch_of = dataset.sample_batches()
    obs = dataset.obs
    markers = (sampling.evaluation_markers if sampling else dataset.marker_ids("evaluation"))
    markers = markers[:MAX_EXPLORER_MARKERS]
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
            for marker in markers[:3]:
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

    failures = []
    for result in detectors:
        if result.status != "ok":
            failures.append({"kind": "detector", "unit": result.detector,
                             "status": result.status, "reason": result.reason or ""})
    for candidate_id, result in sorted(corrections.items()):
        if result.status.value != "success":
            failures.append({"kind": "candidate", "unit": candidate_id,
                             "status": result.status.value, "reason": result.reason or ""})
    for record in all_records:
        if record.status in (MetricStatus.FAILED, MetricStatus.NOT_IMPLEMENTED):
            failures.append({"kind": "metric", "unit": record.metric_name,
                             "status": str(record.status), "reason": record.reason or ""})

    context = {
        "title": config.report.title,
        "version": __version__,
        "diagnosis_only": diagnosis_only,
        "config": redactor.strip_paths(config.resolved_dict()),
        "fingerprint": config.analytical_fingerprint(),
        "dataset": dataset.manifest,
        "processing": dataset.processing,
        "samples": samples,
        "design": design,
        "sampling": sampling,
        "detectors": detectors,
        "corrections": corrections,
        "decision": decision,
        "metric_rows": metric_rows,
        "metric_columns": metric_columns,
        "metric_csv": _data_uri_csv(metric_rows, metric_columns),
        "figures": figures,
        "comparison_figures": comparison_figures,
        "composition_figure": composition_figure,
        "failures": failures,
        "registered_metrics": sorted(
            (f"{ns}.{name}", spec) for (ns, name), spec in registered_metrics().items()
        ),
        "deferred_metrics": DEFERRED_METRICS,
        "upstream_not_found": UPSTREAM_METRICS_NOT_FOUND,
        "redacted": config.report.redact_identifiers,
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
