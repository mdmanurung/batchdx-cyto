"""Acquisition-drift diagnostic.

Within-file drift is a distinct failure mode: a signal that slides across the
acquisition axis inside a single file is not a batch offset and is not fixed by
batch-level correction.  Reporting it separately stops a drift problem being mistaken
for a correctable batch effect.

Per sample and marker the diagnostic reports a rank correlation against the
acquisition axis (monotonic trend) and the largest step between consecutive contiguous
blocks relative to the sample's own robust spread (discontinuity).  A per-sample
event-rate diagnostic reports events per unit time along the same axis and flags rate
discontinuities, because an acquisition interruption is a distinct cause from a
gradual signal slide.

The acquisition axis is the FCS Time channel when it was ingested (the obs column
``acquisition_time``); otherwise the event index is used, and every assessment says
which axis it used.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from batchdx.cyto.diagnostics.base import DetectorContext, register_detector
from batchdx.cyto.io.dataset import (
    ACQUISITION_ORDER,
    ACQUISITION_TIME,
    EVENT_INDEX,
    SAMPLE_ID,
    Dataset,
)
from batchdx.cyto.schemas.contracts import DetectorAssessment, DetectorResult
from batchdx.cyto.schemas.enums import MetricStatus

N_BLOCKS = 10
MIN_CELLS = 100

#: Substrings whose presence in the processing manifest's qc_history /
#: transformation / notes counts as a declared bead-normalization or CS&T-style step.
BEAD_KEYWORDS = ("bead", "cs&t", "cst", "calibrat")


def _order_column(dataset: Dataset) -> tuple[str | None, str]:
    """Pick the acquisition axis and state explicitly which axis is used."""
    obs = dataset.adata.obs
    if ACQUISITION_TIME in obs.columns:
        return ACQUISITION_TIME, (
            "the FCS Time channel, ingested as the obs column 'acquisition_time'"
        )
    if ACQUISITION_ORDER in obs.columns:
        return ACQUISITION_ORDER, (
            "the event index stored as obs['acquisition_order']; no FCS Time channel was "
            "ingested, so this axis is event order, not clock time"
        )
    if EVENT_INDEX in obs.columns:
        return EVENT_INDEX, (
            "the original event index; no FCS Time channel was ingested, so this axis is "
            "event order, not clock time"
        )
    return None, "no acquisition time, acquisition order or event index is available"


def _bead_normalization_declared(processing) -> tuple[bool, list[str]]:
    """Structured check of the processing manifest for a declared bead/CS&T step."""
    declared = [str(q) for q in (processing.qc_history or [])]
    if processing.transformation:
        declared.append(str(processing.transformation))
    if processing.notes:
        declared.append(str(processing.notes))
    hits = [d for d in declared if any(k in d.lower() for k in BEAD_KEYWORDS)]
    return bool(hits), hits


class AcquisitionDriftDetector:
    """Monotonic trend, discontinuity and event rate against the acquisition axis."""

    name = "acquisition_drift"
    namespace = "common"
    eligibility_key = "acquisition_drift"

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        order_col, order_note = _order_column(dataset)
        software: list = []
        if order_col is None:
            return DetectorResult(
                detector=self.name, namespace="common", run_id=context.run_id,
                status="ineligible", reason=order_note,
                requirements={"acquisition_order_column": False},
            )

        bead_declared, bead_evidence = _bead_normalization_declared(dataset.processing)
        warnings: list[str] = []
        if not bead_declared:
            warnings.append(
                "no bead-normalization or CS&T-style step is declared in the processing "
                "manifest (qc_history/transformation). Within-file drift is normally "
                "remedied by bead normalization or recalibration before any batch "
                "correction, and cyCombine assumes bead/dead-cell/doublet pre-gating; "
                "without a declared calibration step, a within-file trend may reflect "
                "remediable instrument drift rather than a batch effect."
            )

        obs = dataset.obs
        markers = dataset.marker_ids("evaluation") or dataset.marker_ids()
        values = dataset.expression(features=markers)
        assessments: list[DetectorAssessment] = []

        for sid in sorted(obs[SAMPLE_ID].unique()):
            mask = (obs[SAMPLE_ID] == sid).to_numpy()
            n = int(mask.sum())
            order = obs.loc[mask, order_col].to_numpy(dtype=float)
            block = values[mask, :]
            sort_idx = np.argsort(order, kind="stable")
            block = block[sort_idx, :]

            # -- per-sample event-rate diagnostic -----------------------------------
            rate_id = f"rate:{sid}"
            if n < MIN_CELLS:
                assessments.append(DetectorAssessment(
                    assessment_id=rate_id, detector=self.name, namespace="common",
                    scope="sample", sample_id=str(sid),
                    statistic="event_rate_max_block_step_over_median",
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=f"{n} events, fewer than the {MIN_CELLS} required; {order_note}",
                ))
            elif not np.isfinite(order).all():
                assessments.append(DetectorAssessment(
                    assessment_id=rate_id, detector=self.name, namespace="common",
                    scope="sample", sample_id=str(sid),
                    statistic="event_rate_max_block_step_over_median",
                    status=MetricStatus.UNDEFINED,
                    reason=f"the acquisition axis has missing values in this sample; "
                           f"{order_note}",
                ))
            else:
                t = np.sort(order)
                rates: list[float] = []
                zero_span = False
                for b in np.array_split(np.arange(n), N_BLOCKS):
                    span = float(t[b[-1]] - t[b[0]])
                    if span <= 0:
                        zero_span = True
                        break
                    rates.append(len(b) / span)
                med = float(np.median(rates)) if rates else 0.0
                if zero_span or not rates or med <= 0:
                    assessments.append(DetectorAssessment(
                        assessment_id=rate_id, detector=self.name, namespace="common",
                        scope="sample", sample_id=str(sid),
                        statistic="event_rate_max_block_step_over_median",
                        status=MetricStatus.UNDEFINED,
                        reason="a contiguous block of the acquisition axis has zero span "
                               "(or the axis is constant), so the event rate is undefined",
                    ))
                else:
                    step = float(np.max(np.abs(np.diff(rates))) / med)
                    assessments.append(DetectorAssessment(
                        assessment_id=rate_id, detector=self.name, namespace="common",
                        scope="sample", sample_id=str(sid),
                        statistic="event_rate_max_block_step_over_median",
                        value=step, flag=step > context.drift_rate_threshold,
                        threshold=context.drift_rate_threshold,
                        status=MetricStatus.OK, evidence_class="descriptive",
                        reason=f"largest step between the event rates of {N_BLOCKS} "
                               f"contiguous time blocks, relative to the median rate; "
                               f"{order_note}",
                    ))

            if n < MIN_CELLS:
                for marker in markers:
                    assessments.append(DetectorAssessment(
                        assessment_id=f"drift:{sid}:{marker}", detector=self.name,
                        namespace="common", scope="sample_marker", sample_id=str(sid),
                        marker=marker, statistic="spearman_rho_vs_acquisition_order",
                        status=MetricStatus.INSUFFICIENT_DATA,
                        reason=f"{n} events, fewer than the {MIN_CELLS} required",
                    ))
                continue

            for j, marker in enumerate(markers):
                column = block[:, j]
                rho, pvalue = stats.spearmanr(np.arange(n), column)
                if np.isfinite(rho):
                    assessments.append(DetectorAssessment(
                        assessment_id=f"drift:{sid}:{marker}", detector=self.name,
                        namespace="common", scope="sample_marker", sample_id=str(sid),
                        marker=marker, statistic="spearman_rho_vs_acquisition_order",
                        value=float(rho),
                        flag=abs(float(rho)) > context.drift_rho_threshold,
                        threshold=context.drift_rho_threshold,
                        status=MetricStatus.OK, evidence_class="descriptive",
                        reason=f"p={float(pvalue):.3g}; {order_note}",
                    ))
                else:
                    # A constant (or fully tied) marker has no defined rank trend.  An
                    # undefined statistic is never recorded as a measured zero.
                    assessments.append(DetectorAssessment(
                        assessment_id=f"drift:{sid}:{marker}", detector=self.name,
                        namespace="common", scope="sample_marker", sample_id=str(sid),
                        marker=marker, statistic="spearman_rho_vs_acquisition_order",
                        value=None, status=MetricStatus.UNDEFINED,
                        reason=f"the marker is constant in this sample, so the rank trend "
                               f"is undefined; {order_note}",
                    ))

                blocks = np.array_split(column, N_BLOCKS)
                medians = np.array([np.median(b) for b in blocks if len(b)])
                spread = float(stats.median_abs_deviation(column, scale="normal"))
                if spread <= 0 or len(medians) < 2:
                    assessments.append(DetectorAssessment(
                        assessment_id=f"jump:{sid}:{marker}", detector=self.name,
                        namespace="common", scope="sample_marker", sample_id=str(sid),
                        marker=marker, statistic="max_block_step_over_mad",
                        status=MetricStatus.UNDEFINED,
                        reason="the marker has zero robust spread in this sample",
                    ))
                    continue
                step = float(np.max(np.abs(np.diff(medians))) / spread)
                assessments.append(DetectorAssessment(
                    assessment_id=f"jump:{sid}:{marker}", detector=self.name,
                    namespace="common", scope="sample_marker", sample_id=str(sid),
                    marker=marker, statistic="max_block_step_over_mad",
                    value=step, status=MetricStatus.OK, evidence_class="descriptive",
                    reason=f"largest step between {N_BLOCKS} contiguous blocks, in robust SD units",
                ))

        flagged = sorted({
            f"{a.sample_id}/{a.marker}" for a in assessments
            if a.flag and a.statistic == "spearman_rho_vs_acquisition_order"
        })
        if flagged:
            warnings.append(
                f"{len(flagged)} sample-marker combination(s) show a monotonic within-file "
                "trend. Within-batch drift is not addressed by batch-level correction; "
                "the standard remediations are bead normalization or recalibration at "
                "acquisition time, and exclusion of the affected samples."
            )

        return DetectorResult(
            detector=self.name, namespace="common", run_id=context.run_id, status="ok",
            assessments=assessments,
            requirements={
                "acquisition_order_column": order_col,
                "acquisition_axis": order_note,
                "min_cells": MIN_CELLS,
                "bead_normalization_declared": bead_declared,
                "bead_normalization_evidence": bead_evidence,
            },
            parameters={
                "n_blocks": N_BLOCKS,
                "rho_threshold": context.drift_rho_threshold,
                "rate_step_threshold": context.drift_rate_threshold,
            },
            warnings=warnings, software=software,
        )


register_detector(AcquisitionDriftDetector())
