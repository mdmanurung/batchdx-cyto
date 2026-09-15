"""Acquisition-drift diagnostic.

Within-file drift is a distinct failure mode: a signal that slides across acquisition
order inside a single file is not a batch offset and is not fixed by batch-level
correction.  Reporting it separately stops a drift problem being mistaken for a
correctable batch effect.

Per sample and marker the diagnostic reports a rank correlation against acquisition
order (monotonic trend) and the largest step between consecutive contiguous blocks
relative to the sample's own robust spread (discontinuity).
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from batchdx.cyto.diagnostics.base import DetectorContext, register_detector
from batchdx.cyto.io.dataset import ACQUISITION_ORDER, EVENT_INDEX, SAMPLE_ID, Dataset
from batchdx.cyto.schemas.contracts import DetectorAssessment, DetectorResult
from batchdx.cyto.schemas.enums import MetricStatus

N_BLOCKS = 10
MIN_CELLS = 100


def _order_column(dataset: Dataset) -> tuple[str | None, str]:
    obs = dataset.adata.obs
    if ACQUISITION_ORDER in obs.columns:
        return ACQUISITION_ORDER, "acquisition order / time channel recorded at ingestion"
    if EVENT_INDEX in obs.columns:
        return EVENT_INDEX, (
            "original event index used as a proxy for acquisition order; no time channel "
            "was recorded, so a genuine time axis may differ"
        )
    return None, "no acquisition order or event index is available"


class AcquisitionDriftDetector:
    """Monotonic trend and discontinuity against acquisition order, per sample and marker."""

    name = "acquisition_drift"
    namespace = "common"
    eligibility_key = "acquisition_drift"

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        order_col, order_note = _order_column(dataset)
        software = []
        if order_col is None:
            return DetectorResult(
                detector=self.name, namespace="common", run_id=context.run_id,
                status="ineligible", reason=order_note,
                requirements={"acquisition_order_column": False},
            )

        obs = dataset.obs
        markers = dataset.marker_ids("evaluation") or dataset.marker_ids()
        values = dataset.expression(features=markers)
        assessments: list[DetectorAssessment] = []
        warnings: list[str] = []

        for sid in sorted(obs[SAMPLE_ID].unique()):
            mask = (obs[SAMPLE_ID] == sid).to_numpy()
            n = int(mask.sum())
            order = obs.loc[mask, order_col].to_numpy(dtype=float)
            block = values[mask, :]
            sort_idx = np.argsort(order, kind="stable")
            block = block[sort_idx, :]

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
                rho = float(rho) if np.isfinite(rho) else 0.0
                assessments.append(DetectorAssessment(
                    assessment_id=f"drift:{sid}:{marker}", detector=self.name,
                    namespace="common", scope="sample_marker", sample_id=str(sid),
                    marker=marker, statistic="spearman_rho_vs_acquisition_order",
                    value=rho, flag=abs(rho) > context.drift_rho_threshold,
                    threshold=context.drift_rho_threshold,
                    status=MetricStatus.OK, evidence_class="descriptive",
                    reason=f"p={float(pvalue):.3g}; {order_note}",
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
                "trend. Within-batch drift is not addressed by batch-level correction."
            )

        return DetectorResult(
            detector=self.name, namespace="common", run_id=context.run_id, status="ok",
            assessments=assessments,
            requirements={"acquisition_order_column": order_col, "min_cells": MIN_CELLS},
            parameters={"n_blocks": N_BLOCKS, "rho_threshold": context.drift_rho_threshold},
            warnings=warnings, software=software,
        )


register_detector(AcquisitionDriftDetector())
