"""CytoBatchFlagR detection adapter (published as CytoBatchFlagR; formerly cytoFlagR).

Upstream is distributed as a set of R scripts under GPL-3, not as an installable
package, and its own entry point is interactive (``rstudioapi::selectDirectory`` and
Shiny helpers).  The licence decision recorded in
``docs/decisions/0003-cytobatchflagr-integration.md`` is therefore to integrate at
arm's length: the upstream scripts are ``source()``d *unmodified* from a checkout the
user supplies, in a separate R process, and the analysis functions are called directly
with values taken from configuration.  Nothing upstream is copied into this package
and nothing is refactored, so batchdx-cyto does not become a derivative work.

The three IQR metrics (negative median intensity, positive median intensity, percent
positive) are computed per control sample per marker across batches, and a batch is
flagged when a value falls beyond 1.5 IQR of the quartiles.  That makes a control
sample in every batch a hard requirement, and it makes the rule degenerate below four
batches; both constraints are encoded in the Stage 3 eligibility manifest.

Two failure modes are handled here rather than by the upstream code, both of which
would otherwise be reported as an absence of batch effects:

* ``neg_mfi_df_all`` and ``pos_mfi_df_all`` default to ``min_cells = TRUE`` and silently
  drop any control with fewer than 100 cells on that side of a marker threshold;
  ``percent_pos_df_all`` has no such parameter and drops a control the same way by
  emitting no row for it.  The wrapper passes ``min_cells = FALSE`` where the argument
  exists and this module instead requires every batch to contribute a value built on at
  least :data:`MIN_CELLS_PER_CONTROL_MARKER` cells, so a dropped control becomes a stated
  ``ineligible`` rather than a missing row that quietly drops the per-marker value count
  below the four the Tukey fence needs.
* The upstream ``library()`` calls are unconditional, so one absent package aborts the
  whole ``source()``.  The wrapper preflights the full set and reports the absence, and
  this module maps that to ``not_available`` with the remedy, rather than letting a raw
  R error surface as ``status="failed"``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from batchdx.cyto.adapters.r_runner import (
    RNotAvailable,
    RScriptFailed,
    r_version,
    run_r_script,
)
from batchdx.cyto.diagnostics.base import DetectorContext, eligibility_for, register_detector
from batchdx.cyto.io.binary import write_matrix
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.schemas.contracts import DetectorAssessment, DetectorResult, SoftwareVersion
from batchdx.cyto.schemas.enums import MetricStatus, ReferenceRole

R_SCRIPT = "cytobatchflagr_detect.R"
REQUIRED_UPSTREAM_FILE = "R/CytoBatchFlagR_functions.R"

#: The three IQR metrics, in the order ``iqr_dataframe`` emits them.  Upstream relabels
#: the ``type`` column to these strings in the CSV it writes but not in the frame it
#: returns, so this is the mapping the coverage check has to undo.
_IQR_TYPES = ("-MFI", "+MFI", "%pos")

#: Upstream's own minimum cell count, adopted here rather than invented.  With
#: ``min_cells = TRUE`` (the default) ``neg_mfi_df_all`` and ``pos_mfi_df_all`` drop any
#: control holding fewer than this many cells on one side of a marker threshold
#: (functions.R:1084 ``small <- names(freq[freq < 100])``, applied at :1090-1102).  The
#: wrapper passes ``min_cells = FALSE`` and this module re-applies the same number as an
#: explicit, reported precondition: the median of a 20-cell population is not the
#: per-control signal the IQR rule compares across batches, and dropping the row silently
#: would leave the fence comparing three values instead of four and therefore unable to
#: fire on a batch that really did shift.
MIN_CELLS_PER_CONTROL_MARKER = 100


class CytoBatchFlagRDetector:
    """Non-interactive wrapper around the upstream IQR and EMD assessments."""

    name = "cytobatchflagr"
    namespace = "native"
    eligibility_key = "cytobatchflagr"

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        eligibility = eligibility_for(context.design, self.eligibility_key)
        software: list[SoftwareVersion] = [r_version(context.rscript)]

        if eligibility is not None and not eligibility.eligible:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="ineligible", reason=eligibility.reason, eligibility=eligibility,
                software=software,
                requirements={"control_in_every_batch": True, "min_batches": 4},
            )

        checkout = Path(context.upstream_dir) / "CytoBatchFlagR" if context.upstream_dir else None
        if checkout is None or not (checkout / REQUIRED_UPSTREAM_FILE).exists():
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="not_available",
                reason=(
                    "no CytoBatchFlagR checkout was found. Upstream ships R scripts rather "
                    "than an installable package, so this adapter needs a checkout at the "
                    "commit pinned in config/upstream_lock.yaml; run scripts/fetch_upstream.sh "
                    "and set upstream.checkout_dir. Absent controls or absent code yield an "
                    "explicit limitation, never a fabricated flag."
                ),
                eligibility=eligibility, software=software,
                requirements={"upstream_file": REQUIRED_UPSTREAM_FILE},
            )

        controls = {
            sid: rec for sid, rec in dataset.samples.items()
            if rec.reference_role in (ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR)
        }
        if not controls:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="ineligible",
                reason=(
                    "no sample declares a measured reference role; the three IQR metrics are "
                    "defined per control sample per marker across batches"
                ),
                eligibility=eligibility, software=software,
            )

        out_dir = Path(context.out_dir) / self.name
        out_dir.mkdir(parents=True, exist_ok=True)

        thresholds = _threshold_table(dataset, context)
        if thresholds is None:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="ineligible",
                reason=(
                    "no fixed positivity thresholds are available. The published evaluation "
                    "reports five experienced analysts choosing five different thresholds for "
                    "the same markers, so this adapter refuses to pick them itself: supply "
                    "evaluation.gating_thresholds in the declared coordinate system."
                ),
                eligibility=eligibility, software=software,
            )

        obs = dataset.obs
        control_cells = obs.loc[obs[SAMPLE_ID].isin(controls.keys())]
        markers = sorted(thresholds["feature_id"].unique())
        cells = control_cells[CELL_ID].astype(str).tolist()
        matrix_path = out_dir / "controls.f64"
        write_matrix(matrix_path, dataset.expression(features=markers, cell_ids=cells),
                     row_ids=cells, col_names=markers)

        batch_of = dataset.sample_batches()
        params = {
            "upstream_dir": str(checkout),
            "input_matrix": str(matrix_path),
            "out_dir": str(out_dir),
            "markers": markers,
            "thresholds": [float(thresholds.loc[thresholds["feature_id"] == m, "threshold"].iloc[0])
                           for m in markers],
            "control": control_cells[SAMPLE_ID].map(
                {sid: (rec.reference_id or sid) for sid, rec in controls.items()}
            ).astype(str).tolist(),
            "batch": [str(batch_of.get(str(s))) for s in control_cells[SAMPLE_ID]],
        }

        try:
            r_result = run_r_script(R_SCRIPT, params, out_dir, rscript=context.rscript,
                                    log_path=out_dir / "cytobatchflagr.log")
        except RNotAvailable as exc:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="not_available", reason=str(exc), eligibility=eligibility,
                software=software,
            )
        except RScriptFailed as exc:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="failed",
                reason=f"the CytoBatchFlagR wrapper failed: {exc}",
                eligibility=eligibility, software=software, warnings=[exc.log[-2000:]],
            )

        payload = r_result.payload
        if payload.get("status") == "not_available":
            # The wrapper's own preflight found the environment unable to satisfy
            # upstream's unconditional library() block.  That is unavailability, not a
            # failed detection, and it is the difference between "no flag" and "no
            # answer": R was never reached.
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="not_available", reason=str(payload.get("reason") or ""),
                eligibility=eligibility, software=software,
                requirements={"missing_packages": list(payload.get("missing_packages") or [])},
            )

        coverage = _incomplete_panels(payload)
        if coverage:
            # Every batch must contribute a value built on at least
            # MIN_CELLS_PER_CONTROL_MARKER cells, or the fence is being applied to fewer
            # than four values and cannot fire.  A short panel is an unanswerable
            # question, not a negative one.
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="ineligible",
                reason=(
                    "the IQR metrics are not estimable for every batch: "
                    + "; ".join(
                        f"{n_short} of {n_batches} batches hold fewer than "
                        f"{MIN_CELLS_PER_CONTROL_MARKER} cells on the side the {mtype} metric "
                        f"is computed from, for marker {marker!r} in control {control!r} "
                        f"({', '.join(batches)})"
                        for mtype, marker, control, batches, n_batches, n_short in coverage
                    )
                    + ". The fence would then be evaluated on fewer than four values and could "
                    "not fire, which would present a batch effect as an absence of one. Upstream's "
                    "own helpers drop such a control silently (min_cells; functions.R:1084), so "
                    "the wrapper disables that and the shortfall is reported here instead. Raise "
                    "cells per control, or narrow the marker set to the markers a control of this "
                    "size can support on both sides of its threshold."
                ),
                eligibility=eligibility, software=software,
                tables={k: str(v) for k, v in _mapping(payload, "tables").items()},
                requirements={"control_in_every_batch": True, "min_batches": 4,
                              "min_cells_per_control_marker":
                                  MIN_CELLS_PER_CONTROL_MARKER},
            )

        assessments = [
            DetectorAssessment(
                assessment_id=f"iqr:{row['type']}:{row['control']}:{row['marker']}:{row['batch']}",
                detector=self.name, namespace="native", scope="control_marker_batch",
                marker=str(row["marker"]), batch=str(row["batch"]),
                sample_id=str(row["control"]), statistic=str(row["type"]),
                value=None if row.get("value") is None else float(row["value"]),
                flag=bool(row.get("batch_effect", 0)),
                threshold=1.5, status=MetricStatus.OK,
                evidence_class="technical_reference",
                reason=(
                    "flagged when the value falls beyond 1.5 IQR of the quartiles taken "
                    "across batches for this control and marker"
                ),
            )
            for row in payload.get("iqr", [])
        ]

        return DetectorResult(
            detector=self.name, namespace="native", run_id=context.run_id, status="ok",
            eligibility=eligibility, assessments=assessments,
            tables={k: str(v) for k, v in _mapping(payload, "tables").items()},
            plots={k: str(v) for k, v in _mapping(payload, "plots").items()},
            parameters={"thresholds": dict(zip(params["markers"], params["thresholds"],
                                               strict=True)),
                        "iqr_multiplier": 1.5,
                        "control_interpretation": "measured technical replicate or anchor only",
                        "min_cells": False,
                        "min_cells_per_control_marker": MIN_CELLS_PER_CONTROL_MARKER},
            software=software + [
                SoftwareVersion(name="CytoBatchFlagR", source="git checkout",
                                commit=payload.get("commit"),
                                available=True,
                                detail="sourced unmodified; GPL-3, arm's-length invocation")
            ],
            warnings=list(payload.get("warnings", [])),
            requirements={"control_in_every_batch": True, "fixed_thresholds": True,
                          "min_cells_per_control_marker": MIN_CELLS_PER_CONTROL_MARKER},
        )


def _mapping(payload: dict, key: str) -> dict:
    """A payload field that is a mapping by contract, read defensively.

    jsonlite serialises an empty R ``list()`` as ``[]`` rather than ``{}``, and the
    wrapper is a separate program whose serialiser this package does not control.  An
    absent map therefore arrives as either, and both mean the same thing: nothing.
    """

    return dict(payload.get(key) or {})


def _incomplete_panels(payload: dict) -> list[tuple[str, str, str, list[str], int, int]]:
    """Control x batch x metric x marker panels a batch is missing from, with the reason.

    ``iqr_metric_function`` groups by control and applies the 1.5 IQR fence to that
    control's per-batch values for one metric and marker.  A batch that contributes no
    value is a batch the fence is not applied to, and below four values it cannot fire at
    all, so a shortfall has to be reported rather than absorbed.  Each entry is
    ``(metric_type, marker, control, batches, n_batches, n_short)`` where ``batches`` are
    the ones whose cell count is under :data:`MIN_CELLS_PER_CONTROL_MARKER` on the side
    that metric needs.
    """

    coverage = payload.get("coverage") or []
    if not coverage:
        return []
    frame = pd.DataFrame(coverage)
    required = {"control", "batch", "marker", "n_negative", "n_positive"}
    if not required <= set(frame.columns):
        return []

    # Which side of the threshold each metric is computed from.  "-MFI" is the median of
    # the cells below the threshold and "+MFI" the median of those above; "%pos" is a
    # proportion over all cells, so it is only undefined when the control is empty.
    sides = {"-MFI": "n_negative", "+MFI": "n_positive", "%pos": "n_cells"}

    out: list[tuple[str, str, str, list[str], int, int]] = []
    for (control, marker), group in frame.groupby(["control", "marker"], sort=True):
        batches = sorted({str(b) for b in group["batch"]})
        for mtype in _IQR_TYPES:
            side = sides[mtype]
            short = sorted(
                str(r["batch"]) for _, r in group.iterrows()
                if int(r[side]) < MIN_CELLS_PER_CONTROL_MARKER
            )
            if short:
                out.append((mtype, str(marker), str(control), short, len(batches), len(short)))
    return out


def _threshold_table(dataset: Dataset, context: DetectorContext) -> pd.DataFrame | None:
    """Fixed thresholds in the declared coordinate system, or ``None`` if unavailable."""

    path = getattr(context, "threshold_path", None)
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep)
    if not {"feature_id", "threshold"} <= set(df.columns):
        return None
    known = set(dataset.adata.var_names)
    return df.loc[df["feature_id"].isin(known), ["feature_id", "threshold"]].drop_duplicates(
        subset=["feature_id"]
    )


register_detector(CytoBatchFlagRDetector())
