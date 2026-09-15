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
                reason=f"the CytoBatchFlagR wrapper exited with status {exc.returncode}",
                eligibility=eligibility, software=software, warnings=[exc.log[-2000:]],
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
            for row in r_result.payload.get("iqr", [])
        ]

        return DetectorResult(
            detector=self.name, namespace="native", run_id=context.run_id, status="ok",
            eligibility=eligibility, assessments=assessments,
            tables={k: str(v) for k, v in r_result.payload.get("tables", {}).items()},
            plots={k: str(v) for k, v in r_result.payload.get("plots", {}).items()},
            parameters={"thresholds": dict(zip(params["markers"], params["thresholds"],
                                               strict=True)),
                        "iqr_multiplier": 1.5,
                        "control_interpretation": "measured technical replicate or anchor only"},
            software=software + [
                SoftwareVersion(name="CytoBatchFlagR", source="git checkout",
                                commit=r_result.payload.get("commit"),
                                available=True,
                                detail="sourced unmodified; GPL-3, arm's-length invocation")
            ],
            warnings=list(r_result.payload.get("warnings", [])),
            requirements={"control_in_every_batch": True, "fixed_thresholds": True},
        )


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
