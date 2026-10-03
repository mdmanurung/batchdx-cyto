"""cyCombine detection adapters.

Two upstream entry points are wrapped against the pinned API:

``detect_batch_effect_express``
    quick per-marker EMD across batch pairs, density plots and an MDS plot.

``detect_batch_effect``
    SOM-based clustering plus per-cluster EMD and cluster-percentage outlier tests.
    The audited source restricts that analysis to ``3 <= n_batches <= 30`` and emits
    UMAPs alone outside it, so a two-batch cross-cohort design is declared *ineligible*
    for this detector with a machine-readable reason rather than run and reported as an
    uninformative result.

Both upstream functions communicate through messages and plot files, not return
values.  The structured tables here therefore come from a direct call to the exported
``cyCombine::compute_emd``, never from parsing console output or a rendered image.
"""

from __future__ import annotations

from pathlib import Path

from batchdx.cyto.adapters.r_runner import (
    RNotAvailable,
    RScriptFailed,
    r_package_version,
    r_version,
    run_r_script,
)
from batchdx.cyto.diagnostics.base import DetectorContext, eligibility_for, register_detector
from batchdx.cyto.io.binary import write_matrix
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.schemas.contracts import DetectorAssessment, DetectorResult
from batchdx.cyto.schemas.enums import MetricStatus

R_SCRIPT = "cycombine_detect.R"


class _CyCombineDetector:
    namespace = "native"

    def __init__(self, name: str, mode: str, eligibility_key: str):
        self.name = name
        self.mode = mode
        self.eligibility_key = eligibility_key

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        eligibility = eligibility_for(context.design, self.eligibility_key)
        software = [
            r_version(context.rscript),
            r_package_version("cyCombine", context.rscript),
            r_package_version("emdist", context.rscript),
        ]

        if eligibility is not None and not eligibility.eligible:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="ineligible", reason=eligibility.reason, eligibility=eligibility,
                software=software,
                requirements={"batch_count_rule": eligibility.requirements_met},
            )

        cycombine = next(s for s in software if s.name == "cyCombine")
        if not cycombine.available:
            return DetectorResult(
                detector=self.name, namespace="native", run_id=context.run_id,
                status="not_available",
                reason=(
                    "the cyCombine R package is not installed in this environment "
                    f"({cycombine.detail}). No flags are produced; an unavailable "
                    "dependency is never reported as an absence of batch effects."
                ),
                eligibility=eligibility, software=software,
            )

        out_dir = Path(context.out_dir) / self.name
        out_dir.mkdir(parents=True, exist_ok=True)
        obs = dataset.obs
        markers = dataset.marker_ids("evaluation") or dataset.marker_ids()
        cells = (
            context.sampling.selected_cell_ids
            if context.sampling is not None else obs[CELL_ID].astype(str).tolist()
        )
        matrix_path = out_dir / "input.f64"
        write_matrix(
            matrix_path,
            dataset.expression(features=markers, cell_ids=cells),
            row_ids=cells, col_names=markers,
        )
        batch_of = dataset.sample_batches()
        sample_of = dict(zip(obs[CELL_ID].astype(str), obs[SAMPLE_ID].astype(str), strict=True))
        batches = [str(batch_of.get(sample_of[c])) for c in cells]
        samples = [sample_of[c] for c in cells]
        labels = (
            [context.sampling.labels_by_cell[c] for c in cells]
            if context.sampling is not None else None
        )

        params = {
            "mode": self.mode,
            "input_matrix": str(matrix_path),
            "out_dir": str(out_dir),
            "markers": markers,
            "batch": batches,
            # cyCombine's detector groups by a column literally named "sample"
            # (detect_batch_effect.R:200, hard-coded with no rename parameter), so
            # the per-cell sample identity has to cross the boundary alongside the
            # per-cell batch labels.  Passing it is a requirement of the pinned API.
            "sample": samples,
            "label": labels,
            "seed": context.seed,
            "bin_size": 0.1,
        }
        try:
            r_result = run_r_script(R_SCRIPT, params, out_dir, rscript=context.rscript,
                                    log_path=out_dir / "cycombine_detect.log")
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
                reason=f"the R detector exited with status {exc.returncode}",
                eligibility=eligibility, software=software,
                warnings=[exc.log[-2000:]],
            )

        assessments: list[DetectorAssessment] = []
        for row in r_result.payload.get("emd", []):
            assessments.append(DetectorAssessment(
                assessment_id=(
                    f"emd:{row.get('population', 'all')}:{row['marker']}:"
                    f"{row['batch_a']}-{row['batch_b']}"
                ),
                detector=self.name, namespace="native", scope="batch_pair_marker",
                marker=str(row["marker"]), population=row.get("population"),
                statistic="cycombine_emd", value=float(row["value"]),
                status=MetricStatus.OK, evidence_class="descriptive",
                reason=(
                    "binned Earth Mover's Distance as computed by cyCombine::compute_emd "
                    "(bin units, hist density weights)"
                ),
            ))
        for row in r_result.payload.get("cluster_percentage", []):
            assessments.append(DetectorAssessment(
                assessment_id=f"clusterpct:{row['population']}:{row['batch']}",
                detector=self.name, namespace="native", scope="batch_population",
                population=str(row["population"]), batch=str(row["batch"]),
                statistic="cluster_percentage", value=float(row["value"]),
                status=MetricStatus.OK, evidence_class="descriptive",
                reason="percentage of the batch's cells in this SOM node",
            ))

        return DetectorResult(
            detector=self.name, namespace="native", run_id=context.run_id, status="ok",
            eligibility=eligibility, assessments=assessments,
            plots={k: str(v) for k, v in r_result.payload.get("plots", {}).items()},
            tables={k: str(v) for k, v in r_result.payload.get("tables", {}).items()},
            parameters=params, software=software,
            warnings=list(r_result.payload.get("warnings", [])),
            requirements={
                "note": (
                    "upstream detect_batch_effect* return no structured value; tables here "
                    "are produced by calling cyCombine::compute_emd directly"
                )
            },
        )


register_detector(_CyCombineDetector("cycombine_express", "express", "cycombine_express"))
register_detector(_CyCombineDetector("cycombine_cluster", "cluster", "cycombine_cluster"))
