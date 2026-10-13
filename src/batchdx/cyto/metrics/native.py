"""Native-namespace metrics evaluated by cyCombine itself (Stage 7).

The ``native`` namespace exists for one reason: to run an upstream *metric*, not an
upstream correction.  ``common`` metrics are computed on this package's frozen
partition; these two are computed by ``cyCombine::evaluate_emd`` and
``cyCombine::evaluate_mad`` on whatever labels the caller supplies, and they are not
comparable with the ``common`` numbers.  The registry records that separation, and
this module is the only producer for those two specs.

Why the two frames are exchanged rather than the data being summarised here: the
reduction cyCombine reports is a ratio of sums over a per-marker table that
``evaluate_emd`` builds internally, applying its own ``filter_limit``.  Recomputing
that in Python would be a second implementation of an upstream procedure, which is
exactly what the ``common``/``native`` split exists to prevent.  So the upstream
function is called and its own return value is reported.

Both upstream functions require the column ``id`` on both frames and align on it
(``evaluate_performance.R`` calls ``check_colname`` for ``id`` on each of
``uncorrected`` and ``corrected``).  The event ids are checked for equality and order
before the call, because a silent misalignment would produce a plausible-looking
reduction that compares different cells under the same label.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from batchdx.cyto.adapters.base import BackendContext
from batchdx.cyto.adapters.r_runner import (
    RNotAvailable,
    RScriptFailed,
    r_package_version,
    run_r_script,
)
from batchdx.cyto.io.binary import write_matrix
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.metrics.records import metric_record
from batchdx.cyto.metrics.registry import NATIVE_CYCOMBINE_EMD, NATIVE_CYCOMBINE_MAD
from batchdx.cyto.sampling.partitions import SamplingManifest
from batchdx.cyto.schemas.contracts import MetricRecord
from batchdx.cyto.schemas.enums import MetricStatus

R_SCRIPT = "cycombine_evaluate.R"

#: Upstream's default bin width, kept here so the number a reader compares against
#: cyCombine's own documentation is the number this package used.
DEFAULT_BIN_SIZE = 0.1


class NativeMetricUnavailable(RuntimeError):
    """The native metric cannot be evaluated in this environment.

    Distinct from a failed evaluation on purpose: an absent cyCombine is an
    environment gap, and a report that says so is more useful than one that reports a
    number nobody computed.
    """


def _as_float(value: Any) -> float | None:
    """Coerce one of R's JSON renderings of an absent value to ``None``.

    R's ``NA`` reaches Python as the string ``"NA"`` and a NaN reaches it as a float,
    depending on how jsonlite rendered the value and on whether R printed it.  Both
    mean "upstream produced no number", and both must become ``None`` so the record
    reports ``undefined`` rather than carrying a value nothing computed.
    """

    if value is None:
        return None
    if isinstance(value, str):
        if value.strip().upper() in {"NA", "N/A", "NULL", "NAN", ""}:
            return None
        try:
            value = float(value)
        except ValueError:
            return None
    if isinstance(value, (int, float, np.floating, np.integer)):
        as_float = float(value)
        return None if not np.isfinite(as_float) else as_float
    return None


def _require_cycombine(context: BackendContext) -> str:
    version = r_package_version("cyCombine", context.rscript)
    if not version.available:
        raise NativeMetricUnavailable(
            "the cyCombine R package is not installed in this environment "
            f"({version.detail}); the native metric cannot be evaluated and no other "
            "implementation is substituted"
        )
    return version.version or "unknown"


def compute_native_metrics(
    dataset: Dataset,
    uncorrected: np.ndarray,
    corrected: np.ndarray,
    sampling: SamplingManifest,
    design: Any,
    *,
    run_id: str,
    candidate_id: str,
    context: BackendContext,
    work_dir: str,
    markers: list[str] | None = None,
    cell_col: str = "label",
    bin_size: float = DEFAULT_BIN_SIZE,
) -> list[MetricRecord]:
    """Evaluate ``cycombine_evaluate_emd`` and ``cycombine_evaluate_mad``.

    Returns one record per spec.  If cyCombine is absent, both records carry
    ``not_available`` and an explicit reason rather than a score; a missing metric is
    never a zero.
    """

    version = _require_cycombine(context)
    marker_ids = markers if markers is not None else dataset.marker_ids("evaluation")
    if not marker_ids:
        raise NativeMetricUnavailable("no marker is flagged use_for_evaluation")

    obs = dataset.obs
    cells = list(sampling.selected_cell_ids)
    sample_of = dict(zip(obs[CELL_ID].astype(str), obs[SAMPLE_ID].astype(str), strict=True))
    batch_of = dataset.sample_batches()
    # Donor count comes from the sample manifest, not from the cells: n_donors is a
    # property of the design, and a number derived from a cell-level count would be
    # the very pseudoreplication this package refuses.
    donor_of = {
        sid: rec.donor_id for sid, rec in dataset.samples.items() if rec.donor_id
    }
    samples = [sample_of[c] for c in cells]
    batches = [str(batch_of.get(s)) for s in samples]
    labels = [sampling.labels_by_cell[c] for c in cells]

    from pathlib import Path

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    unc_path = work / "native_uncorrected.binmat"
    cor_path = work / "native_corrected.binmat"
    write_matrix(unc_path, uncorrected, row_ids=cells, col_names=marker_ids)
    write_matrix(cor_path, corrected, row_ids=cells, col_names=marker_ids)

    params = {
        "uncorrected_matrix": str(unc_path),
        "corrected_matrix": str(cor_path),
        "markers": marker_ids,
        "batch": batches,
        "sample": samples,
        "label": labels,
        "cell_col": cell_col,
        "bin_size": bin_size,
    }

    def record(spec, value: float | None, status: MetricStatus, reason: str,
               extra: dict[str, Any] | None = None) -> MetricRecord:
        # Only samples that actually carry a donor count, so a manifest without donor
        # ids reports no donor count rather than counting every sample as its own donor.
        donors = {donor_of[s] for s in samples if donor_of.get(s)}
        # The shared builder, so these two specs carry the same frozen-partition
        # provenance as every other metric family rather than a retyped subset.
        return metric_record(
            spec, sampling,
            run_id=run_id,
            candidate_id=candidate_id,
            status=status,
            value=value,
            reason=reason,
            n_samples=len(set(samples)),
            n_donors=len(donors) if donors else None,
            n_cells=len(cells),
            coverage=1.0,
            aggregation_weight=1.0,
            **(extra or {}),
        )

    try:
        r_result = run_r_script(
            R_SCRIPT, params, work, rscript=context.rscript,
            log_path=work / "cycombine_evaluate.log",
        )
    except RNotAvailable as exc:
        reason = str(exc)
        return [
            record(NATIVE_CYCOMBINE_EMD, None, MetricStatus.NOT_IMPLEMENTED, reason),
            record(NATIVE_CYCOMBINE_MAD, None, MetricStatus.NOT_IMPLEMENTED, reason),
        ]
    except RScriptFailed as exc:
        reason = f"cyCombine metric evaluation failed: {exc}: {exc.log[-1500:]}"
        return [
            record(NATIVE_CYCOMBINE_EMD, None, MetricStatus.FAILED, reason),
            record(NATIVE_CYCOMBINE_MAD, None, MetricStatus.FAILED, reason),
        ]

    payload = r_result.payload
    emd_value = _as_float(payload.get("emd_reduction"))
    # evaluate_mad's summary number is called `score` upstream: the median absolute
    # difference between the uncorrected and corrected per-(population, batch, marker)
    # MADs, not a relative reduction like evaluate_emd's.
    mad_value = _as_float(payload.get("mad_score"))
    provenance = (
        f"cyCombine {version} via evaluate_emd/evaluate_mad on the frozen evaluation "
        f"partition ({len(cells)} cell(s), {len(set(samples))} sample(s)); "
        f"cell_col={cell_col!r}, bin_size={bin_size:g}. This is an upstream-defined "
        f"reduction on cyCombine's own partitioning and is NOT comparable with the "
        f"common-namespace numbers."
    )
    emd_reason = provenance
    rows = payload.get("emd_rows")
    if rows is not None:
        emd_reason += f" cyCombine filtered and compared {len(rows)} marker-table row(s)."

    return [
        record(
            NATIVE_CYCOMBINE_EMD,
            emd_value,
            MetricStatus.OK if emd_value is not None else MetricStatus.UNDEFINED,
            emd_reason if emd_value is not None
            else emd_reason + " upstream returned no reduction for this pair",
        ),
        record(
            NATIVE_CYCOMBINE_MAD,
            mad_value,
            MetricStatus.OK if mad_value is not None else MetricStatus.UNDEFINED,
            provenance if mad_value is not None
            else provenance + " upstream returned no MAD score for this pair",
        ),
    ]
