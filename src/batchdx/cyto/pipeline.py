"""Stage entry points shared by the Snakemake rules and the command-line interface.

Each function reads its inputs from disk and writes one artifact, so the workflow's
declared dependencies are real dependencies: deleting an output forces recomputation,
and no downstream step infers success from a directory existing.
"""

from __future__ import annotations

import json
import platform
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from batchdx.cyto._version import __version__
from batchdx.cyto.adapters import (
    BackendContext,
    BackendError,
    BackendNotAvailable,
    get_backend,
)
from batchdx.cyto.config import Config, load_config
from batchdx.cyto.decisions import assess_candidates
from batchdx.cyto.design import audit_design
from batchdx.cyto.diagnostics import DetectorContext, get_detector
from batchdx.cyto.io import load_dataset, write_candidate_output
from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.metrics import compute_common_metrics, compute_reductions
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.schemas import export_json_schemas
from batchdx.cyto.schemas.contracts import (
    CorrectionResult,
    DatasetManifest,
    DesignAudit,
    MetricResult,
    SamplingManifest,
    StageStatus,
)
from batchdx.cyto.schemas.enums import CorrectionStatus, MetricStatus

DETECTOR_FLAGS = {
    "acquisition_drift": "acquisition_drift",
    "composition": "composition",
    "cycombine_express": "cycombine_express",
    "cycombine_cluster": "cycombine_cluster",
    "cytobatchflagr": "cytobatchflagr",
}


@dataclass(frozen=True)
class RunLayout:
    """Where every artifact of one run lives."""

    root: Path

    @property
    def resolved_config(self) -> Path:
        return self.root / "config.resolved.json"

    @property
    def dataset_manifest(self) -> Path:
        return self.root / "manifests" / "dataset_manifest.json"

    @property
    def schemas_dir(self) -> Path:
        return self.root / "manifests" / "schemas"

    @property
    def design_audit(self) -> Path:
        return self.root / "design" / "design_audit.json"

    @property
    def sampling_manifest(self) -> Path:
        return self.root / "manifests" / "sampling_manifest.json"

    @property
    def decision(self) -> Path:
        return self.root / "decisions" / "decision_result.json"

    @property
    def exploratory_overrides(self) -> Path:
        """The record of every design gate this run bypassed by configuration."""

        return self.root / "decisions" / "exploratory_overrides.json"

    @property
    def reductions(self) -> Path:
        return self.root / "metrics" / "reductions.json"

    @property
    def report(self) -> Path:
        return self.root / "report" / "report.html"

    @property
    def status(self) -> Path:
        return self.root / "status" / "run_status.json"

    def detector_result(self, name: str) -> Path:
        return self.root / "diagnostics" / name / "detector_result.json"

    def candidate_dir(self, candidate_id: str) -> Path:
        return self.root / "candidates" / candidate_id

    def correction_result(self, candidate_id: str) -> Path:
        return self.candidate_dir(candidate_id) / "correction_result.json"

    def corrected_h5ad(self, candidate_id: str) -> Path:
        return self.candidate_dir(candidate_id) / "corrected.h5ad"

    def metric_result(self, candidate_id: str) -> Path:
        return self.root / "metrics" / candidate_id / "metric_result.json"


def environment_report() -> dict:
    return {
        "batchdx_cyto": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }


def _write_json(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return path


def _read_exploratory_overrides(layout: RunLayout) -> list[dict[str, str]]:
    """Every recorded design-gate bypass in this run, oldest first.

    Read here rather than inferred, so the decision layer can carry the fact that a
    candidate ran on a design the audit called ineligible.  A run that bypassed nothing
    has no file and reports no overrides, which is the same answer as a run whose only
    overrides are absent from disk.
    """

    path = layout.exploratory_overrides
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    return [dict(entry) for entry in payload.get("overrides", [])]


def _record_exploratory_override(layout: RunLayout, entry: dict[str, str]) -> str:
    """Add one bypass to the audit trail and return the sentence to warn with.

    Keyed by candidate so that re-running ``correct`` for the same candidate replaces its
    record instead of appending a second, contradictory one: the file answers "was this
    candidate run under an override in this run", not "how many times did the stage
    execute".
    """

    existing = [
        e for e in _read_exploratory_overrides(layout)
        if e.get("candidate_id") != entry["candidate_id"]
    ]
    existing.append(entry)
    _write_json(layout.exploratory_overrides, {"overrides": existing})
    return (
        f"EXPLORATORY OVERRIDE: {entry['backend']} ran on a design its own eligibility "
        f"rule rejects ({entry['overrode']}), because {entry['declared_by']}. Its output "
        "is exploratory and cannot support a recommendation."
    )


# --------------------------------------------------------------------------------------
# Stages
# --------------------------------------------------------------------------------------


def stage_validate(config_path: str | Path, base_dir: str | Path, layout: RunLayout) -> Path:
    """Validate the configuration and the input, then freeze the dataset manifest."""

    config = load_config(config_path)
    dataset = load_dataset(config, base_dir)
    dataset.assert_baseline_unchanged()

    _write_json(layout.resolved_config, {
        "config": config.resolved_dict(),
        "analytical_fingerprint": config.analytical_fingerprint(),
        "environment": environment_report(),
    })
    export_json_schemas(layout.schemas_dir)
    dataset.manifest.write(layout.dataset_manifest)
    return layout.dataset_manifest


def _load(config_path: str | Path, base_dir: str | Path) -> tuple[Config, Dataset]:
    config = load_config(config_path)
    return config, load_dataset(config, base_dir)


def stage_design(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Audit the design and emit the eligibility manifest."""

    config, dataset = _load(config_path, base_dir)
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
    design.write(layout.design_audit)
    return layout.design_audit


def stage_partition(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Freeze evaluation labels, sampling and the common comparison set."""

    config, dataset = _load(config_path, base_dir)
    design = DesignAudit.read(layout.design_audit)
    labels, tier, confidence = resolve_labels(dataset, config, base_dir)
    manifest = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    manifest.write(layout.sampling_manifest)
    return layout.sampling_manifest


def stage_diagnose(
    detector_name: str, config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Run one detector and write its result, including a scientific skip."""

    config, dataset = _load(config_path, base_dir)
    design = DesignAudit.read(layout.design_audit)
    sampling = (
        SamplingManifest.read(layout.sampling_manifest)
        if layout.sampling_manifest.exists() else None
    )
    out_dir = layout.root / "diagnostics"
    context = DetectorContext(
        run_id=config.project.id, design=design, sampling=sampling, out_dir=str(out_dir),
        seed=config.project.seed, rscript=config.upstream.rscript,
        upstream_dir=config.upstream.checkout_dir,
        threshold_path=(
            str(Path(base_dir) / config.evaluation.gating_thresholds)
            if config.evaluation.gating_thresholds else None
        ),
    )
    result = get_detector(detector_name).run(dataset, context)
    result.write(layout.detector_result(detector_name))
    return layout.detector_result(detector_name)


def stage_correct(
    candidate_id: str, config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Run one correction candidate from the immutable baseline."""

    config, dataset = _load(config_path, base_dir)
    design = DesignAudit.read(layout.design_audit)
    sampling = (
        SamplingManifest.read(layout.sampling_manifest)
        if layout.sampling_manifest.exists() else None
    )
    spec = next((c for c in config.correction.candidates if c.id == candidate_id), None)
    if spec is None:
        raise KeyError(f"candidate {candidate_id!r} is not in the configuration")

    from batchdx.cyto.adapters.base import failure_result, make_spec, resolve_parameters

    context = BackendContext(
        design=design, sampling=sampling,
        work_dir=str(layout.candidate_dir(candidate_id) / "work"),
        seed=config.project.seed, rscript=config.upstream.rscript,
        upstream_dir=config.upstream.checkout_dir,
    )
    backend = get_backend(spec.backend)

    eligibility = next(
        (e for e in design.eligibility if e.name == spec.backend and e.kind == "backend"), None
    )
    override_note: str | None = None
    if eligibility is not None and not eligibility.eligible:
        if config.design.allow_exploratory_override and eligibility.override_allowed:
            # The gate is bypassed, not satisfied, and the bypass is recorded before the
            # backend runs so that a crash downstream cannot erase the fact that it was
            # attempted.  Two places carry it: this audit file, which stage_decide turns
            # into ``DecisionResult.overrides`` for the report, and a warning on the
            # correction result itself, so a reader of that artifact alone still sees it.
            override_note = _record_exploratory_override(layout, {
                "candidate_id": candidate_id,
                "backend": spec.backend,
                "overrode": f"{eligibility.reason} (reason code: {eligibility.reason_code})",
                "declared_by": "design.allow_exploratory_override = true",
                "effect": (
                    "the backend was run on a design its own eligibility rule rejects. "
                    "Its output is exploratory: it is scored and displayed, and it cannot "
                    "support a recommendation"
                ),
            })
        else:
            record = make_spec(candidate_id, spec.backend, None, dict(spec.parameters),
                               config.project.seed, dataset, [])
            failure_result(
                record, CorrectionStatus.INELIGIBLE,
                f"{eligibility.reason} (reason code: {eligibility.reason_code})",
            ).write(layout.correction_result(candidate_id))
            return layout.correction_result(candidate_id)

    def finish(result: CorrectionResult) -> Path:
        if override_note is not None:
            result.warnings = [*result.warnings, override_note]
        return result.write(layout.correction_result(candidate_id))

    try:
        matrix, result = backend.run(dataset, candidate_id, dict(spec.parameters), context)
    except BackendNotAvailable as exc:
        # Ordered before the BackendError branch because it is a subclass, and it must
        # be: the environment-gap and method-failure states are what
        # docs/decisions/0006-candidate-results-and-failure-handling.md requires a
        # report to keep apart, and they are told apart by type rather than by
        # matching words in a message.
        record = make_spec(candidate_id, spec.backend, None, dict(spec.parameters),
                           config.project.seed, dataset, [])
        finish(failure_result(record, CorrectionStatus.NOT_AVAILABLE, str(exc),
                              software=backend.availability(context)))
        return layout.correction_result(candidate_id)
    except BackendError as exc:
        try:
            resolved = resolve_parameters(backend.capabilities, dict(spec.parameters))
        except BackendError:
            resolved = dict(spec.parameters)
        record = make_spec(candidate_id, spec.backend, None, resolved,
                           config.project.seed, dataset, [])
        finish(failure_result(record, CorrectionStatus.FAILED, str(exc),
                              software=backend.availability(context)))
        return layout.correction_result(candidate_id)

    path = write_candidate_output(
        dataset, matrix, layout.corrected_h5ad(candidate_id),
        candidate_id=candidate_id, output_scale=result.output_scale,
    )
    result.output_path = str(path)
    finish(result)
    return layout.correction_result(candidate_id)


def _read_candidate_matrix(
    dataset: Dataset, layout: RunLayout, candidate_id: str, result: CorrectionResult
) -> np.ndarray | None:
    if result.status is not CorrectionStatus.SUCCESS:
        return None
    path = layout.corrected_h5ad(candidate_id)
    if not path.exists():
        return None
    import anndata as ad

    adata = ad.read_h5ad(path)
    layer = str(adata.uns.get("expression_layer", "corrected"))
    matrix = np.asarray(adata.layers[layer], dtype=np.float64)
    if list(map(str, adata.obs["cell_id"])) != list(map(str, dataset.adata.obs["cell_id"])):
        raise ValueError(
            f"candidate {candidate_id!r}: the stored output does not carry the baseline's "
            "cell ids in order"
        )
    return matrix


def stage_metrics(
    candidate_id: str, config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Score one candidate on the frozen common partition."""

    config, dataset = _load(config_path, base_dir)
    design = DesignAudit.read(layout.design_audit)
    sampling = SamplingManifest.read(layout.sampling_manifest)
    correction = CorrectionResult.read(layout.correction_result(candidate_id))
    matrix = _read_candidate_matrix(dataset, layout, candidate_id, correction)

    if matrix is None:
        MetricResult(
            run_id=config.project.id, candidate_id=candidate_id, metric_group="common",
            status="ineligible",
            reason=(
                f"candidate {candidate_id!r} has status {correction.status.value}: "
                f"{correction.reason or 'no corrected output was produced'}. "
                "A missing metric is not a zero score."
            ),
        ).write(layout.metric_result(candidate_id))
        return layout.metric_result(candidate_id)

    result = compute_common_metrics(
        dataset, matrix, sampling, design,
        run_id=config.project.id, candidate_id=candidate_id,
    )

    # The native namespace is optional by construction: it calls upstream R and cannot
    # run without cyCombine.  A failure here is recorded as a status on those two specs
    # and never removes the common-namespace records, so an environment gap costs the
    # native numbers and nothing else.  Only the identity candidate is compared against
    # itself, because a native reduction is a ratio of the corrected and uncorrected
    # states and the identity is the only candidate whose uncorrected state is itself.
    if candidate_id == _identity_candidate_id(config):
        try:
            from batchdx.cyto.metrics.native import compute_native_metrics

            result.records.extend(compute_native_metrics(
                dataset, dataset.expression(), matrix, sampling, design,
                run_id=config.project.id, candidate_id=candidate_id,
                context=BackendContext(
                    design=design, sampling=sampling,
                    work_dir=str(layout.candidate_dir(candidate_id) / "native"),
                    seed=config.project.seed, rscript=config.upstream.rscript,
                    upstream_dir=config.upstream.checkout_dir,
                ),
                work_dir=str(layout.candidate_dir(candidate_id) / "native"),
            ))
        except Exception as exc:  # noqa: BLE001 - a status is the point
            from batchdx.cyto.metrics.records import metric_record
            from batchdx.cyto.metrics.registry import NATIVE_CYCOMBINE_EMD

            # One record, carrying the same fields every other metric record carries, so
            # the failure is visible in the same place a value would have been and cannot
            # be mistaken for an omitted metric.
            result.records.append(metric_record(
                NATIVE_CYCOMBINE_EMD, sampling,
                run_id=config.project.id,
                candidate_id=candidate_id,
                status=MetricStatus.FAILED,
                reason=f"the native-namespace metrics could not be produced: {exc}",
                n_samples=0, n_cells=0, coverage=0.0, aggregation_weight=1.0,
            ))

    result.write(layout.metric_result(candidate_id))
    return layout.metric_result(candidate_id)


def _identity_candidate_id(config: Config) -> str:
    for candidate in config.correction.candidates:
        if candidate.backend == "none":
            return candidate.id
    raise KeyError("the configuration has no identity candidate")


def stage_decide(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Aggregate, compute reductions, and assign a state to every candidate."""

    config = load_config(config_path)
    design = DesignAudit.read(layout.design_audit)
    sampling = SamplingManifest.read(layout.sampling_manifest)

    # The sample table is re-read here only to answer the held-out question, which is a
    # question about which samples sit in which split.  It is not part of the frozen
    # manifest, and passing it keeps the decision layer free of any dataset handle.
    _, dataset = _load(config_path, base_dir)

    corrections: dict[str, CorrectionResult] = {}
    metrics: dict[str, MetricResult] = {}
    for candidate in config.correction.candidates:
        path = layout.correction_result(candidate.id)
        if path.exists():
            corrections[candidate.id] = CorrectionResult.read(path)
        mpath = layout.metric_result(candidate.id)
        if mpath.exists():
            metrics[candidate.id] = MetricResult.read(mpath)

    identity = _identity_candidate_id(config)
    reductions = compute_reductions(
        metrics, sampling, run_id=config.project.id, identity_candidate_id=identity
    )
    _write_json(layout.reductions, {"records": [r.model_dump(mode="json") for r in reductions]})

    decision = assess_candidates(
        design, sampling, corrections, metrics, reductions,
        samples=dict(dataset.samples),
        run_id=config.project.id, policy_name=config.evaluation.policy,
        identity_candidate_id=identity,
        overrides=_read_exploratory_overrides(layout),
    )
    decision.write(layout.decision)
    return layout.decision


def stage_status(config_path: str | Path, layout: RunLayout) -> Path:
    """Summarise what ran, what was skipped, and why.  Never inferred from directories."""

    config = load_config(config_path)
    entries: list[StageStatus] = []

    for name, flag in DETECTOR_FLAGS.items():
        requested = getattr(config.detectors, flag)
        path = layout.detector_result(name)
        if not requested:
            entries.append(StageStatus(stage="diagnose", unit=name, status="skipped",
                                       reason="not requested in the configuration"))
        elif not path.exists():
            entries.append(StageStatus(stage="diagnose", unit=name, status="failed",
                                       reason="the detector produced no result file"))
        else:
            from batchdx.cyto.schemas.contracts import DetectorResult

            result = DetectorResult.read(path)
            mapped = {"ok": "ok", "ineligible": "ineligible", "failed": "failed",
                      "not_available": "blocked"}[result.status]
            entries.append(StageStatus(stage="diagnose", unit=name, status=mapped,
                                       reason=result.reason or "completed"))

    for candidate in config.correction.candidates:
        path = layout.correction_result(candidate.id)
        if not path.exists():
            entries.append(StageStatus(stage="correct", unit=candidate.id, status="failed",
                                       reason="no correction result file was written"))
            continue
        result = CorrectionResult.read(path)
        mapped = {
            CorrectionStatus.SUCCESS: "ok",
            CorrectionStatus.INELIGIBLE: "ineligible",
            CorrectionStatus.FAILED: "failed",
            CorrectionStatus.NOT_AVAILABLE: "blocked",
        }[result.status]
        entries.append(StageStatus(stage="correct", unit=candidate.id, status=mapped,
                                   reason=result.reason or "completed"))

    _write_json(layout.status, {
        "run_id": config.project.id,
        "analytical_fingerprint": config.analytical_fingerprint(),
        "environment": environment_report(),
        "entries": [e.model_dump(mode="json") for e in entries],
    })
    return layout.status


def stage_report(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout, diagnosis_only: bool = False
) -> Path:
    """Render the offline HTML report from the structured artifacts already on disk."""

    from batchdx.cyto.reporting import render_report

    return render_report(config_path, base_dir, layout, diagnosis_only=diagnosis_only)


def read_dataset_manifest(layout: RunLayout) -> DatasetManifest:
    return DatasetManifest.read(layout.dataset_manifest)
