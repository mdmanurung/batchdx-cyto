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
    OutputValidationError,
    failure_result,
    get_backend,
    make_spec,
    resolve_parameters,
)
from batchdx.cyto.config import Config, ConfigError, load_config
from batchdx.cyto.decisions import assess_candidates
from batchdx.cyto.design import audit_design
from batchdx.cyto.diagnostics import DetectorContext, get_detector
from batchdx.cyto.io import load_dataset, write_candidate_output
from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.metrics import (
    compute_biological_contrasts,
    compute_common_metrics,
    compute_reductions,
)
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
from batchdx.cyto.schemas.enums import CorrectionStatus

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


# --------------------------------------------------------------------------------------
# Stages
# --------------------------------------------------------------------------------------


def stage_validate(config_path: str | Path, base_dir: str | Path, layout: RunLayout) -> Path:
    """Validate the configuration and the input, then freeze the dataset manifest."""

    config = load_config(config_path)
    dataset = load_dataset(config, base_dir)
    dataset.assert_baseline_unchanged()

    # A declared upstream lock is a real input: if it is named, it must exist and
    # name the pinned commits the adapters were verified against. An unread lock
    # would let a configuration claim a pin it cannot enforce.
    if config.upstream.lock:
        lock_path = Path(base_dir) / config.upstream.lock
        if not lock_path.exists():
            raise ConfigError(
                f"upstream.lock names {config.upstream.lock!r} but the file does not "
                f"exist under {base_dir}; run scripts/fetch_upstream.sh or remove the key"
            )
        try:
            import yaml

            lock_data = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            raise ConfigError(
                f"upstream.lock {config.upstream.lock!r} is not valid YAML: {exc}"
            ) from exc
        if not isinstance(lock_data, dict) or not lock_data:
            raise ConfigError(
                f"upstream.lock {config.upstream.lock!r} must be a non-empty mapping of "
                "tool name to pinned commit"
            )

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


def _write_failure_candidate(
    dataset: Dataset, layout: RunLayout, candidate_id: str, record: CorrectionSpec,
    status: CorrectionStatus, message: str, context: BackendContext,
    backend: CorrectionBackend,
) -> Path:
    """Record a failed/ineligible candidate and still write its declared output.

    The ``correct`` rule declares ``corrected.h5ad`` as an output; for a candidate
    that never produced a correction we write the unmodified baseline so the
    declared output always exists. The result records the failure status, so no
    consumer mistakes this for a correction.
    """

    failure_result(record, status, message, software=backend.availability(context)).write(
        layout.correction_result(candidate_id)
    )
    write_candidate_output(
        dataset, dataset.expression(), layout.corrected_h5ad(candidate_id),
        candidate_id=candidate_id,
    )
    return layout.correction_result(candidate_id)


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
    waived: list[dict[str, str]] = []
    if eligibility is not None and not eligibility.eligible:
        if config.design.allow_exploratory_override and eligibility.override_allowed:
            # The waiver is recorded on the correction result and on the decision, so
            # an exploratory run is never indistinguishable from a fully eligible one.
            waived = [{
                "candidate_id": candidate_id,
                "backend": spec.backend,
                "reason_code": eligibility.reason_code or "unknown",
                "reason": eligibility.reason or "",
                "allowed_by": "design.allow_exploratory_override",
            }]
        else:
            record = make_spec(candidate_id, spec.backend, None, dict(spec.parameters),
                               config.project.seed, dataset, [])
            return _write_failure_candidate(
                dataset, layout, candidate_id, record, CorrectionStatus.INELIGIBLE,
                f"{eligibility.reason} (reason code: {eligibility.reason_code})",
                context, backend,
            )

    try:
        matrix, result = backend.run(dataset, candidate_id, dict(spec.parameters), context)
    except BackendNotAvailable as exc:
        # Availability is a property of the environment, decided by the exception
        # type, never by sniffing the error text.
        try:
            resolved = resolve_parameters(backend.capabilities, dict(spec.parameters))
        except BackendError:
            resolved = dict(spec.parameters)
        record = make_spec(candidate_id, spec.backend, None, resolved,
                           config.project.seed, dataset, [])
        return _write_failure_candidate(
            dataset, layout, candidate_id, record, CorrectionStatus.NOT_AVAILABLE,
            str(exc), context, backend,
        )
    except (BackendError, OutputValidationError, FileNotFoundError) as exc:
        # OutputValidationError: the backend produced a non-finite or reshaped
        # matrix; FileNotFoundError: the R script never wrote its output. Both are
        # candidate failures and must be recorded, not crash the whole workflow.
        try:
            resolved = resolve_parameters(backend.capabilities, dict(spec.parameters))
        except BackendError:
            resolved = dict(spec.parameters)
        record = make_spec(candidate_id, spec.backend, None, resolved,
                           config.project.seed, dataset, [])
        return _write_failure_candidate(
            dataset, layout, candidate_id, record, CorrectionStatus.FAILED,
            str(exc), context, backend,
        )

    path = write_candidate_output(
        dataset, matrix, layout.corrected_h5ad(candidate_id),
        candidate_id=candidate_id, output_scale=result.output_scale,
    )
    result.output_path = str(path)
    result.waived_eligibility = waived
    if waived:
        result.warnings = list(result.warnings) + [
            "eligibility constraint(s) waived for this candidate: "
            + "; ".join(f"{w['reason_code']} ({w['reason']})" for w in waived)
        ]
    result.write(layout.correction_result(candidate_id))
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
    if config.design.biological_contrasts:
        result.records.extend(compute_biological_contrasts(
            dataset, dataset.expression(), matrix, sampling, design,
            config.design.biological_contrasts,
            run_id=config.project.id, candidate_id=candidate_id,
            donor_factor=config.design.donor,
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

    # Surface any waived eligibility constraints on the decision record, so the
    # override audit trail promised by the assessment layer exists in real runs.
    overrides = [
        {"candidate_id": w["candidate_id"], "reason_code": w["reason_code"],
         "reason": w["reason"], "allowed_by": w["allowed_by"]}
        for c in corrections.values() for w in (c.waived_eligibility or [])
    ]

    decision = assess_candidates(
        design, sampling, corrections, metrics, reductions,
        run_id=config.project.id, policy_name=config.evaluation.policy,
        identity_candidate_id=identity, overrides=overrides,
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
