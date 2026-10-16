"""Stage entry points shared by the Snakemake rules and the command-line interface.

Each function reads its inputs from disk and writes one artifact, so the workflow's
declared dependencies are real dependencies: deleting an output forces recomputation,
and no downstream step infers success from a directory existing.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from batchdx.cyto._version import __version__
from batchdx.cyto.adapters import (
    BackendContext,
    BackendError,
    BackendIneligible,
    BackendNotAvailable,
    CorrectionBackend,
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
    CandidateSpecRecord,
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

    @property
    def baseline_dir(self) -> Path:
        """The frozen preprocessed baseline and everything it was written with."""
        return self.root / "baseline"

    @property
    def preprocess_baseline(self) -> Path:
        return self.baseline_dir / "preprocessed.h5ad"

    @property
    def removal_log(self) -> Path:
        """The per-event removal log of the preprocessing step (preprocess/removals.py)."""
        return self.baseline_dir / "removed_events.csv"

    @property
    def baseline_processing_manifest(self) -> Path:
        """The processing manifest of the baseline, not of the raw input."""
        return self.baseline_dir / "processing.json"

    @property
    def preprocess_summary(self) -> Path:
        """What each declared preprocessing step did, as written by the stage."""
        return self.baseline_dir / "preprocess_summary.json"


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
    # When preprocessing is enabled the dataset manifest is frozen by stage_preprocess
    # from the reloaded baseline, so the counts the design audit sees are the post-QC
    # ones.  Writing the raw-input manifest here as well would put two contradictory
    # manifests behind one path, and the last writer would win silently.
    if not config.preprocessing.enabled:
        dataset.manifest.write(layout.dataset_manifest)
    return layout.dataset_manifest

#: The expression layer the frozen baseline carries.  One name here, read by the
#: derived configuration below, so the layer is never re-typed at a call site.
BASELINE_LAYER = "preprocessed"


def _baseline_config(config: Config, layout: RunLayout) -> Config:
    """The configuration that reads the frozen baseline instead of the raw input.

    The baseline is an h5ad with its own processing manifest, and it is already
    transformed: the derived configuration declares no further preprocessing, because
    ``_validate_processing`` (io/dataset.py) refuses a run that would transform twice.
    That refusal is the contract working -- the baseline manifest records the
    transformation that was applied, and the derived configuration records that nothing
    further was.
    """

    derived = config.model_copy(deep=True)
    derived.input.format = "h5ad"
    derived.input.path = str(layout.preprocess_baseline)
    derived.input.expression_layer = BASELINE_LAYER
    derived.input.processing_manifest = str(layout.baseline_processing_manifest)
    derived.preprocessing.enabled = False
    derived.preprocessing.transformation = None
    return derived


def _load(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout | None = None
) -> tuple[Config, Dataset]:
    config = load_config(config_path)
    if layout is not None and config.preprocessing.enabled and layout.preprocess_baseline.exists():
        return config, load_dataset(_baseline_config(config, layout), base_dir)
    return config, load_dataset(config, base_dir)


def stage_preprocess(config_path: str | Path, base_dir: str | Path, layout: RunLayout) -> Path:
    """Apply the declared preprocessing chain to the raw input and freeze the baseline.

    The chain is the one the configuration declares -- preliminary 2-D gating
    (``preprocess/gating.py``), the per-panel asinh transform (``preprocess/transform.py``,
    with the cofactor for each panel taken from ``preprocessing.cofactor_by_panel``),
    and the PeacoQC pass (``preprocess/peacoqc.py``) when a settings file declares one --
    applied in that order, because a gate is evaluated on the events the previous step
    kept and PeacoQC is evaluated on the transformed values the run will analyse.

    Four artifacts are written under the run's ``baseline/`` directory: the preprocessed
    baseline itself (immutable, the input every later stage loads), the per-event
    removal log (``preprocess/removals.py``), the processing manifest *of the baseline*
    (declared scale, compensated/unmixed null with a note, transformation provenance),
    and a summary of what each step did.  The dataset manifest is then frozen from the
    reloaded baseline, so cell identity and the per-sample counts the design audit sees
    are the post-QC ones and not the raw ones.
    """

    import anndata as ad
    import pandas as pd

    from batchdx.cyto.io.fcs import read_fcs
    from batchdx.cyto.io.hashing import file_sha256
    from batchdx.cyto.io.manifests import (
        read_marker_table,
        read_processing_manifest,
        read_sample_table,
        write_removed_events,
    )
    from batchdx.cyto.preprocess.gating import apply_gates, gates_from_config
    from batchdx.cyto.preprocess.peacoqc import run_peacoqc, settings_from_config
    from batchdx.cyto.preprocess.removals import STEP_GATE, STEP_QC, empty_removal_log
    from batchdx.cyto.preprocess.transform import asinh_transform, resolve_convention

    config = load_config(config_path)
    if not config.preprocessing.enabled:
        raise ConfigError(
            "stage_preprocess: preprocessing.enabled is false; the stage exists only for "
            "a configuration that declares an upstream preprocessing step"
        )
    if config.input.format != "fcs":
        raise ConfigError(
            f"preprocessing: the declared input format is {config.input.format!r}; the "
            "preprocessing stage applies the declared chain to raw events, so it "
            "requires FCS input"
        )
    base = Path(base_dir)
    samples = {r.sample_id: r for r in read_sample_table(base / config.input.samples)}
    markers = read_marker_table(base / config.input.markers)
    raw_processing = read_processing_manifest(base / config.input.processing_manifest)
    marker_channels = [m.channel for m in markers]
    if any(c is None for c in marker_channels):
        unnamed = [
            m.feature_id for m, c in zip(markers, marker_channels, strict=True) if c is None
        ]
        raise ConfigError(
            f"marker manifest: feature(s) {unnamed} declare no channel, which FCS input "
            "requires; a marker without a channel cannot be read from raw events"
        )

    gate_table = gates_from_config(config, base_dir)
    qc_settings = settings_from_config(config, base_dir)
    convention = (
        resolve_convention(config.preprocessing.transformation)
        if config.preprocessing.transformation not in (None, "none") else None
    )
    cofactors = config.preprocessing.cofactor_by_panel
    panels = sorted({m.panel_id for m in markers if m.panel_id})
    if convention is not None:
        unpanelled = [m.feature_id for m in markers if not m.panel_id]
        if unpanelled:
            raise ConfigError(
                f"marker manifest: feature(s) {unpanelled} declare no panel_id; under "
                "the declared transformation every marker needs a panel, because a "
                "marker outside every declared panel would silently keep its raw scale "
                "in a baseline whose manifest says it was transformed"
            )
        missing = [p for p in panels if p not in cofactors]
        if missing:
            raise ConfigError(
                f"preprocessing.cofactor_by_panel: no cofactor is declared for panel(s) "
                f"{missing}; there is no universal cofactor, so the transform is refused "
                "rather than applied with a guessed one"
            )

    raw_dir = base / config.input.path
    blocks: list[np.ndarray] = []
    obs_rows: list[pd.DataFrame] = []
    removals: list[pd.DataFrame] = []
    gate_records: list[dict] = []
    qc_records: list[dict] = []
    qc_warnings: list[str] = []
    n_input_total = 0

    for sid, rec in samples.items():
        if not rec.source_file:
            raise ConfigError(
                f"sample {sid!r}: source_file is required for FCS input and is never "
                "inferred from file names"
            )
        path = raw_dir / rec.source_file
        contents = read_fcs(path)
        events = contents.events
        n = contents.n_events
        n_input_total += n
        index = np.arange(n, dtype=np.int64)
        ids = np.repeat(str(sid), n)
        keep = np.ones(n, dtype=bool)
        if gate_table is not None:
            result = apply_gates(
                gate_table, events, contents.channel_names,
                sample_ids=ids, original_event_index=index, step=STEP_GATE,
            )
            keep = result.keep
            if len(result.removed):
                removals.append(result.removed)
            gate_records.extend({
                "step": STEP_GATE, "detail": (
                    f"gate {o.gate.name} ({o.gate.kind} over "
                    f"{'/'.join(o.gate.channels)}): kept {o.kept} of {o.n_input} events "
                    f"({o.removed} removed)"
                ),
            } for o in result.outcomes)
        kept = events[keep]
        kept_index = index[keep]
        if convention is not None:
            transformed = np.array(kept, copy=True)
            for panel in panels:
                cols = [
                    contents.channel_names.index(str(markers[j].channel))
                    for j in range(len(markers)) if markers[j].panel_id == panel
                ]
                transformed[:, cols] = asinh_transform(
                    transformed[:, cols], float(cofactors[panel])
                )
            kept = transformed
        if qc_settings is not None:
            qc = run_peacoqc(
                kept, contents.channel_names, sample_ids=ids[keep],
                original_event_index=kept_index, settings=qc_settings,
                work_dir=layout.baseline_dir / "peacoqc" / str(sid),
                rscript=config.upstream.rscript,
            )
            if len(qc.removed):
                removals.append(qc.removed)
                removed_here = qc.removed["original_event_index"].to_numpy(np.int64)
                alive = ~np.isin(kept_index, removed_here)
                kept, kept_index = kept[alive], kept_index[alive]
            qc_records.extend({
                "step": STEP_QC, "detail": (
                    f"sample {o.sample_id}: kept {o.n_input - o.n_removed} of "
                    f"{o.n_input} events ({o.n_removed} removed); criteria applied: "
                    f"{', '.join(o.criteria_applied) or 'none'}"
                ),
            } for o in qc.samples)
            qc_warnings.extend(qc.warnings)

        cols = [contents.channel_names.index(str(c)) for c in marker_channels]
        blocks.append(kept[:, cols])
        checksum = file_sha256(path)
        obs_rows.append(pd.DataFrame({
            "sample_id": np.repeat(str(sid), kept.shape[0]),
            "original_event_index": kept_index,
            "acquisition_order": kept_index,
            "source_checksum": np.repeat(checksum, kept.shape[0]),
        }))

    matrix = np.vstack(blocks) if blocks else np.zeros((0, len(marker_channels)))
    obs = pd.concat(obs_rows, ignore_index=True)
    adata = ad.AnnData(
        X=np.zeros((matrix.shape[0], matrix.shape[1]), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index([m.feature_id for m in markers], dtype=object)),
    )
    adata.layers[BASELINE_LAYER] = matrix
    layout.preprocess_baseline.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(layout.preprocess_baseline)

    log = pd.concat(removals, ignore_index=True) if removals else empty_removal_log()
    write_removed_events(log, layout.removal_log)
    by_step = log.groupby("step").size().to_dict() if len(log) else {}

    scale = (
        "asinh (per-panel cofactors: "
        + ", ".join(f"{p}={float(cofactors[p]):g}" for p in panels) + ")"
        if convention is not None else str(raw_processing.scale)
    )
    processing = {
        "compensated": None,
        "unmixed": None,
        "compensation_source": None,
        "qc_history": [r["detail"] for r in gate_records + qc_records],
        "transformation": convention.name if convention is not None else None,
        "transformation_parameters": (
            {f"cofactor[{p}]": float(cofactors[p]) for p in panels}
            if convention is not None else {}
        ),
        "scale": scale,
        "units": (
            "arcsinh units (per-panel cofactors)" if convention is not None
            else raw_processing.units
        ),
        "corrected": False,
        "notes": (
            "written by batchdx.cyto.pipeline.stage_preprocess from the raw input "
            f"({config.input.path}); compensated/unmixed are null because the raw "
            f"processing manifest ({config.input.processing_manifest}) declares "
            f"compensated={raw_processing.compensated}, unmixed={raw_processing.unmixed} "
            "and this stage performs neither, so it does not re-assert them"
        ),
    }
    layout.baseline_processing_manifest.parent.mkdir(parents=True, exist_ok=True)
    layout.baseline_processing_manifest.write_text(
        json.dumps(processing, indent=2) + "\n", encoding="utf-8"
    )

    _write_json(layout.preprocess_summary, {
        "n_events_input": n_input_total,
        "n_events_kept": int(matrix.shape[0]),
        "n_events_removed": n_input_total - int(matrix.shape[0]),
        "n_samples": len(samples),
        "steps": gate_records + qc_records,
        "removals_by_step": {str(k): int(v) for k, v in by_step.items()},
        "declared": {
            "gating": config.preprocessing.gating,
            "transformation": config.preprocessing.transformation,
            "cofactor_by_panel": {str(k): float(v) for k, v in cofactors.items()},
            "peacoqc": config.preprocessing.peacoqc,
        },
        "processing": processing,
        "warnings": qc_warnings,
    })

    # The manifest is frozen from the reloaded baseline, not from the in-memory frame,
    # so cell identity and the per-sample counts are exactly what every later stage
    # will load -- the same guarantee stage_validate gives the raw input.
    dataset = load_dataset(_baseline_config(config, layout), base_dir)
    dataset.manifest.write(layout.dataset_manifest)
    return layout.dataset_manifest


def stage_design(
    config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Audit the design and emit the eligibility manifest."""

    config, dataset = _load(config_path, base_dir, layout)
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

    config, dataset = _load(config_path, base_dir, layout)
    design = DesignAudit.read(layout.design_audit)
    labels, tier, confidence = resolve_labels(dataset, config, base_dir)
    manifest = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    manifest.write(layout.sampling_manifest)
    return layout.sampling_manifest


def stage_diagnose(
    detector_name: str, config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Run one detector and write its result, including a scientific skip."""

    config, dataset = _load(config_path, base_dir, layout)
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
    dataset: Dataset, layout: RunLayout, candidate_id: str,
    record: CandidateSpecRecord, status: CorrectionStatus, message: str,
    context: BackendContext, backend: CorrectionBackend,
    extra_warnings: list[str] | None = None,
) -> Path:
    """Record a failed/ineligible candidate and still write its declared output.

    The ``correct`` rule declares ``corrected.h5ad`` as an output; for a candidate
    that never produced a correction we write the unmodified baseline so the
    declared output always exists. The result records the failure status, so no
    consumer mistakes this for a correction. ``extra_warnings`` carries notes that
    belong on the result whatever the outcome -- an exploratory-override waiver was
    granted before the backend ran, so a candidate that then failed is still a
    candidate whose result must say it ran under a bypassed design gate.
    """

    warnings: list[str] = []
    try:
        software = backend.availability(context)
    except (BackendError, ValueError, OSError, subprocess.TimeoutExpired) as probe_error:
        # A candidate that failed must still be recorded even if the availability
        # probe itself is what is broken; the probe failure is preserved as a warning
        # and never replaces the candidate's own reason.
        software = []
        warnings.append(f"software availability could not be recorded: {probe_error}")
    result = failure_result(record, status, message, software=software)
    result.warnings = warnings + (extra_warnings or [])
    result.write(layout.correction_result(candidate_id))
    write_candidate_output(
        dataset, dataset.expression(), layout.corrected_h5ad(candidate_id),
        candidate_id=candidate_id,
    )
    return layout.correction_result(candidate_id)


def stage_correct(
    candidate_id: str, config_path: str | Path, base_dir: str | Path, layout: RunLayout
) -> Path:
    """Run one correction candidate from the immutable baseline."""

    config, dataset = _load(config_path, base_dir, layout)
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

    # A retry invalidates its previous success before anything can fail or be cancelled.
    for path in (layout.correction_result(candidate_id), layout.corrected_h5ad(candidate_id),
                 layout.metric_result(candidate_id)):
        path.unlink(missing_ok=True)
    log_path = layout.root / "logs" / f"correct_{candidate_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"Starting candidate {candidate_id!r} ({spec.backend}).\n", encoding="utf-8")

    eligibility = next(
        (e for e in design.eligibility if e.name == spec.backend and e.kind == "backend"), None
    )
    waived: list[dict[str, str]] = []
    override_note: str | None = None
    if eligibility is not None and not eligibility.eligible:
        if config.design.allow_exploratory_override and eligibility.override_allowed:
            # The gate is bypassed, not satisfied, and the bypass is recorded before the
            # backend runs so that a crash downstream cannot erase the fact that it was
            # attempted.  Three places carry it: this audit file, which stage_decide
            # turns into ``DecisionResult.overrides`` for the report; the waiver on the
            # correction result, so a reader of that artifact alone still sees it; and
            # the warning below, in the same words for both.
            waived = [{
                "candidate_id": candidate_id,
                "backend": spec.backend,
                "reason_code": eligibility.reason_code or "unknown",
                "reason": eligibility.reason or "",
                "allowed_by": "design.allow_exploratory_override",
            }]
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
            return _write_failure_candidate(
                dataset, layout, candidate_id, record, CorrectionStatus.INELIGIBLE,
                f"{eligibility.reason} (reason code: {eligibility.reason_code})",
                context, backend,
            )

    # Only expected candidate execution/output errors are reportable candidate
    # failures; shared input loading above and cancellation remain whole-run failures.
    # OutputValidationError (a ValueError): the backend produced a non-finite or
    # reshaped matrix.  FileNotFoundError (an OSError): the R script never wrote its
    # output.  JSONDecodeError (a ValueError): the R script's result payload is
    # unparseable.  TimeoutExpired: the R script overran its budget.
    candidate_errors = (
        BackendError, OutputValidationError, ValueError, OSError, subprocess.TimeoutExpired,
    )

    def _record_for_exc(exc: Exception) -> CandidateSpecRecord:
        try:
            resolved = resolve_parameters(backend.capabilities, dict(spec.parameters))
        except BackendError:
            resolved = dict(spec.parameters)
        return make_spec(candidate_id, spec.backend, None, resolved,
                         config.project.seed, dataset, [])

    try:
        matrix, result = backend.run(dataset, candidate_id, dict(spec.parameters), context)
    except BackendNotAvailable as exc:
        # Availability is a property of the environment, decided by the exception
        # type, never by sniffing the error text.
        with log_path.open("a", encoding="utf-8") as log:
            log.write(traceback.format_exc())
        return _write_failure_candidate(
            dataset, layout, candidate_id, _record_for_exc(exc),
            CorrectionStatus.NOT_AVAILABLE, str(exc), context, backend,
            extra_warnings=[override_note] if override_note else None,
        )
    except BackendIneligible as exc:
        # A design-coverage refusal is a statement about the data, recorded as
        # INELIGIBLE rather than NOT_AVAILABLE: the package is installed, the split
        # simply cannot support the fit.  The branch precedes the generic
        # candidate_errors branch, which would otherwise file the same refusal as
        # FAILED because BackendIneligible is a BackendError.
        with log_path.open("a", encoding="utf-8") as log:
            log.write(traceback.format_exc())
        return _write_failure_candidate(
            dataset, layout, candidate_id, _record_for_exc(exc),
            CorrectionStatus.INELIGIBLE, str(exc), context, backend,
            extra_warnings=[override_note] if override_note else None,
        )
    except candidate_errors as exc:
        with log_path.open("a", encoding="utf-8") as log:
            log.write(traceback.format_exc())
        return _write_failure_candidate(
            dataset, layout, candidate_id, _record_for_exc(exc),
            CorrectionStatus.FAILED, str(exc), context, backend,
            extra_warnings=[override_note] if override_note else None,
        )

    if result.status is not CorrectionStatus.SUCCESS:
        # The backend reported its own failure and returned no matrix; record that
        # failure rather than exporting whatever came back.
        return _write_failure_candidate(
            dataset, layout, candidate_id, result.candidate,
            result.status, result.reason or "the backend reported failure",
            context, backend,
            extra_warnings=[override_note] if override_note else None,
        )

    try:
        path = write_candidate_output(
            dataset, matrix, layout.corrected_h5ad(candidate_id),
            candidate_id=candidate_id, output_scale=result.output_scale,
        )
    except candidate_errors as exc:
        # A partial export is removed: a file that is half an H5AD must never be
        # readable as a corrected output.  The failure candidate writes the
        # unmodified baseline in its place, so the declared output exists and the
        # result records that it is not a correction.
        with log_path.open("a", encoding="utf-8") as log:
            log.write(traceback.format_exc())
        layout.corrected_h5ad(candidate_id).unlink(missing_ok=True)
        return _write_failure_candidate(
            dataset, layout, candidate_id, _record_for_exc(exc),
            CorrectionStatus.FAILED, str(exc), context, backend,
            extra_warnings=[override_note] if override_note else None,
        )

    result.output_path = str(path)
    result.waived_eligibility = waived
    if waived:
        result.warnings = list(result.warnings) + [
            f"EXPLORATORY OVERRIDE: {w['backend']} ran on a design its own eligibility "
            f"rule rejects ({w['reason'] or 'unknown'}; reason code: {w['reason_code']}), "
            "because design.allow_exploratory_override = true. Its output is exploratory "
            "and cannot support a recommendation."
            for w in waived
        ]
    result.write(layout.correction_result(candidate_id))
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"Status: {result.status.value}. {result.reason or 'completed'}\n")
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

    config, dataset = _load(config_path, base_dir, layout)
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

    # The sample table is re-read here only to answer the held-out question, which is a
    # question about which samples sit in which split.  It is not part of the frozen
    # manifest, and passing it keeps the decision layer free of any dataset handle.
    _, dataset = _load(config_path, base_dir, layout)

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
    # override audit trail promised by the assessment layer exists in real runs. Read
    # from the audit file rather than rebuilt from the correction results: the file is
    # written before the backend runs, so a candidate that crashed mid-correction is
    # still accounted for here.
    overrides = _read_exploratory_overrides(layout)

    decision = assess_candidates(
        design, sampling, corrections, metrics, reductions,
        samples=dict(dataset.samples),
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
