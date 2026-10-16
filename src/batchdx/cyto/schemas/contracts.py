"""The seven cross-stage contracts.

Every record is versioned with :data:`CONTRACT_SCHEMA_VERSION` and forbids unknown
keys, so a record written by a newer build fails loudly instead of being partially
read.  JSON Schemas for the same models are exported for the R side by
:func:`export_json_schemas`, which is what makes the Python/R boundary round-trip
testable rather than merely conventional.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from batchdx.cyto._version import CONTRACT_SCHEMA_VERSION
from batchdx.cyto.schemas.enums import (
    ComparisonSupport,
    CorrectionStatus,
    DecisionState,
    LabelTier,
    MarkerRole,
    MetricDirection,
    MetricStatus,
    ReferenceRole,
)


class _Contract(BaseModel):
    """Base class: strict, versioned, and stable under JSON round-trip."""

    model_config = ConfigDict(extra="forbid", use_enum_values=False, validate_assignment=True)

    schema_version: int = Field(default=CONTRACT_SCHEMA_VERSION)

    def to_json(self, indent: int = 2) -> str:
        return self.model_dump_json(indent=indent)

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json() + "\n", encoding="utf-8")
        return path

    @classmethod
    def read(cls, path: str | Path):
        return cls.model_validate_json(Path(path).read_text(encoding="utf-8"))


# --------------------------------------------------------------------------------------
# Provenance building blocks
# --------------------------------------------------------------------------------------


class ProvenanceRef(_Contract):
    """Pointer from a derived record back to the artifact it was computed from."""

    kind: str
    identifier: str
    content_hash: str | None = None
    path: str | None = None


class SoftwareVersion(_Contract):
    """One resolved dependency, as observed at run time (never assumed)."""

    name: str
    version: str | None = None
    source: str | None = None
    commit: str | None = None
    available: bool = True
    detail: str | None = None


class StageStatus(_Contract):
    """Machine-readable status for any unit of work that can legitimately not run."""

    stage: str
    unit: str
    status: Literal["ok", "skipped", "ineligible", "failed", "blocked", "deferred"]
    reason: str
    detail: str | None = None


# --------------------------------------------------------------------------------------
# 1. DatasetManifest
# --------------------------------------------------------------------------------------


class SampleRecord(_Contract):
    """One measured sample.  Identity levels are kept distinct on purpose."""

    sample_id: str
    donor_id: str | None = None
    biospecimen_id: str | None = None
    aliquot_id: str | None = None
    batch: str | None = None
    technical_factors: dict[str, str] = Field(default_factory=dict)
    group: str | None = None
    stimulation: str | None = None
    timepoint: str | None = None
    reference_id: str | None = None
    reference_role: ReferenceRole = ReferenceRole.NONE
    n_cells: int | None = None
    source_file: str | None = None
    source_checksum: str | None = None


class MarkerRecord(_Contract):
    """One feature of the panel, with the three independent masks kept separate."""

    feature_id: str
    channel: str | None = None
    marker: str | None = None
    panel_id: str | None = None
    role: MarkerRole = MarkerRole.OTHER
    use_for_clustering: bool = False
    use_for_correction: bool = False
    use_for_evaluation: bool = False


class ProcessingManifest(_Contract):
    """What has already been done to the values, declared rather than guessed."""

    compensated: bool | None = None
    unmixed: bool | None = None
    compensation_source: str | None = None
    qc_history: list[str] = Field(default_factory=list)
    transformation: str | None = None
    transformation_parameters: dict[str, float] = Field(default_factory=dict)
    scale: str = "unknown"
    units: str | None = None
    corrected: bool = False
    notes: str | None = None


class DatasetManifest(_Contract):
    """Immutable description of the preprocessed baseline used by every candidate."""

    dataset_id: str
    n_cells: int
    n_features: int
    expression_layer: str
    samples: list[SampleRecord]
    markers: list[MarkerRecord]
    processing: ProcessingManifest
    cell_index_columns: list[str] = Field(
        default_factory=lambda: ["cell_id", "sample_id", "source_checksum", "original_event_index"]
    )
    baseline_hash: str | None = None
    provenance: list[ProvenanceRef] = Field(default_factory=list)

    def marker_ids(self, mask: str | None = None) -> list[str]:
        """Feature ids, optionally restricted to one of the three declared masks."""
        attr = {"clustering": "use_for_clustering",
                "correction": "use_for_correction",
                "evaluation": "use_for_evaluation"}
        if mask is None:
            return [m.feature_id for m in self.markers]
        if mask not in attr:
            raise KeyError(f"unknown marker mask {mask!r}; expected one of {sorted(attr)}")
        return [m.feature_id for m in self.markers if getattr(m, attr[mask])]


# --------------------------------------------------------------------------------------
# 2. DesignAudit
# --------------------------------------------------------------------------------------


class ComparisonAssessment(_Contract):
    """One batch-pair comparison, optionally within one evaluation population."""

    comparison_id: str
    batch_a: str
    batch_b: str
    population: str | None = None
    support: ComparisonSupport
    reason: str
    n_samples_a: int = 0
    n_samples_b: int = 0
    n_shared_reference_ids: int = 0
    shared_reference_ids: list[str] = Field(default_factory=list)
    overlapping_biology: list[str] = Field(default_factory=list)


class MethodEligibility(_Contract):
    """Whether one method or detector may run, with a machine-readable reason."""

    name: str
    kind: Literal["detector", "backend", "metric"]
    eligible: bool
    reason_code: str
    reason: str
    requirements_met: dict[str, bool] = Field(default_factory=dict)
    override_allowed: bool = True


class DesignAudit(_Contract):
    """Sample-level design assessment; the gate in front of every recommendation."""

    dataset_id: str
    correction_factor: str
    nested_within: str | None = None
    protected_factors: list[str] = Field(default_factory=list)
    batches: list[str] = Field(default_factory=list)
    n_batches: int = 0
    #: protected factor -> batch -> level -> number of samples
    contingency: dict[str, dict[str, dict[str, int]]] = Field(default_factory=dict)
    design_rank: int | None = None
    design_full_rank: bool | None = None
    confounded_factors: list[str] = Field(default_factory=list)
    residual_replication: dict[str, int] = Field(default_factory=dict)
    reference_coverage: dict[str, list[str]] = Field(default_factory=dict)
    reference_graph_components: list[list[str]] = Field(default_factory=list)
    below_design_floor: bool = False
    design_floor_reason: str | None = None
    design_floor_origin: str | None = None
    comparisons: list[ComparisonAssessment] = Field(default_factory=list)
    eligibility: list[MethodEligibility] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: list[ProvenanceRef] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# 3. SamplingManifest
# --------------------------------------------------------------------------------------


class RepresentationSettings(_Contract):
    """The four pinned settings recorded on every quantitative metric record."""

    space: Literal["marker"] = "marker"
    marker_mask: str = "evaluation"
    scaling: Literal["none", "zscore_reference"] = "none"
    distance: Literal["euclidean"] = "euclidean"
    k: int = 30


class SamplingManifest(_Contract):
    """Frozen evaluation cells, labels and comparison pairs shared by all candidates."""

    dataset_id: str
    seed: int
    unit: str
    max_cells_per_sample: int | None
    label_source: LabelTier
    label_confidence: str
    label_column: str
    labels_by_cell: dict[str, str] = Field(default_factory=dict)
    selected_cell_ids: list[str] = Field(default_factory=list)
    per_sample_counts: dict[str, int] = Field(default_factory=dict)
    per_sample_original_counts: dict[str, int] = Field(default_factory=dict)
    population_counts: dict[str, int] = Field(default_factory=dict)
    population_coverage_warnings: list[str] = Field(default_factory=list)
    comparison_pairs: list[tuple[str, str]] = Field(default_factory=list)
    evaluation_markers: list[str] = Field(default_factory=list)
    representation: RepresentationSettings = Field(default_factory=RepresentationSettings)
    aggregation_weights: dict[str, float] = Field(default_factory=dict)
    splits: dict[str, list[str]] = Field(default_factory=dict)
    #: The shares the split actually used, normalised, and where they were declared.
    #: ``None`` on both means the manifest predates the field or the split was built
    #: without one; a manifest written by this build always populates both, so a
    #: reader can disagree with the proportions explicitly rather than inherit them.
    split_proportions: dict[str, float] | None = None
    split_proportion_origin: str | None = None
    provenance: list[ProvenanceRef] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# 4. DetectorResult
# --------------------------------------------------------------------------------------


class DetectorAssessment(_Contract):
    """One statement a detector makes.  Never a calibrated probability."""

    assessment_id: str
    detector: str
    namespace: Literal["native", "common"]
    scope: str
    marker: str | None = None
    batch: str | None = None
    sample_id: str | None = None
    population: str | None = None
    statistic: str
    value: float | None = None
    flag: bool | None = None
    threshold: float | None = None
    status: MetricStatus = MetricStatus.OK
    reason: str | None = None
    evidence_class: Literal[
        "technical_reference", "conditional_biological", "descriptive"
    ] = "descriptive"


class DetectorResult(_Contract):
    """Everything one detector produced, including why parts of it did not run."""

    detector: str
    detector_version: str | None = None
    namespace: Literal["native", "common"]
    run_id: str
    status: Literal["ok", "ineligible", "failed", "not_available"]
    reason: str | None = None
    eligibility: MethodEligibility | None = None
    assessments: list[DetectorAssessment] = Field(default_factory=list)
    tables: dict[str, str] = Field(default_factory=dict)
    plots: dict[str, str] = Field(default_factory=dict)
    requirements: dict[str, Any] = Field(default_factory=dict)
    parameters: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    software: list[SoftwareVersion] = Field(default_factory=list)
    provenance: list[ProvenanceRef] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# 5. CorrectionResult
# --------------------------------------------------------------------------------------


class CandidateSpecRecord(_Contract):
    """Canonical, hashable identity of one correction candidate."""

    candidate_id: str
    backend: str
    backend_version: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    seed: int | None = None
    correction_scope: dict[str, list[str]] = Field(default_factory=dict)
    input_hashes: dict[str, str] = Field(default_factory=dict)
    canonical_signature: str | None = None


class CorrectionResult(_Contract):
    """Outcome of one candidate, with success/ineligible/failed kept distinct."""

    candidate: CandidateSpecRecord
    status: CorrectionStatus
    reason: str | None = None
    output_path: str | None = None
    output_hash: str | None = None
    output_scale: str | None = None
    n_cells: int | None = None
    markers_changed: list[str] = Field(default_factory=list)
    markers_unchanged: list[str] = Field(default_factory=list)
    samples_changed: list[str] = Field(default_factory=list)
    populations_changed: list[str] = Field(default_factory=list)
    magnitude_by_marker: dict[str, float] = Field(default_factory=dict)
    magnitude_by_sample: dict[str, float] = Field(default_factory=dict)
    magnitude_by_population: dict[str, float] = Field(default_factory=dict)
    validation: dict[str, bool] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    software: list[SoftwareVersion] = Field(default_factory=list)
    resources: dict[str, float] = Field(default_factory=dict)
    provenance: list[ProvenanceRef] = Field(default_factory=list)
    #: Eligibility constraint(s) waived by an exploratory override, with the reason
    #: code and the configuration flag that allowed the waiver.  Empty unless the run
    #: overrode an ineligibility; a waived constraint is always visible here.
    waived_eligibility: list[dict[str, str]] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# 6. MetricResult / MetricRecord
# --------------------------------------------------------------------------------------


class MetricRecord(_Contract):
    """One metric value.  The fields below are mandatory by the cross-stage contract."""

    run_id: str
    candidate_id: str
    metric_namespace: Literal["native", "common", "compatibility"]
    metric_name: str
    metric_version: str
    stage: str
    marker: str | None = None
    population: str | None = None
    comparison_id: str | None = None
    value: float | None = None
    direction: MetricDirection = MetricDirection.NEUTRAL
    scale: str = "unitless"
    status: MetricStatus = MetricStatus.OK
    reason: str | None = None
    n_samples: int | None = None
    n_donors: int | None = None
    n_cells: int | None = None
    coverage: float | None = None
    aggregation_weight: float = 1.0
    uncertainty_low: float | None = None
    uncertainty_high: float | None = None
    uncertainty_kind: str | None = None
    label_tier: LabelTier | None = None
    label_confidence: str | None = None
    representation: RepresentationSettings | None = None
    provenance: list[ProvenanceRef] = Field(default_factory=list)


class MetricResult(_Contract):
    """All metric records produced for one candidate by one metric group."""

    run_id: str
    candidate_id: str
    metric_group: str
    records: list[MetricRecord] = Field(default_factory=list)
    status: Literal["ok", "partial", "failed", "ineligible"] = "ok"
    reason: str | None = None
    warnings: list[str] = Field(default_factory=list)
    software: list[SoftwareVersion] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# 7. DecisionResult
# --------------------------------------------------------------------------------------


class CandidateAssessment(_Contract):
    """Per-candidate verdict with the evidence that produced it."""

    candidate_id: str
    state: DecisionState
    reason: str
    guardrails_failed: list[str] = Field(default_factory=list)
    evidence: dict[str, float | None] = Field(default_factory=dict)
    coverage: dict[str, float | None] = Field(default_factory=dict)
    comparable_evidence_key: str | None = None


class DecisionResult(_Contract):
    """The recommendation layer.  Abstention is a first-class outcome."""

    run_id: str
    policy: str
    policy_settings: dict[str, Any] = Field(default_factory=dict)
    threshold_origins: dict[str, str] = Field(default_factory=dict)
    assessments: list[CandidateAssessment] = Field(default_factory=list)
    recommended_candidate_id: str | None = None
    recommendation_withheld_reason: str | None = None
    pareto_front: list[str] = Field(default_factory=list)
    pareto_comparable: bool = False
    overrides: list[dict[str, str]] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# Schema export for the R boundary
# --------------------------------------------------------------------------------------


def all_contracts() -> dict[str, type[_Contract]]:
    """The seven contracts named in the cross-stage section of the specification."""

    return {
        "DatasetManifest": DatasetManifest,
        "DesignAudit": DesignAudit,
        "SamplingManifest": SamplingManifest,
        "DetectorResult": DetectorResult,
        "CorrectionResult": CorrectionResult,
        "MetricResult": MetricResult,
        "DecisionResult": DecisionResult,
    }


def export_json_schemas(out_dir: str | Path) -> dict[str, Path]:
    """Write one JSON Schema per contract, for validation from R or any other reader."""

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, model in all_contracts().items():
        path = out / f"{name}.schema.json"
        path.write_text(json.dumps(model.model_json_schema(), indent=2) + "\n", encoding="utf-8")
        written[name] = path
    return written
