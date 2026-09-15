"""Strict configuration schema (``schema_version: 2``).

Unknown keys are rejected everywhere.  Backend parameters left unspecified are
resolved to the pinned documented defaults recorded in
:mod:`batchdx.cyto.adapters` and serialised into the resolved configuration, so a
run never depends on an upstream default silently changing.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from batchdx.cyto._version import CONFIG_SCHEMA_VERSION


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class ProjectConfig(_Strict):
    id: str
    seed: int = 473


class InputConfig(_Strict):
    format: Literal["h5ad", "fcs"]
    path: str
    expression_layer: str
    samples: str
    markers: str
    processing_manifest: str
    #: For ``format: fcs`` the path is a directory and this pattern selects files.
    fcs_glob: str = "*.fcs"


class PreprocessingConfig(_Strict):
    enabled: bool = False
    gating: str | None = None
    transformation: Literal["asinh", "logicle", "none"] | None = None
    cofactor_by_panel: dict[str, float] = Field(default_factory=dict)
    peacoqc: bool = False
    time_channel: str | None = None

    @model_validator(mode="after")
    def _require_declared_cofactor(self) -> PreprocessingConfig:
        if self.enabled and self.transformation == "asinh" and not self.cofactor_by_panel:
            raise ValueError(
                "preprocessing.transformation: 'asinh' requires "
                "preprocessing.cofactor_by_panel; there is no universal cofactor"
            )
        return self


class DesignConfig(_Strict):
    correction_batch: str
    nested_within: str | None = None
    protect: list[str] = Field(default_factory=list)
    donor: str | None = None
    reference_identity: str | None = None
    reference_role: str | None = None
    allow_exploratory_override: bool = False
    #: Design floor.  Both values are configuration with a declared origin, not
    #: hidden constants; see ``design_floor_origin`` in the emitted DesignAudit.
    minimum_batches: int = 3
    require_shared_reference: bool = True
    minimum_batches_origin: str = (
        "package default; declared in config/ and reported in design_audit.json"
    )


class RepresentationConfig(_Strict):
    space: Literal["marker"] = "marker"
    marker_mask: Literal["evaluation", "clustering", "correction"] = "evaluation"
    scaling: Literal["none", "zscore_reference"] = "none"
    distance: Literal["euclidean"] = "euclidean"
    k: int = 30


class SamplingConfig(_Strict):
    unit: str = "sample_id"
    max_cells_per_sample: int | None = 5000
    min_cells_per_population: int = 50


class EvaluationConfig(_Strict):
    label_source: Literal["annotation", "gating", "clustering"]
    labels: str | None = None
    gating_thresholds: str | None = None
    clustering_k: int = 8
    representation: RepresentationConfig = Field(default_factory=RepresentationConfig)
    sampling: SamplingConfig = Field(default_factory=SamplingConfig)
    partitions: Literal["frozen_common"] = "frozen_common"
    uncertainty_unit: str | None = None
    policy: Literal["conservative", "permissive"] = "conservative"

    @model_validator(mode="after")
    def _require_label_inputs(self) -> EvaluationConfig:
        if self.label_source == "annotation" and not self.labels:
            raise ValueError("evaluation.labels is required when label_source is 'annotation'")
        if self.label_source == "gating" and not self.gating_thresholds:
            raise ValueError(
                "evaluation.gating_thresholds is required when label_source is 'gating'"
            )
        return self


class CandidateConfig(_Strict):
    id: str
    backend: str
    parameters: dict[str, Any] = Field(default_factory=dict)


class CorrectionConfig(_Strict):
    control_policy: Literal["measured_only", "allow_pseudo_controls"] = "measured_only"
    candidates: list[CandidateConfig]

    @field_validator("candidates")
    @classmethod
    def _unique_ids(cls, v: list[CandidateConfig]) -> list[CandidateConfig]:
        seen: set[str] = set()
        for c in v:
            if c.id in seen:
                raise ValueError(f"duplicate candidate id {c.id!r}")
            seen.add(c.id)
        if not any(c.backend == "none" for c in v):
            raise ValueError(
                "correction.candidates must include an identity candidate "
                "(backend: none); 'no correction' is a valid result and is always benchmarked"
            )
        return v


class ReportConfig(_Strict):
    offline: bool = True
    redact_identifiers: bool = True
    title: str = "batchdx-cyto report"


class ExportsConfig(_Strict):
    h5ad: bool = True
    fcs: bool = False


class DetectorsConfig(_Strict):
    """Which detectors to attempt.  Availability is resolved at run time."""

    acquisition_drift: bool = True
    composition: bool = True
    cycombine_express: bool = False
    cycombine_cluster: bool = False
    cytobatchflagr: bool = False


class UpstreamConfig(_Strict):
    """Where pinned upstream checkouts live, for adapters that shell out to R."""

    lock: str | None = None
    checkout_dir: str | None = None
    rscript: str = "Rscript"


class Config(_Strict):
    schema_version: int
    project: ProjectConfig
    input: InputConfig
    preprocessing: PreprocessingConfig = Field(default_factory=PreprocessingConfig)
    design: DesignConfig
    evaluation: EvaluationConfig
    correction: CorrectionConfig
    detectors: DetectorsConfig = Field(default_factory=DetectorsConfig)
    upstream: UpstreamConfig = Field(default_factory=UpstreamConfig)
    report: ReportConfig = Field(default_factory=ReportConfig)
    exports: ExportsConfig = Field(default_factory=ExportsConfig)

    @field_validator("schema_version")
    @classmethod
    def _check_version(cls, v: int) -> int:
        if v != CONFIG_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version {v} is not supported; this build reads "
                f"schema_version {CONFIG_SCHEMA_VERSION}"
            )
        return v

    # -- analytical identity ---------------------------------------------------------

    def analytical_fingerprint(self) -> str:
        """Hash of the scientifically material settings only.

        Cosmetic reporting settings are excluded so that changing a report title does
        not invalidate correction outputs.
        """
        payload = self.model_dump(mode="json", exclude={"report"})
        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()

    def resolved_dict(self) -> dict[str, Any]:
        return self.model_dump(mode="json")


class ConfigError(ValueError):
    """Raised for an invalid configuration, with the offending keys named."""


def load_config(path: str | Path) -> Config:
    """Read and strictly validate a configuration file."""

    p = Path(path)
    if not p.exists():
        raise ConfigError(f"configuration file not found: {p}")
    try:
        raw = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:  # pragma: no cover - passthrough
        raise ConfigError(f"{p}: not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{p}: top level of the configuration must be a mapping")
    return parse_config(raw, source=str(p))


def parse_config(raw: dict[str, Any], source: str = "<dict>") -> Config:
    from pydantic import ValidationError

    try:
        return Config.model_validate(raw)
    except ValidationError as exc:
        lines = [f"{source}: configuration is invalid"]
        for err in exc.errors():
            loc = ".".join(str(x) for x in err["loc"])
            lines.append(f"  - {loc}: {err['msg']}")
        raise ConfigError("\n".join(lines)) from exc
