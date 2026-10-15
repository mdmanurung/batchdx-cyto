"""Strict configuration schema (``schema_version: 2``).

Unknown keys are rejected everywhere.  So are keys this build cannot honour, in two
distinct ways, because a configuration that is accepted and then ignored is the worst
outcome available: the run looks configured and is not.

* :data:`UNSUPPORTED_KEYS` names keys this build removed entirely.  Setting one is
  refused with the reason and the thing to use instead.
* :data:`SUPPORTED_VALUES` names keys that stay, but only for the values the package
  can actually deliver.  Any other value is refused by name rather than accepted and
  dropped on the floor.  Where a key has exactly one meaningful value the schema pins
  it, so the key records an invariant instead of advertising a choice that does not
  exist.

``tests/unit/test_config_keys.py`` walks every model below and fails the build if a leaf
is neither read by a non-schema module nor listed, with a justification, in the test's
``DECLARATIVE_ONLY`` allow-list.  That test is what stops this drifting back.

Backend parameters left unspecified are resolved to the pinned documented defaults
recorded in :mod:`batchdx.cyto.adapters` and serialised into the resolved configuration,
so a run never depends on an upstream default silently changing.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from batchdx.cyto._version import CONFIG_SCHEMA_VERSION

#: Keys this build does not support.  Each was once accepted by the schema and then
#: never read, so a user who set one got a run that looked configured and was not.  The
#: value is the reason, and is shown in the validation error.
#:
#: ``preprocessing.gating`` and ``preprocessing.peacoqc`` were in this table until
#: Tranche 3 gave them readers: they are now paths, read by
#: ``batchdx.cyto.preprocess.gating`` and ``batchdx.cyto.preprocess.peacoqc``, and they
#: left it on the same change that made them mean something.  Their names were kept
#: rather than renamed, so a user who still has the old value gets a type error naming
#: the key instead of the vaguer "unknown key" a rename would have produced.
UNSUPPORTED_KEYS: dict[str, str] = {
    "preprocessing.time_channel": (
        "the acquisition-drift detector reads the acquisition order recorded at "
        "ingestion and falls back to the original event index; it never reads a channel "
        "named by configuration, so declaring one here would change nothing. Record the "
        "time channel in the marker manifest instead"
    ),
}

#: Dotted key path -> (the values this build honours, why every other value is refused).
#: A key is kept, rather than removed, when a configuration that states it should agree
#: with one that omits it; what it must not do is accept a second value and not act on
#: it.
SUPPORTED_VALUES: dict[str, tuple[tuple[Any, ...], str]] = {
    "design.reference_identity": (
        ("reference_id", None),
        "the sample manifest reader fixes the reference-identity column at "
        "'reference_id' (batchdx.cyto.io.manifests), and the design audit reads "
        "SampleRecord.reference_id; naming a different column would be accepted and "
        "ignored. Declare the reference in the 'reference_id' column of samples.tsv",
    ),
    "design.reference_role": (
        ("reference_role", None),
        "the sample manifest reader fixes the reference-role column at "
        "'reference_role' (batchdx.cyto.io.manifests); naming a different column would "
        "be accepted and ignored. Declare the role in the 'reference_role' column of "
        "samples.tsv, and note that 'biological_pseudo_control' is never counted as "
        "measured technical evidence whatever this key says",
    ),
    "evaluation.representation.scaling": (
        ("none",),
        "this build applies no scaling to the evaluation representation: every "
        "common-namespace metric is computed per marker on the declared coordinate "
        "system, and nothing in the package or in batchdx.cyto.metrics reads a joint "
        "space in which markers could be standardised. The key is still recorded on "
        "every metric record, so a future standardising metric has a declared setting "
        "to honour; 'zscore_reference' additionally names a reference population that "
        "no configuration key selects",
    ),
    "evaluation.sampling.unit": (
        ("sample_id",),
        "evaluation cells are drawn per sample: the frozen sampling manifest records "
        "per-sample counts, aggregation weights and the cap, and every downstream "
        "metric aggregates per sample, so a donor- or biospecimen-level draw could not "
        "be recorded truthfully. Pool the samples of one donor upstream, or draw per "
        "sample and report the donor structure through design.donor",
    ),
    "evaluation.uncertainty_unit": (
        (None, "sample_id", "donor_id", "biospecimen_id", "aliquot_id"),
        "this build declares it but computes no unit-level uncertainty: "
        "uncertainty_low/uncertainty_high/uncertainty_kind are null on every metric "
        "record, because donor-level uncertainty belongs to Stage 8 biological "
        "preservation, which is not implemented (docs/blocked.md). The value is "
        "surfaced as a coverage warning in the sampling manifest and in the report",
    ),
    "correction.control_policy": (
        ("measured_only",),
        "'allow_pseudo_controls' would have to make a biological_pseudo_control count as "
        "technical evidence, which docs/usage.md and the reference-role contract forbid "
        "whatever the configuration says; the design audit and the detectors read the "
        "role from the sample manifest and never consult this key. Only measured "
        "technical replicates and anchors are technical evidence in this build",
    ),
    "upstream.lock": (
        (None, "config/upstream_lock.yaml"),
        "the pinned manifest is read by the Stage 0 audit (scripts/audit_upstream.py, "
        "tests/integration/test_stage0_audit.py) and by scripts/fetch_upstream.sh, out "
        "of band; no stage reads it while a run is executing, so pointing it elsewhere "
        "would not change a run. Keep the one path so the run and the audit cannot "
        "disagree about what is pinned",
    ),
    "report.offline": (
        (True,),
        "the only renderer in this build is one self-contained offline HTML report; "
        "there is no online, served or partially-remote output path to switch on, so "
        "'false' would name a mode that does not exist",
    ),
    "exports.h5ad": (
        (True,),
        "every successful candidate writes corrected.h5ad, because the metrics rule "
        "reads that file and would otherwise have nothing to score; 'false' would leave "
        "the run unable to proceed",
    ),
    # exports.fcs is deliberately absent: the report stage calls the exporter when the
    # key is true, so both values are supported and none is pinned.
}


def _check_supported(dotted: str, value: Any) -> None:
    """Refuse a value this build cannot deliver, naming the key and the reason."""

    allowed, why = SUPPORTED_VALUES[dotted]
    if value in allowed:
        return
    accepted = ", ".join("unset" if v is None else repr(v) for v in allowed)
    raise ValueError(
        f"{dotted}: {value!r} is not supported by this build, because {why}. "
        f"Accepted: {accepted}."
    )


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


#: The one default seed for the whole package.  A run's seed is a single source:
#: :attr:`ProjectConfig.seed`, and every other default that used to repeat the literal
#: ``473`` is derived from this name — the synthetic generator
#: (:mod:`batchdx.cyto.simulate`), the shipped examples
#: (:mod:`batchdx.cyto.examples`) and the detector context
#: (:class:`batchdx.cyto.diagnostics.base.DetectorContext`).
#:
#: It lives here, in the lightest module in the package, and not in
#: ``batchdx.cyto.simulate``: that module imports anndata, which costs ~2.9 s, and a
#: schema module must not drag a data generator in to name a constant.
#:
#: The seed is material, not cosmetic: ``Config.analytical_fingerprint`` excludes only
#: ``report``, so two runs that differ only in ``project.seed`` are correctly held to be
#: different analyses.  ``tests/scientific/test_multiseed.py`` asserts that.
DEFAULT_SEED = 473


class ProjectConfig(_Strict):
    id: str
    seed: int = DEFAULT_SEED


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
    """A declaration about the data, and — since Tranche 3 — a declaration of the
    upstream step that produced it.

    ``enabled`` and ``transformation`` are read exactly once, to refuse a configuration
    that claims a transformation the processing manifest contradicts, and
    ``cofactor_by_panel`` is read only to refuse a universal default.  The *run* still
    never transforms events: it reads an already-preprocessed baseline, and
    ``_validate_processing`` in ``io.dataset`` refuses a configuration that asks it to.

    ``gating`` and ``peacoqc`` are **paths to the upstream step's own declaration**, not
    instructions to this run.  They are read by
    :func:`batchdx.cyto.preprocess.gating.gates_from_config` and
    :func:`batchdx.cyto.preprocess.peacoqc.settings_from_config`, which apply the gate
    table and run PeacoQC over the raw events *before* a baseline exists.  A run that
    never calls those functions ignores them, exactly as it ignores
    ``cofactor_by_panel``; that is stated here so the key is not mistaken for a switch
    the correction workflow reads.
    """

    enabled: bool = False
    # "logicle" was removed from this Literal on 2026-09-30.  It was a declaration, not a
    # capability: this build implements no logicle, io/export.py has no tested inverse
    # for one, and the pinned reference implementation contains none either
    # (`grep -i logicle` over .upstream/spiscy returns zero matches).  Accepting it
    # validated a value and bought nothing, which is exactly the accepted-and-ignored
    # failure this module exists to prevent.  preprocess.transform.LOGICLE keeps the
    # refusal reason reachable in code.
    transformation: Literal["asinh", "none"] | None = None
    #: Declared per-panel cofactors.  Never applied: the cofactor that the FCS
    #: exporter inverts is the one in the processing manifest, which is what a run
    #: actually uses.
    cofactor_by_panel: dict[str, float] = Field(default_factory=dict)
    #: Path to a two-dimensional gate table in the reference implementation's YAML
    #: shape.  ``None`` means no gate is declared, which is a statement about the
    #: upstream step, not an instruction to skip one that exists.
    gating: str | None = None
    #: Path to a PeacoQC settings file.  ``None`` means no PeacoQC pass is declared.
    peacoqc: str | None = None

    @field_validator("gating", "peacoqc")
    @classmethod
    def _a_declared_path_is_a_path(cls, v: str | None) -> str | None:
        """Refuse a blank path rather than accept one that resolves to nothing.

        ``gating: ""`` would otherwise validate and then name a directory. An empty
        string is not "no gate"; it is a path to the current working directory, and
        reading a gate table from there is not what anyone meant.
        """

        if v is not None and not str(v).strip():
            raise ValueError(
                "this key names a file. An empty string is a path to the current working "
                "directory, not an absence: leave the key out entirely to declare that "
                "there is no such step"
            )
        return v

    @model_validator(mode="after")
    def _require_declared_cofactor(self) -> PreprocessingConfig:
        if self.enabled and self.transformation == "asinh" and not self.cofactor_by_panel:
            raise ValueError(
                "preprocessing.transformation: 'asinh' requires "
                "preprocessing.cofactor_by_panel; there is no universal cofactor"
            )
        return self


class BiologicalContrastConfig(_Strict):
    """One prespecified biological contrast that correction must not destroy.

    A contrast is a difference between two levels of a sample-level factor (usually
    the biological ``group``) on one declared marker.  The tolerance is the maximum
    change in that difference a candidate may cause for the contrast to count as
    preserved.  Both the tolerance and *where it came from* are configuration: a
    tolerance without a declared origin is rejected, because an unexplained cutoff is
    not a predeclared tolerance.
    """

    contrast_id: str
    marker: str
    group_a: str
    group_b: str
    #: Sample-level factor whose levels ``group_a``/``group_b`` are.  Must resolve on
    #: the sample manifest; ``group`` is the usual choice.
    group_factor: str = "group"
    #: Maximum tolerated change in the contrast's effect size (in units of the declared
    #: expression scale), after correction relative to before.
    preservation_tolerance: float
    #: Declared origin of ``preservation_tolerance`` (e.g. a preanalysis document, a
    #: published minimal clinically important difference, or an explicit pilot).
    preservation_tolerance_origin: str

    @field_validator("preservation_tolerance")
    @classmethod
    def _positive_tolerance(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("preservation_tolerance must be strictly positive")
        return v

    @field_validator("preservation_tolerance_origin")
    @classmethod
    def _origin_declared(cls, v: str) -> str:
        if not v.strip():
            raise ValueError(
                "preservation_tolerance_origin must name where the tolerance comes from; "
                "a threshold without a declared origin is not predeclared"
            )
        return v


class DesignConfig(_Strict):
    correction_batch: str
    nested_within: str | None = None
    protect: list[str] = Field(default_factory=list)
    donor: str | None = None
    #: Names of the sample-manifest columns that carry reference identity and role.
    #: Pinned, not free: see ``SUPPORTED_VALUES``.
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
    #: Prespecified biological contrasts that every candidate is assessed against.
    #: Empty by default: no contrast is invented that the analyst did not declare.
    biological_contrasts: list[BiologicalContrastConfig] = Field(default_factory=list)

    @field_validator("biological_contrasts")
    @classmethod
    def _unique_contrast_ids(
        cls, v: list[BiologicalContrastConfig]
    ) -> list[BiologicalContrastConfig]:
        seen: set[str] = set()
        for c in v:
            if c.contrast_id in seen:
                raise ValueError(f"duplicate biological contrast id {c.contrast_id!r}")
            seen.add(c.contrast_id)
        return v

    @model_validator(mode="after")
    def _check_reference_columns(self) -> DesignConfig:
        _check_supported("design.reference_identity", self.reference_identity)
        _check_supported("design.reference_role", self.reference_role)
        return self


class RepresentationConfig(_Strict):
    """The representation contract recorded on every quantitative metric record.

    ``marker_mask`` is the only setting that selects data: it chooses which markers are
    frozen into the evaluation partition.  The rest are pinned in this build — there is
    no other space, no scaling applied to the values, one distance, and no code path
    that builds a k-nearest-neighbour graph over them.
    """

    space: Literal["marker"] = "marker"
    marker_mask: Literal["evaluation", "clustering", "correction"] = "evaluation"
    scaling: Literal["none"] = "none"
    distance: Literal["euclidean"] = "euclidean"
    k: int = 30

    @field_validator("scaling", mode="before")
    @classmethod
    def _scaling_is_pinned(cls, v: Any) -> Any:
        _check_supported("evaluation.representation.scaling", v)
        return v


class SamplingConfig(_Strict):
    #: The unit the draw is balanced over.  Pinned to ``sample_id``; see
    #: ``SUPPORTED_VALUES`` for why a donor- or biospecimen-level unit cannot be
    #: honoured truthfully by the frozen manifest.
    unit: str = "sample_id"
    max_cells_per_sample: int | None = 5000
    min_cells_per_population: int = 50

    @field_validator("unit", mode="before")
    @classmethod
    def _unit_is_pinned(cls, v: Any) -> Any:
        _check_supported("evaluation.sampling.unit", v)
        return v


class EvaluationConfig(_Strict):
    label_source: Literal["annotation", "gating", "clustering"]
    labels: str | None = None
    gating_thresholds: str | None = None
    clustering_k: int = 8
    representation: RepresentationConfig = Field(default_factory=RepresentationConfig)
    sampling: SamplingConfig = Field(default_factory=SamplingConfig)
    partitions: Literal["frozen_common"] = "frozen_common"
    #: Declared, never computed.  Kept and surfaced rather than removed because the
    #: metric contract reserves ``uncertainty_kind`` for it; see
    #: ``SUPPORTED_VALUES``.
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
        if self.uncertainty_unit is not None:
            _check_supported("evaluation.uncertainty_unit", self.uncertainty_unit)
        return self


class CandidateConfig(_Strict):
    id: str
    backend: str
    parameters: dict[str, Any] = Field(default_factory=dict)

    @field_validator("id")
    @classmethod
    def _candidate_directory_name(cls, v: str) -> str:
        if not v or v in {".", ".."} or "/" in v or "\\" in v:
            raise ValueError("candidate id must be a non-empty directory name without path separators")
        return v

    @field_validator("backend")
    @classmethod
    def _known_backend(cls, v: str) -> str:
        from batchdx.cyto.adapters import available_backends

        backends = available_backends()
        if v not in backends:
            raise ValueError(
                f"unknown correction backend {v!r}; available backends: {sorted(backends)}"
            )
        return v


class CorrectionConfig(_Strict):
    control_policy: Literal["measured_only"] = "measured_only"
    candidates: list[CandidateConfig]

    @field_validator("control_policy", mode="before")
    @classmethod
    def _control_policy_is_pinned(cls, v: Any) -> Any:
        _check_supported("correction.control_policy", v)
        return v

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
    #: Pinned to true: the only renderer emits one self-contained offline report.
    offline: bool = True
    redact_identifiers: bool = True
    title: str = "batchdx-cyto report"

    @field_validator("offline", mode="before")
    @classmethod
    def _offline_is_pinned(cls, v: Any) -> Any:
        _check_supported("report.offline", v)
        return v


class ExportsConfig(_Strict):
    #: Pinned: see ``SUPPORTED_VALUES``.  The corrected H5AD is not optional.  The FCS
    #: exporter is called by the report stage when ``fcs`` is true, so both values work.
    h5ad: bool = True
    fcs: bool = False

    @model_validator(mode="after")
    def _exports_are_pinned(self) -> ExportsConfig:
        _check_supported("exports.h5ad", self.h5ad)
        return self


class DetectorsConfig(_Strict):
    """Which detectors to attempt.  Availability is resolved at run time.

    Read through ``pipeline.DETECTOR_FLAGS`` by ``getattr`` in the Snakefile and in
    ``stage_status``, so the detector roster is static and comes from the validated
    configuration rather than from whichever result folders exist.
    """

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

    @field_validator("lock", mode="before")
    @classmethod
    def _lock_is_pinned(cls, v: Any) -> Any:
        _check_supported("upstream.lock", v)
        return v


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

    @model_validator(mode="before")
    @classmethod
    def _refuse_removed_keys(cls, raw: Any) -> Any:
        """Name a key this build removed, rather than reporting it as unknown."""

        if isinstance(raw, dict):
            for dotted, why in UNSUPPORTED_KEYS.items():
                section, _, leaf = dotted.partition(".")
                block = raw.get(section)
                if isinstance(block, dict) and leaf in block:
                    raise ValueError(
                        f"{dotted}: this build does not support this key, because {why}. "
                        f"It was removed rather than left accepted-and-ignored."
                    )
        return raw

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
            msg = err["msg"]
            # A cross-field rule reports no location; the message already names the key.
            if err["loc"]:
                lines.append(f"  - {'.'.join(str(x) for x in err['loc'])}: {msg}")
            else:
                lines.append(f"  {msg.removeprefix('Value error, ')}")
        raise ConfigError("\n".join(lines)) from exc
