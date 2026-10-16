"""The correction backend interface and registry.

Every backend declares, up front, what metadata it needs, which reference layouts it
supports, whether it can separate fitting from application, and which parameters it
accepts.  Unsupported arguments are rejected rather than ignored, and a runtime
failure never returns the unchanged input as a successful correction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np

from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.io.hashing import canonical_hash
from batchdx.cyto.schemas.contracts import (
    CandidateSpecRecord,
    CorrectionResult,
    DesignAudit,
    SamplingManifest,
    SoftwareVersion,
)
from batchdx.cyto.schemas.enums import CorrectionStatus


class BackendError(RuntimeError):
    """A correction backend failed or refused to run."""


class BackendNotAvailable(BackendError):
    """The backend is absent or unusable in this environment, rather than broken.

    This is a distinct type, not a message convention, because
    ``docs/decisions/0006-candidate-results-and-failure-handling.md`` requires
    ``not_available`` and ``failed`` to stay distinguishable.  The two mean
    different things to a reader of a report: "this method was not installed here"
    is an environment fact that resolves by installing it, while "this method ran
    and broke" is a statement about the method.  A backend that raises
    ``BackendError`` is recorded as ``failed``; only this subclass is recorded as
    ``not_available``.

    Substring matching on the message cannot make that distinction: any real
    failure whose text happens to mention a missing package would be reported as
    a missing package, which is how a genuine defect gets recorded as an
    environment gap.
    """


class BackendIneligible(BackendError):
    """The backend is installed and works, but this design cannot support a fit.

    A design-coverage refusal is a statement about the DATA, not about the
    installation: "no measured reference falls in the training split" and "the frozen
    partition is not arm-balanced" are properties of the dataset and the split, and
    they resolve by declaring more biospecimens or re-partitioning -- never by
    installing a package.  A missing R package is :class:`BackendNotAvailable`, an
    environment fact that resolves by installing it.

    The two are distinct exception types, not message conventions, because substring
    matching on the message cannot make the distinction: any real failure whose text
    happens to mention a training split would be recorded as a coverage gap, which is
    how a genuine defect gets recorded as a data limitation.
    """


@dataclass(frozen=True)
class BackendCapabilities:
    """What a backend needs and what it promises."""

    name: str
    required_sample_fields: tuple[str, ...] = ()
    reference_layouts: tuple[str, ...] = ()
    supports_covariates: bool = False
    supports_marker_mask: bool = False
    input_scale: str = "transformed"
    output_scale: str = "transformed"
    supports_fit_apply: bool = False
    fit_apply_verified: bool = False
    resource_profile: str = "memory: all relevant cells held jointly"
    #: Parameter name -> pinned default.  Defaults are serialised into the resolved
    #: configuration so a run does not silently inherit an upstream change.
    accepted_parameters: dict[str, Any] = field(default_factory=dict)
    #: Parameter name -> why it is refused.
    rejected_parameters: dict[str, str] = field(default_factory=dict)
    notes: str = ""


@dataclass
class BackendContext:
    """Everything a backend may read, beyond the dataset itself."""

    design: DesignAudit
    sampling: SamplingManifest | None
    work_dir: str
    seed: int
    rscript: str = "Rscript"
    upstream_dir: str | None = None
    threads: int = 1


class CorrectionBackend(Protocol):
    """Protocol implemented by every correction backend."""

    capabilities: BackendCapabilities

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        """Resolve the software this backend needs, without assuming it is present."""

    def run(
        self,
        dataset: Dataset,
        candidate_id: str,
        parameters: dict[str, Any],
        context: BackendContext,
    ) -> tuple[np.ndarray, CorrectionResult]:
        """Correct and return ``(matrix, result)``; raise :class:`BackendError` on failure."""


_REGISTRY: dict[str, CorrectionBackend] = {}


def register_backend(backend: CorrectionBackend) -> CorrectionBackend:
    _REGISTRY[backend.capabilities.name] = backend
    return backend


def get_backend(name: str) -> CorrectionBackend:
    if name not in _REGISTRY:
        raise BackendError(
            f"unknown correction backend {name!r}; available backends: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[name]


def available_backends() -> dict[str, BackendCapabilities]:
    return {name: b.capabilities for name, b in _REGISTRY.items()}


def resolve_parameters(
    capabilities: BackendCapabilities, requested: dict[str, Any]
) -> dict[str, Any]:
    """Merge requested parameters onto the pinned defaults, rejecting anything unsupported."""

    refused = [k for k in requested if k in capabilities.rejected_parameters]
    if refused:
        details = "; ".join(f"{k}: {capabilities.rejected_parameters[k]}" for k in refused)
        raise BackendError(f"backend {capabilities.name!r} refuses parameter(s) {refused}. {details}")
    unknown = [
        k for k in requested
        if k not in capabilities.accepted_parameters and k not in capabilities.rejected_parameters
    ]
    if unknown:
        raise BackendError(
            f"backend {capabilities.name!r} does not accept parameter(s) {unknown}; "
            f"accepted parameters are {sorted(capabilities.accepted_parameters)}"
        )
    resolved = dict(capabilities.accepted_parameters)
    resolved.update(requested)
    return resolved


def candidate_signature(
    candidate_id: str,
    backend: str,
    backend_version: str | None,
    parameters: dict[str, Any],
    seed: int | None,
    input_hashes: dict[str, str],
) -> str:
    """Stable identity for a candidate under a canonical parameter serialisation."""

    return canonical_hash({
        "candidate_id": candidate_id,
        "backend": backend,
        "backend_version": backend_version,
        "parameters": parameters,
        "seed": seed,
        "input_hashes": input_hashes,
    })


def make_spec(
    candidate_id: str,
    backend: str,
    backend_version: str | None,
    parameters: dict[str, Any],
    seed: int | None,
    dataset: Dataset,
    correction_markers: list[str],
) -> CandidateSpecRecord:
    input_hashes = {"baseline": dataset.baseline_hash()}
    return CandidateSpecRecord(
        candidate_id=candidate_id,
        backend=backend,
        backend_version=backend_version,
        parameters=parameters,
        seed=seed,
        correction_scope={"markers": list(correction_markers)},
        input_hashes=input_hashes,
        canonical_signature=candidate_signature(
            candidate_id, backend, backend_version, parameters, seed, input_hashes
        ),
    )


def failure_result(
    spec: CandidateSpecRecord,
    status: CorrectionStatus,
    reason: str,
    software: list[SoftwareVersion] | None = None,
) -> CorrectionResult:
    """A failed or ineligible candidate, with no output path and no fabricated values."""

    return CorrectionResult(
        candidate=spec, status=status, reason=reason, software=software or [],
        validation={"produced_output": False},
    )
