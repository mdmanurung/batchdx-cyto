"""The detector interface.

A detector returns assessments, plots, provenance, requirements and warnings.  It
distinguishes descriptive batch association from technical-reference evidence, and it
never converts a count of flags into a probability: correlated metrics firing together
are one observation seen three ways, not three independent confirmations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.schemas.contracts import DesignAudit, DetectorResult, SamplingManifest


@dataclass
class DetectorContext:
    """Shared inputs for every detector."""

    run_id: str
    design: DesignAudit
    sampling: SamplingManifest | None
    out_dir: str
    seed: int = 473
    rscript: str = "Rscript"
    upstream_dir: str | None = None
    #: Fixed positivity thresholds in the declared coordinate system (Stage 4).
    threshold_path: str | None = None
    #: |Spearman rho| above which a within-file trend is reported as monotonic drift.
    #: Engineering default, not a validated cutoff; reported with every assessment.
    drift_rho_threshold: float = 0.2
    #: Relative step between consecutive block event rates (over the median block rate)
    #: above which an event-rate discontinuity is reported.  Engineering default, not
    #: a validated cutoff; reported with every assessment.
    drift_rate_threshold: float = 3.0
    #: Nominal level for the composition test.  Reported, never used to declare a
    #: technical cause on its own.
    composition_alpha: float = 0.05


class Detector(Protocol):
    """Protocol implemented by every detector."""

    name: str
    namespace: str
    eligibility_key: str

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        """Produce a :class:`DetectorResult`; never raise for a scientific skip."""


_DETECTORS: dict[str, Detector] = {}


def register_detector(detector: Detector) -> Detector:
    _DETECTORS[detector.name] = detector
    return detector


def get_detector(name: str) -> Detector:
    if name not in _DETECTORS:
        raise KeyError(f"unknown detector {name!r}; available: {sorted(_DETECTORS)}")
    return _DETECTORS[name]


def available_detectors() -> dict[str, Detector]:
    return dict(_DETECTORS)


def eligibility_for(design: DesignAudit, key: str):
    """Look up the eligibility entry the design audit produced for a detector."""

    for entry in design.eligibility:
        if entry.name == key:
            return entry
    return None
