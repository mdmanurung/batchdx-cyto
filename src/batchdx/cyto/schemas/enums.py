"""Closed vocabularies used across contracts.

Every enum here is part of the on-disk contract.  Downstream code must branch on
these values rather than on the presence of a file or directory.
"""

from __future__ import annotations

from enum import Enum


class StrEnum(str, Enum):
    """``str`` mixin enum so values serialise as plain strings in JSON."""

    def __str__(self) -> str:  # pragma: no cover - trivial
        return str(self.value)


class LabelTier(StrEnum):
    """Provenance tier of the evaluation labels (Stage 4).

    The tier is carried on every metric record so that a reader can tell whether a
    score rests on user annotation or on clustering of the very data whose batch
    structure is under test.
    """

    ANNOTATION = "annotation"
    GATING = "gating"
    CLUSTERING = "clustering"


#: Confidence attached to each label tier.  Clustering-derived labels are explicitly
#: flagged as potentially batch-driven and may never on their own support a
#: recommendation (see :mod:`batchdx.cyto.decisions`).
LABEL_TIER_CONFIDENCE: dict[LabelTier, str] = {
    LabelTier.ANNOTATION: "high",
    LabelTier.GATING: "medium",
    LabelTier.CLUSTERING: "low",
}


class ComparisonSupport(StrEnum):
    """Classification of a batch comparison by the design audit (Stage 3)."""

    REFERENCE_SUPPORTED = "reference_supported"
    BIOLOGICALLY_STRATIFIED_CONDITIONAL = "biologically_stratified_conditional"
    DESCRIPTIVE_ONLY = "descriptive_only"
    NON_IDENTIFIABLE = "non_identifiable"


class MetricStatus(StrEnum):
    """Execution status of a single metric record.

    ``UNDEFINED`` is reserved for a mathematically undefined value (for example a
    relative reduction with a near-zero denominator).  It is never replaced by zero
    and never silently dropped from a summary.
    """

    OK = "ok"
    INELIGIBLE = "ineligible"
    INSUFFICIENT_DATA = "insufficient_data"
    UNDEFINED = "undefined"
    FAILED = "failed"
    NOT_IMPLEMENTED = "not_implemented"


class MetricDirection(StrEnum):
    """Which direction of a metric is the better one, if either."""

    LOWER_IS_BETTER = "lower_is_better"
    HIGHER_IS_BETTER = "higher_is_better"
    NEUTRAL = "neutral"


class CorrectionStatus(StrEnum):
    """Outcome of running one correction candidate."""

    SUCCESS = "success"
    INELIGIBLE = "ineligible"
    FAILED = "failed"
    NOT_AVAILABLE = "not_available"


class DecisionState(StrEnum):
    """Candidate-level assessment states (Stage 9).

    These are deliberately separate from :class:`MetricStatus`: a candidate whose
    metrics all executed can still be ``insufficient_evidence``.
    """

    NOT_NEEDED = "not_needed"
    ELIGIBLE_FOR_REVIEW = "eligible_for_review"
    RECOMMENDED_UNDER_POLICY = "recommended_under_policy"
    REJECTED_DISTORTION = "rejected_distortion"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    NON_IDENTIFIABLE = "non_identifiable"
    BELOW_DESIGN_FLOOR = "below_design_floor"
    EXECUTION_FAILED = "execution_failed"


class MarkerRole(StrEnum):
    """Declared role of a feature in the panel."""

    LINEAGE = "lineage"
    FUNCTIONAL = "functional"
    TECHNICAL = "technical"
    VIABILITY = "viability"
    BARCODE = "barcode"
    OTHER = "other"


class ReferenceRole(StrEnum):
    """Role a sample plays in the reference layout."""

    NONE = "none"
    TECHNICAL_REPLICATE = "technical_replicate"
    ANCHOR = "anchor"
    BIOLOGICAL_PSEUDO_CONTROL = "biological_pseudo_control"
