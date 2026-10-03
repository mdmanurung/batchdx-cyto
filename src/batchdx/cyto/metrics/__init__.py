"""Metric registry, common-partition metrics and derived comparisons (Stages 7 and 8)."""

from batchdx.cyto.metrics.common import MIN_GROUP_CELLS, compute_common_metrics
from batchdx.cyto.metrics.deltas import compute_reductions
from batchdx.cyto.metrics.openproblems import compute_openproblems_metrics
from batchdx.cyto.metrics.preservation import (
    ContrastPlan,
    ContrastRefusal,
    ContrastSpec,
    compute_preservation_metrics,
    resolve_contrast,
)
from batchdx.cyto.metrics.registry import (
    DEFERRED_METRICS,
    OPENPROBLEMS_METRIC_NAMES,
    REIMPLEMENTED_UPSTREAM_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    MetricSpec,
    get_metric,
    register_metric,
    registered_metrics,
)

__all__ = [
    "DEFERRED_METRICS",
    "MIN_GROUP_CELLS",
    "OPENPROBLEMS_METRIC_NAMES",
    "REIMPLEMENTED_UPSTREAM_METRICS",
    "ContrastPlan",
    "ContrastRefusal",
    "ContrastSpec",
    "MetricSpec",
    "UPSTREAM_METRICS_NOT_FOUND",
    "compute_common_metrics",
    "compute_openproblems_metrics",
    "compute_preservation_metrics",
    "compute_reductions",
    "get_metric",
    "register_metric",
    "registered_metrics",
    "resolve_contrast",
]
