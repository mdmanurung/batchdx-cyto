"""Metric registry, common-partition metrics and derived comparisons (Stage 7)."""

from batchdx.cyto.metrics.common import MIN_GROUP_CELLS, compute_common_metrics
from batchdx.cyto.metrics.deltas import compute_reductions
from batchdx.cyto.metrics.registry import (
    DEFERRED_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    MetricSpec,
    get_metric,
    register_metric,
    registered_metrics,
)

__all__ = [
    "DEFERRED_METRICS",
    "MIN_GROUP_CELLS",
    "MetricSpec",
    "UPSTREAM_METRICS_NOT_FOUND",
    "compute_common_metrics",
    "compute_reductions",
    "get_metric",
    "register_metric",
    "registered_metrics",
]
