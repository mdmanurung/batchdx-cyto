"""Metric registry, common-partition metrics and derived comparisons (Stage 7)."""

from batchdx.cyto.metrics.biology import compute_biological_contrasts
from batchdx.cyto.metrics.common import MIN_GROUP_CELLS, compute_common_metrics
from batchdx.cyto.metrics.deltas import compute_reductions
from batchdx.cyto.metrics.registry import (
    BIOLOGICAL_CONTRAST,
    DEFERRED_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    MetricSpec,
    get_metric,
    register_metric,
    registered_metrics,
)
from batchdx.cyto.metrics.uncertainty import (
    attach_emd_uncertainty,
    donor_block_bootstrap_emd,
    restricted_permutation_test,
)

__all__ = [
    "BIOLOGICAL_CONTRAST",
    "DEFERRED_METRICS",
    "MIN_GROUP_CELLS",
    "MetricSpec",
    "UPSTREAM_METRICS_NOT_FOUND",
    "attach_emd_uncertainty",
    "compute_biological_contrasts",
    "compute_common_metrics",
    "compute_reductions",
    "donor_block_bootstrap_emd",
    "get_metric",
    "register_metric",
    "registered_metrics",
    "restricted_permutation_test",
]
