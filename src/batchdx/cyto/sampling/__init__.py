"""Frozen evaluation populations, deterministic sampling and comparison pairs (Stage 4)."""

from batchdx.cyto.sampling.partitions import (
    GatingRule,
    build_sampling_manifest,
    read_gating_thresholds,
    resolve_labels,
)

__all__ = [
    "GatingRule",
    "build_sampling_manifest",
    "read_gating_thresholds",
    "resolve_labels",
]
