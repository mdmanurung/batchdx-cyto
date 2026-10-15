"""The one place a :class:`MetricRecord` is stamped with frozen-partition provenance.

Every metric, in every namespace, has to carry the tier the evaluation labels came
from, their confidence, and the representation settings they were computed under.
A metric that silently omitted them could not be read next to one that did not, so
the stamping lives here rather than being retyped per metric family.
"""

from __future__ import annotations

from batchdx.cyto.metrics.registry import MetricSpec
from batchdx.cyto.schemas.contracts import MetricRecord, SamplingManifest


def metric_record(spec: MetricSpec, sampling: SamplingManifest, **kwargs) -> MetricRecord:
    """Build a record for ``spec``, carrying the frozen partition's provenance."""

    return MetricRecord(
        metric_namespace=spec.namespace,
        metric_name=spec.name,
        metric_version=spec.version,
        stage=spec.stage,
        direction=spec.direction,
        scale=spec.scale,
        label_tier=sampling.label_source,
        label_confidence=sampling.label_confidence,
        representation=sampling.representation,
        **kwargs,
    )
