"""Validation-fixture correction backends.

These are **not** user correction methods and are deliberately not registered by the
package itself: they exist so that the scientific suite can state exactly what a
metric or a guardrail is supposed to catch.

``fixture_offset_removal``
    a genuine, simple correction: subtract each batch's offset from the measured
    reference samples.  On the synthetic paired-reference fixture this is the right
    answer, so metrics must show it improving alignment without distorting anything.

``fixture_collapse``
    replaces every value with the population's median taken across all batches.  Every
    batch becomes the same point mass, so alignment is perfect and all spread is
    destroyed; a system that rewards alignment alone will recommend it, which is
    exactly the failure the guardrails must prevent.

``fixture_joint_scramble``
    permutes each marker independently within each (population, batch) group.  Every
    marginal distribution is preserved *exactly* while within-population co-expression
    is destroyed, so only a joint-structure metric can detect it.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from batchdx.cyto.adapters.base import (
    BackendCapabilities,
    BackendContext,
    BackendError,
    make_spec,
    register_backend,
    resolve_parameters,
)
from batchdx.cyto.adapters.validate import (
    changed_samples,
    correction_magnitude,
    validate_output,
)
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.io.hashing import matrix_hash
from batchdx.cyto.schemas.contracts import CorrectionResult, SoftwareVersion
from batchdx.cyto.schemas.enums import CorrectionStatus, ReferenceRole

FIXTURE_NOTE = (
    "validation fixture: an artificially destructive or idealised control used by the "
    "test suite. Never offer this as a user correction method."
)


def _finish(
    dataset: Dataset,
    matrix: np.ndarray,
    candidate_id: str,
    backend_name: str,
    parameters: dict[str, Any],
    context: BackendContext,
    markers: list[str],
) -> CorrectionResult:
    checks = validate_output(dataset, matrix, markers, candidate_id=candidate_id)
    labels = context.sampling.labels_by_cell if context.sampling else None
    by_marker, by_sample, by_population, changed, unchanged = correction_magnitude(
        dataset, matrix, labels=labels
    )
    spec = make_spec(candidate_id, backend_name, "fixture", parameters, context.seed,
                     dataset, markers)
    return CorrectionResult(
        candidate=spec, status=CorrectionStatus.SUCCESS,
        output_scale=dataset.processing.scale, output_hash=matrix_hash(matrix),
        n_cells=dataset.n_cells, markers_changed=changed, markers_unchanged=unchanged,
        samples_changed=changed_samples(dataset, matrix),
        populations_changed=sorted(by_population),
        magnitude_by_marker=by_marker, magnitude_by_sample=by_sample,
        magnitude_by_population=by_population, validation=checks,
        warnings=[FIXTURE_NOTE],
        software=[SoftwareVersion(name="test-fixture", version="1", source="tests")],
    )


def _batch_array(dataset: Dataset) -> np.ndarray:
    batch_of = dataset.sample_batches()
    return dataset.obs[SAMPLE_ID].map(batch_of).astype(str).to_numpy()


def _label_array(dataset: Dataset, context: BackendContext) -> np.ndarray:
    if context.sampling is None:
        return np.full(dataset.n_cells, "all", dtype=object)
    labels = context.sampling.labels_by_cell
    return dataset.obs[CELL_ID].map(lambda c: labels.get(c, "unselected")).to_numpy()


class OffsetRemovalBackend:
    """Subtract the offset each batch's measured reference shows against the grand mean."""

    capabilities = BackendCapabilities(
        name="fixture_offset_removal",
        required_sample_fields=("batch", "reference_id", "reference_role"),
        reference_layouts=("measured_reference_in_every_batch",),
        supports_marker_mask=True,
        supports_fit_apply=True,
        fit_apply_verified=False,
        accepted_parameters={},
        notes=FIXTURE_NOTE,
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [SoftwareVersion(name="test-fixture", version="1", source="tests")]

    def run(self, dataset, candidate_id, parameters, context):
        resolved = resolve_parameters(self.capabilities, parameters)
        markers = dataset.marker_ids("correction")
        idx = dataset.feature_index(markers)
        matrix = dataset.expression()
        batches = _batch_array(dataset)

        reference_samples = [
            sid for sid, rec in dataset.samples.items()
            if rec.reference_role in (ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR)
        ]
        if not reference_samples:
            raise BackendError(
                "fixture_offset_removal needs a measured reference sample in every batch"
            )
        is_reference = dataset.obs[SAMPLE_ID].isin(reference_samples).to_numpy()

        grand = matrix[np.ix_(is_reference, idx)].mean(axis=0)
        for batch in sorted(set(batches)):
            mask = (batches == batch) & is_reference
            if not mask.any():
                raise BackendError(f"batch {batch!r} contains no measured reference sample")
            offset = matrix[np.ix_(mask, idx)].mean(axis=0) - grand
            rows = np.flatnonzero(batches == batch)
            matrix[np.ix_(rows, idx)] -= offset

        return matrix, _finish(dataset, matrix, candidate_id, self.capabilities.name,
                               resolved, context, markers)


class CollapseBackend:
    """Replace every value with the population's median taken across all batches.

    The result aligns every batch perfectly (each batch becomes the same point mass)
    while destroying all within-population spread, which is exactly the degenerate
    "correction" a system that rewards alignment alone would recommend.
    """

    capabilities = BackendCapabilities(
        name="fixture_collapse",
        supports_marker_mask=True,
        accepted_parameters={},
        notes=FIXTURE_NOTE,
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [SoftwareVersion(name="test-fixture", version="1", source="tests")]

    def run(self, dataset, candidate_id, parameters, context):
        resolved = resolve_parameters(self.capabilities, parameters)
        markers = dataset.marker_ids("correction")
        idx = dataset.feature_index(markers)
        matrix = dataset.expression()
        labels = _label_array(dataset, context)

        for population in sorted(set(labels)):
            rows = np.flatnonzero(labels == population)
            if len(rows) == 0:
                continue
            matrix[np.ix_(rows, idx)] = np.median(matrix[np.ix_(rows, idx)], axis=0)

        return matrix, _finish(dataset, matrix, candidate_id, self.capabilities.name,
                               resolved, context, markers)


class JointScrambleBackend:
    """Permute each marker independently within (population, batch).

    Every marginal distribution is preserved exactly; the joint distribution is not.
    """

    capabilities = BackendCapabilities(
        name="fixture_joint_scramble",
        supports_marker_mask=True,
        accepted_parameters={},
        notes=FIXTURE_NOTE,
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [SoftwareVersion(name="test-fixture", version="1", source="tests")]

    def run(self, dataset, candidate_id, parameters, context):
        resolved = resolve_parameters(self.capabilities, parameters)
        markers = dataset.marker_ids("correction")
        idx = dataset.feature_index(markers)
        matrix = dataset.expression()
        batches = _batch_array(dataset)
        labels = _label_array(dataset, context)
        rng = np.random.default_rng(context.seed)

        for population in sorted(set(labels)):
            for batch in sorted(set(batches)):
                rows = np.flatnonzero((labels == population) & (batches == batch))
                if len(rows) < 2:
                    continue
                for column in idx:
                    matrix[rows, column] = matrix[rng.permutation(rows), column]

        return matrix, _finish(dataset, matrix, candidate_id, self.capabilities.name,
                               resolved, context, markers)


class ExplodingBackend:
    """Always raises, so the workflow's failure path can be exercised."""

    capabilities = BackendCapabilities(
        name="fixture_failure", accepted_parameters={}, notes=FIXTURE_NOTE
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [SoftwareVersion(name="test-fixture", version="1", source="tests")]

    def run(self, dataset, candidate_id, parameters, context):
        raise BackendError("fixture_failure always fails, on purpose")


_REGISTERED = False


def register_fixture_backends() -> None:
    """Register the fixture backends exactly once."""

    global _REGISTERED
    if _REGISTERED:
        return
    for backend in (OffsetRemovalBackend(), CollapseBackend(), JointScrambleBackend(),
                    ExplodingBackend()):
        register_backend(backend)
    _REGISTERED = True
