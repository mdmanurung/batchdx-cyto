"""The ``none`` backend: a genuine no-correction baseline, not a placeholder.

"No correction" is a scientific result in its own right.  This backend returns the
preprocessed baseline bit-for-bit and goes through exactly the same validation,
export and scoring path as every other candidate, so a comparison against it is a
real comparison.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from batchdx.cyto._version import __version__
from batchdx.cyto.adapters.base import (
    BackendCapabilities,
    BackendContext,
    make_spec,
    register_backend,
    resolve_parameters,
)
from batchdx.cyto.adapters.validate import correction_magnitude, validate_output
from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.io.hashing import matrix_hash
from batchdx.cyto.schemas.contracts import CorrectionResult, ProvenanceRef, SoftwareVersion
from batchdx.cyto.schemas.enums import CorrectionStatus


class IdentityBackend:
    """Return the baseline unchanged."""

    capabilities = BackendCapabilities(
        name="none",
        required_sample_fields=(),
        reference_layouts=("any", "none"),
        supports_covariates=False,
        supports_marker_mask=False,
        input_scale="as declared in the processing manifest",
        output_scale="as declared in the processing manifest",
        supports_fit_apply=True,
        fit_apply_verified=True,
        resource_profile="memory: streams; performs no computation",
        accepted_parameters={},
        notes="the identity candidate; always benchmarked so that 'no correction' can win",
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [SoftwareVersion(name="batchdx-cyto", version=__version__, source="python")]

    def run(
        self,
        dataset: Dataset,
        candidate_id: str,
        parameters: dict[str, Any],
        context: BackendContext,
    ) -> tuple[np.ndarray, CorrectionResult]:
        resolved = resolve_parameters(self.capabilities, parameters)
        matrix = dataset.expression()

        spec = make_spec(candidate_id, "none", __version__, resolved, context.seed,
                         dataset, correction_markers=[])
        checks = validate_output(dataset, matrix, corrected_markers=[],
                                 candidate_id=candidate_id)
        by_marker, by_sample, by_population, changed, unchanged = correction_magnitude(
            dataset, matrix,
            labels=context.sampling.labels_by_cell if context.sampling else None,
        )
        checks["identity_numerically_unchanged"] = bool(
            np.array_equal(matrix, dataset.expression())
        )

        result = CorrectionResult(
            candidate=spec,
            status=CorrectionStatus.SUCCESS,
            output_scale=dataset.processing.scale,
            output_hash=matrix_hash(matrix),
            n_cells=dataset.n_cells,
            markers_changed=changed,
            markers_unchanged=unchanged,
            samples_changed=[],
            populations_changed=[],
            magnitude_by_marker=by_marker,
            magnitude_by_sample=by_sample,
            magnitude_by_population=by_population,
            validation=checks,
            software=self.availability(context),
            provenance=[ProvenanceRef(kind="baseline", identifier=dataset.manifest.dataset_id,
                                      content_hash=dataset.baseline_hash())],
        )
        return matrix, result


register_backend(IdentityBackend())
