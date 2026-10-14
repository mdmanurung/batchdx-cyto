"""CytoNorm correction backend — declared, not yet implemented (Tier B).

The Stage 0 audit established the API this adapter would target: ``CytoNorm.train``
takes *FCS file paths*, per-file batch labels, channels and a ``transformList``, and
``CytoNorm.normalize`` applies the stored model, so fitting and application are
genuinely separable.  That file-based interface means a working adapter has to write
per-sample FCS from the baseline, hand CytoNorm a transform list consistent with the
declared processing manifest, and read the normalised files back.

None of that is implemented or tested in this build, and neither the CytoNorm package
nor its Bioconductor dependencies could be installed in the environment where this
code was developed (see ``docs/blocked.md``).  Rather than ship a path that looks
complete, the backend registers itself and reports an explicit deferred status, so a
configuration naming it gets a reason instead of a fabricated correction.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from batchdx.cyto.adapters.base import (
    BackendCapabilities,
    BackendContext,
    BackendError,
    BackendNotAvailable,
    register_backend,
)
from batchdx.cyto.adapters.r_runner import r_package_version, r_version
from batchdx.cyto.io.dataset import Dataset
from batchdx.cyto.schemas.contracts import CorrectionResult, SoftwareVersion

DEFERRED_REASON = (
    "the cytonorm backend is declared but not implemented in this build (Tier B). "
    "CytoNorm.train/CytoNorm.normalize operate on FCS files, which requires a tested "
    "FCS round trip and a transformList consistent with the declared processing "
    "manifest; neither is implemented, and the package could not be installed here. "
    "See docs/implementation_status.md and docs/blocked.md."
)


class CytoNormBackend:
    """Capabilities resolved from the pinned upstream source; execution deferred."""

    capabilities = BackendCapabilities(
        name="cytonorm",
        required_sample_fields=("batch", "reference_id", "reference_role"),
        reference_layouts=("measured_reference_in_every_batch", "control_free"),
        supports_covariates=False,
        supports_marker_mask=True,
        input_scale="transformed, with the transform declared to CytoNorm as a transformList",
        output_scale="same transformed scale as the input",
        supports_fit_apply=True,
        fit_apply_verified=False,
        resource_profile="memory: per-file; the model is fitted from control files only",
        accepted_parameters={
            "nCells": 1_000_000,
            "xdim": 15,
            "ydim": 15,
            "nClus": 10,
            "nQ": 99,
            "scale": False,
        },
        rejected_parameters={
            "aggregate_control": (
                "requiring a measured control in every batch is a batchdx policy choice "
                "that restricts the 'cytonorm' backend to its strongest evidence tier; "
                "it is not an upstream constraint. CytoNorm 2.0 itself supports "
                "control-free normalisation (QuantileNorm.train goal='mean', the "
                "default, or a numeric goal distribution), which the design audit "
                "exposes as the separate, weaker 'cytonorm_control_free' tier. An "
                "artificial aggregate control is additionally outside the validated "
                "path and is not offered as a standard option under either tier."
            ),
        },
        notes=DEFERRED_REASON,
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [
            r_version(context.rscript),
            r_package_version("CytoNorm", context.rscript),
            r_package_version("FlowSOM", context.rscript),
            r_package_version("flowCore", context.rscript),
        ]

    def run(
        self,
        dataset: Dataset,
        candidate_id: str,
        parameters: dict[str, Any],
        context: BackendContext,
    ) -> tuple[np.ndarray, CorrectionResult]:
        raise BackendNotAvailable(DEFERRED_REASON)


register_backend(CytoNormBackend())
