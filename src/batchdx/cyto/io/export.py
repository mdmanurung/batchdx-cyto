"""Safe exports: per-sample FCS and per-candidate H5AD.

Inversion back to instrument units is attempted only where the processing manifest
declares a transformation this build can invert exactly.  Where it cannot, the
transformed values are exported and labelled as such, rather than being presented as
instrument measurements.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.io.fcs import write_fcs
from batchdx.cyto.schemas.contracts import ProcessingManifest

#: Transformations whose inverse is exact and tested.
INVERTIBLE = {"asinh"}


@dataclass(frozen=True)
class InversionPlan:
    """What the exporter will do with the values, stated explicitly."""

    invertible: bool
    scale_label: str
    note: str
    cofactor: float | None = None

    def apply(self, values: np.ndarray) -> np.ndarray:
        if not self.invertible:
            return np.array(values, copy=True)
        return np.sinh(values) * float(self.cofactor)


def plan_inversion(processing: ProcessingManifest) -> InversionPlan:
    """Decide whether the declared transformation can be undone exactly."""

    transformation = (processing.transformation or "none").lower()
    if transformation in INVERTIBLE:
        cofactor = processing.transformation_parameters.get("cofactor")
        if cofactor and float(cofactor) > 0:
            return InversionPlan(
                invertible=True, scale_label="instrument units (asinh inverted)",
                cofactor=float(cofactor),
                note=(
                    f"values were inverted with sinh(x) * {float(cofactor)}, the exact "
                    "inverse of the declared asinh transformation"
                ),
            )
        return InversionPlan(
            invertible=False, scale_label=processing.scale,
            note=(
                "the processing manifest declares an asinh transformation but no cofactor, "
                "so no exact inverse exists; transformed values are exported unchanged"
            ),
        )
    return InversionPlan(
        invertible=False, scale_label=processing.scale,
        note=(
            f"transformation {transformation!r} has no tested exact inverse in this build; "
            "transformed values are exported unchanged and are NOT instrument measurements"
        ),
    )


def export_fcs_per_sample(
    dataset: Dataset,
    matrix: np.ndarray,
    out_dir: str | Path,
    *,
    candidate_id: str,
    invert: bool = True,
) -> dict[str, Path]:
    """Write one FCS per sample with a stable channel order and a keyword policy block."""

    plan = plan_inversion(dataset.processing)
    if not invert:
        plan = InversionPlan(invertible=False, scale_label=dataset.processing.scale,
                             note="inversion was disabled by the caller")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    obs = dataset.obs
    features = list(map(str, dataset.adata.var_names))
    channels = [
        (dataset.markers[f].channel or f) for f in features
    ]
    labels = [dataset.markers[f].marker or f for f in features]

    written: dict[str, Path] = {}
    for sid in sorted(obs[SAMPLE_ID].unique()):
        mask = (obs[SAMPLE_ID] == sid).to_numpy()
        values = plan.apply(matrix[mask, :])
        path = out / f"{sid}__{candidate_id}.fcs"
        write_fcs(
            path, values, channel_names=channels, channel_labels=labels,
            scale=plan.scale_label,
            processing_note=(
                f"batchdx-cyto candidate {candidate_id}; {plan.note}. "
                "Already compensated/unmixed and QCed as declared in the processing "
                "manifest: do NOT compensate or transform these values again."
            ),
            candidate_id=candidate_id, source_sample=str(sid),
            extra_keywords={
                "BATCHDX_N_EVENTS": str(int(mask.sum())),
                "BATCHDX_FIRST_CELL_ID": str(obs.loc[mask, CELL_ID].iloc[0]),
                "BATCHDX_INVERTED": "1" if plan.invertible else "0",
            },
            overwrite=True,
        )
        written[str(sid)] = path
    return written
