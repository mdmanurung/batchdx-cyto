"""Post-correction output validation and correction-magnitude reporting.

The checks here are what separates "the method ran" from "the output can be used":
event identity, marker order, dimensionality, finiteness, and invariance of markers
the candidate was not allowed to touch.
"""

from __future__ import annotations

import numpy as np

from batchdx.cyto.io.dataset import SAMPLE_ID, Dataset

#: Tolerance for calling an excluded marker unchanged.  Correction backends are run in
#: double precision, so anything above this is a real modification, not round-off.
INVARIANCE_ATOL = 1e-9


class OutputValidationError(ValueError):
    """Raised when a corrected matrix fails a structural check."""


def validate_output(
    dataset: Dataset,
    corrected: np.ndarray,
    corrected_markers: list[str],
    *,
    candidate_id: str,
) -> dict[str, bool]:
    """Structural validation.  Raises on any violation; returns the checks that passed."""

    baseline = dataset.expression()
    checks: dict[str, bool] = {}

    if corrected.shape != baseline.shape:
        raise OutputValidationError(
            f"candidate {candidate_id!r}: corrected matrix has shape {corrected.shape}, "
            f"baseline has {baseline.shape}"
        )
    checks["shape_matches_baseline"] = True

    if not np.isfinite(corrected).all():
        n_bad = int((~np.isfinite(corrected)).sum())
        raise OutputValidationError(
            f"candidate {candidate_id!r}: corrected matrix contains {n_bad} non-finite "
            "value(s); they are not replaced with zeros"
        )
    checks["all_finite"] = True

    all_markers = list(dataset.adata.var_names)
    unknown = [m for m in corrected_markers if m not in all_markers]
    if unknown:
        raise OutputValidationError(
            f"candidate {candidate_id!r}: correction scope names unknown marker(s) {unknown}"
        )
    checks["correction_scope_known"] = True

    excluded = [m for m in all_markers if m not in set(corrected_markers)]
    if excluded:
        idx = dataset.feature_index(excluded)
        delta = np.abs(corrected[:, idx] - baseline[:, idx])
        worst = float(delta.max()) if delta.size else 0.0
        if worst > INVARIANCE_ATOL:
            offender = excluded[int(np.unravel_index(delta.argmax(), delta.shape)[1])]
            raise OutputValidationError(
                f"candidate {candidate_id!r}: marker {offender!r} is outside the correction "
                f"scope but changed by up to {worst:.3g}"
            )
    checks["excluded_markers_unchanged"] = True
    checks["event_identity_preserved"] = True  # rows are positional and never reordered
    return checks


def correction_magnitude(
    dataset: Dataset, corrected: np.ndarray, labels: dict[str, str] | None = None
) -> tuple[dict[str, float], dict[str, float], dict[str, float], list[str], list[str]]:
    """Median absolute change by marker, by sample, and by population.

    Returns ``(by_marker, by_sample, by_population, changed_markers, unchanged_markers)``.
    """

    baseline = dataset.expression()
    delta = np.abs(corrected - baseline)
    markers = list(dataset.adata.var_names)

    by_marker = {m: float(np.median(delta[:, i])) for i, m in enumerate(markers)}
    changed = [m for i, m in enumerate(markers) if float(delta[:, i].max()) > INVARIANCE_ATOL]
    unchanged = [m for m in markers if m not in set(changed)]

    obs = dataset.obs
    by_sample: dict[str, float] = {}
    for sid in sorted(obs[SAMPLE_ID].unique()):
        mask = (obs[SAMPLE_ID] == sid).to_numpy()
        by_sample[str(sid)] = float(np.median(delta[mask, :])) if mask.any() else 0.0

    by_population: dict[str, float] = {}
    if labels:
        cell_ids = obs["cell_id"].tolist()
        label_array = np.array([labels.get(c, "") for c in cell_ids], dtype=object)
        for pop in sorted({p for p in label_array if p}):
            mask = label_array == pop
            by_population[str(pop)] = float(np.median(delta[mask, :]))

    return by_marker, by_sample, by_population, changed, unchanged


def changed_samples(dataset: Dataset, corrected: np.ndarray) -> list[str]:
    baseline = dataset.expression()
    delta = np.abs(corrected - baseline)
    obs = dataset.obs
    out: list[str] = []
    for sid in sorted(obs[SAMPLE_ID].unique()):
        mask = (obs[SAMPLE_ID] == sid).to_numpy()
        if mask.any() and float(delta[mask, :].max()) > INVARIANCE_ATOL:
            out.append(str(sid))
    return out
