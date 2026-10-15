"""Upstream preprocessing: gating, the fluorescence transform, and PeacoQC.

This package is deliberately **outside** the run the pipeline performs.  The correction
workflow reads an already-preprocessed baseline and never transforms events; that
contract is stated in ``PreprocessingConfig`` and in
``_validate_processing`` (``io/dataset.py:143-152``), and nothing here changes it.

What this package is for is the step *before* the baseline exists: taking raw FCS,
removing debris, singlets and dead cells on the geometry the reference implementation
uses, applying the arcsinh the reference implementation applies, running PeacoQC, and —
the part that did not exist before — recording **which** events each of those steps
removed and why, in the four columns ``Dataset.removed_events`` has always declared.

Three modules, each refusing rather than approximating:

* :mod:`~batchdx.cyto.preprocess.gating` — two-dimensional ``rectangle`` and ``polygon``
  gates.  Not ``sampling.partitions.GatingRule``, which is one-dimensional and labels
  cells instead of removing them; the module docstring says why one-dimensional gating
  would be a strictly weaker gate under the same name.
* :mod:`~batchdx.cyto.preprocess.transform` — ``asinh(raw / cofactor)``, the convention
  :mod:`batchdx.cyto.io.export` inverts exactly, with the cross-check against the
  exporter executed rather than asserted in prose.
* :mod:`~batchdx.cyto.preprocess.peacoqc` — ``PeacoQC::PeacoQC`` from the installed
  Bioconductor package, keeping the ``!GoodCells`` index the reference implementation
  discards, and reporting an isolation tree that declined to run as a status.
"""

from __future__ import annotations

from batchdx.cyto.preprocess.gating import (
    GATE_KINDS,
    UNIMPLEMENTED_GATE_KINDS,
    Gate,
    GateTable,
    GatingError,
    GatingResult,
    apply_gates,
    gates_from_config,
    points_in_polygon,
    read_gates,
)
from batchdx.cyto.preprocess.peacoqc import (
    CRITERIA,
    IDENTITY_CHANNEL,
    VERIFIED_FORMALS,
    VERIFIED_VERSION,
    PeacoQCError,
    PeacoQCNotAvailable,
    PeacoQCResult,
    PeacoQCSettings,
    SampleQCOutcome,
    peacoqc_availability,
    read_peacoqc_settings,
    run_peacoqc,
    settings_from_config,
)
from batchdx.cyto.preprocess.removals import (
    KNOWN_STEPS,
    REMOVED_EVENT_COLUMNS,
    STEP_GATE,
    STEP_QC,
    RemovalLogError,
    check_removal_log,
    empty_removal_log,
    removal_frame,
    summarise_removals,
)
from batchdx.cyto.preprocess.transform import (
    ASINH,
    ASINH_FORM,
    ASINH_INVERSE_FORM,
    CONVENTIONS,
    LOGICLE,
    SCALED_ASINH_FORM,
    TransformConvention,
    TransformError,
    agree_with_exporter_inverse,
    asinh_invert,
    asinh_round_trip_error,
    asinh_transform,
    resolve_convention,
    scaled_asinh_mismatch,
)

__all__ = [
    "ASINH",
    "ASINH_FORM",
    "ASINH_INVERSE_FORM",
    "CONVENTIONS",
    "CRITERIA",
    "GATE_KINDS",
    "IDENTITY_CHANNEL",
    "KNOWN_STEPS",
    "LOGICLE",
    "REMOVED_EVENT_COLUMNS",
    "SCALED_ASINH_FORM",
    "STEP_GATE",
    "STEP_QC",
    "UNIMPLEMENTED_GATE_KINDS",
    "VERIFIED_FORMALS",
    "VERIFIED_VERSION",
    "Gate",
    "GatingError",
    "GatingResult",
    "GateTable",
    "PeacoQCError",
    "PeacoQCNotAvailable",
    "PeacoQCResult",
    "PeacoQCSettings",
    "RemovalLogError",
    "SampleQCOutcome",
    "TransformConvention",
    "TransformError",
    "agree_with_exporter_inverse",
    "apply_gates",
    "asinh_invert",
    "asinh_round_trip_error",
    "asinh_transform",
    "check_removal_log",
    "empty_removal_log",
    "gates_from_config",
    "peacoqc_availability",
    "points_in_polygon",
    "read_gates",
    "read_peacoqc_settings",
    "removal_frame",
    "resolve_convention",
    "run_peacoqc",
    "scaled_asinh_mismatch",
    "settings_from_config",
    "summarise_removals",
]
