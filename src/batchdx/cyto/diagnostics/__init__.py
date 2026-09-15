"""Detectors: native upstream assessments plus the drift and composition diagnostics."""

# Imported for their registration side effects; each module imports from
# .base directly, so the order here does not matter.
from batchdx.cyto.diagnostics import composition as _composition  # noqa: E402,F401
from batchdx.cyto.diagnostics import cycombine_detect as _cycombine_detect  # noqa: E402,F401
from batchdx.cyto.diagnostics import cytobatchflagr as _cytobatchflagr  # noqa: E402,F401
from batchdx.cyto.diagnostics import drift as _drift  # noqa: E402,F401
from batchdx.cyto.diagnostics.base import (
    Detector,
    DetectorContext,
    available_detectors,
    eligibility_for,
    get_detector,
    register_detector,
)

__all__ = [
    "Detector",
    "DetectorContext",
    "available_detectors",
    "eligibility_for",
    "get_detector",
    "register_detector",
]
