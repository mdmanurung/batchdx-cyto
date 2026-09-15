"""Adapters to upstream tools, plus the native identity backend.

Importing this package registers every backend, so ``available_backends()`` reflects
what the build actually ships rather than what a configuration file hopes for.
"""

# Imported for their registration side effects; each module imports from
# .base directly, so the order here does not matter.
from batchdx.cyto.adapters import cycombine as _cycombine  # noqa: E402,F401
from batchdx.cyto.adapters import cytonorm as _cytonorm  # noqa: E402,F401
from batchdx.cyto.adapters import identity as _identity  # noqa: E402,F401
from batchdx.cyto.adapters.base import (
    BackendCapabilities,
    BackendContext,
    BackendError,
    CorrectionBackend,
    available_backends,
    candidate_signature,
    failure_result,
    get_backend,
    make_spec,
    register_backend,
    resolve_parameters,
)

__all__ = [
    "BackendCapabilities",
    "BackendContext",
    "BackendError",
    "CorrectionBackend",
    "available_backends",
    "candidate_signature",
    "failure_result",
    "get_backend",
    "make_spec",
    "register_backend",
    "resolve_parameters",
]
