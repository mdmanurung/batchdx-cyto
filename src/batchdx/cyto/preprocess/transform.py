"""The fluorescence transform, with its convention stated rather than assumed.

THE CONVENTION
---------------
This build transforms with the **unscaled** arcsinh:

.. code-block:: text

    x = asinh(raw / cofactor)

That is what ``flowVS::transFlowVS`` computes, and it is what the reference
implementation applies (``.upstream/spiscy/workflow/scripts/transformation.R:159``).
The installed function's body is a single line, and its own documentation states it:

    mat[, channels[j]] = asinh(mat[, channels[j]] / cofactors[j])

with the ``Details`` section reading "``z_ij`` would be converted to
``asinh(z_ij/cofactors[i])``".  There is no multiplication by the cofactor on the way
out.

WHICH IS EXACTLY WHAT :mod:`batchdx.cyto.io.export` INVERTS
------------------------------------------------------------
:func:`asinh_invert` here is ``sinh(x) * cofactor``, which is the identity map's exact
inverse of ``asinh(raw / cofactor)``:

.. code-block:: text

    sinh(asinh(raw / c)) * c == raw

:func:`agree_with_exporter_inverse` executes that cross-check against the real
``io.export.InversionPlan.apply``, so the two code paths cannot drift apart silently.

THE CONVENTION THAT IS *NOT* THIS ONE
--------------------------------------
The scaled variant

.. code-block:: text

    x = asinh(raw / cofactor) * cofactor

is a different function, and it is **not** what anything here inverts.  Applying
``sinh(x) * cofactor`` to a value produced by the scaled variant does not recover
``raw``; at cofactor 5 and raw 1e6 the relative error is about 2.6e22.  The distinction
is recorded here because it is the trap in Tranche 3 item 3.2: reading the exporter's
``sinh(x) * cofactor`` as the inverse of a *scaled* asinh would be wrong, and the
exporter's own note would then be false.  It is not false: the exporter inverts the
unscaled form, which is the one the reference implementation produces.

LOGICLE
-------
This build implements no logicle transform, has no tested inverse for one, and the
reference implementation contains no logicle either (``grep -i logicle`` over the pinned
SPiSCy checkout returns zero matches).  ``preprocessing.transformation`` nevertheless
offers ``"logicle"`` in its ``Literal`` (``config.py:186``), so a configuration naming it
validates and then buys nothing.  That is the accepted-and-ignored outcome this package
exists to prevent, and it is left in place here only because removing a ``Literal`` member
is a central change to a shipped schema; see the work-package report.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: The forward form, as a string, for the labels this package writes into FCS keywords
#: and into a processing manifest.  Stated explicitly so a reader never has to infer
#: which of the two arcsinh conventions a ``scale:`` label is describing.
ASINH_FORM = "asinh(raw / cofactor)"

#: The inverse of :data:`ASINH_FORM`, and therefore the form
#: :mod:`batchdx.cyto.io.export` already implements.
ASINH_INVERSE_FORM = "sinh(x) * cofactor"

#: The other convention, named so that a value written on this scale is never mistaken
#: for one on :data:`ASINH_FORM`.
SCALED_ASINH_FORM = "asinh(raw / cofactor) * cofactor"


class TransformError(ValueError):
    """Raised when a transform cannot be applied to the values it was given."""


@dataclass(frozen=True)
class TransformConvention:
    """One named coordinate system, and whether this build can undo it exactly."""

    name: str
    forward: str
    inverse: str | None
    note: str

    @property
    def invertible(self) -> bool:
        return self.inverse is not None


#: The one convention this build implements.  ``cofactor`` is a per-channel value, and
#: there is no universal one: ``PreprocessingConfig`` refuses an asinh declaration that
#: does not carry a per-panel cofactor for exactly that reason.
ASINH = TransformConvention(
    name="asinh",
    forward=ASINH_FORM,
    inverse=ASINH_INVERSE_FORM,
    note=(
        "the unscaled arcsinh, identical to flowVS::transFlowVS; the exact inverse is "
        "sinh(x) * cofactor, which batchdx.cyto.io.export.InversionPlan.apply implements"
    ),
)

#: Declared, not implemented.  Listed so that a caller asking for it gets a named
#: refusal with the reason, rather than a silent pass-through of untransformed values.
LOGICLE = TransformConvention(
    name="logicle",
    forward="biexponential (logicle)",
    inverse=None,
    note=(
        "this build implements no logicle transform. The pinned SPiSCy reference "
        "contains no logicle either, and a logicle has no closed-form inverse that this "
        "package has tested, so values cannot be returned to instrument units"
    ),
)

CONVENTIONS = {c.name: c for c in (ASINH, LOGICLE)}


def resolve_convention(name: str) -> TransformConvention:
    """The named convention, or a refusal naming what this build does implement."""

    key = (name or "none").lower()
    if key == "none":
        raise TransformError(
            "no transformation was declared: the values are already in the coordinate "
            "system the processing manifest names, and there is nothing to apply"
        )
    if key not in CONVENTIONS:
        raise TransformError(
            f"transformation {name!r} is not implemented by this build; it implements "
            f"{sorted(CONVENTIONS)}"
        )
    return CONVENTIONS[key]


def _check_cofactor(cofactor: float | np.ndarray) -> float | np.ndarray:
    value = np.asarray(cofactor, dtype=np.float64)
    if not np.isfinite(value).all() or (value <= 0).any():
        raise TransformError(
            f"cofactor must be finite and strictly positive, got {cofactor!r}; there is no "
            "universal cofactor and a zero or negative one has no asinh interpretation"
        )
    return cofactor if np.ndim(cofactor) == 0 else value


def asinh_transform(values: np.ndarray, cofactor: float | np.ndarray) -> np.ndarray:
    """``asinh(raw / cofactor)`` — the convention :mod:`io.export` inverts exactly."""

    c = _check_cofactor(cofactor)
    return np.arcsinh(np.asarray(values, dtype=np.float64) / c)


def asinh_invert(values: np.ndarray, cofactor: float | np.ndarray) -> np.ndarray:
    """``sinh(x) * cofactor`` — the exact inverse of :func:`asinh_transform`."""

    c = _check_cofactor(cofactor)
    return np.sinh(np.asarray(values, dtype=np.float64)) * c


def asinh_round_trip_error(
    values: np.ndarray, cofactor: float | np.ndarray
) -> tuple[float, float]:
    """``(max_absolute, max_relative)`` error of transform-then-invert.

    A missing key is never a zero, so this returns the measured error rather than a
    boolean: a caller can state how exact the round trip was instead of asserting that
    it was.
    """

    c = _check_cofactor(cofactor)
    raw = np.asarray(values, dtype=np.float64)
    if raw.size == 0:
        raise TransformError("cannot measure a round-trip error over an empty array")
    back = asinh_invert(asinh_transform(raw, c), c)
    delta = np.abs(back - raw)
    scale = np.maximum(np.abs(raw), 1.0)
    return float(delta.max()), float((delta / scale).max())


def scaled_asinh_mismatch(
    values: np.ndarray, cofactor: float
) -> float:
    """How badly the *scaled* asinh convention breaks this build's inverse.

    Provided so the difference in item 3.2 is a number that a test asserts, rather than
    a claim in a comment.  It is not a supported transform; it is the thing that must
    never be confused with :func:`asinh_transform`.
    """

    c = _check_cofactor(cofactor)
    raw = np.asarray(values, dtype=np.float64)
    scaled = asinh_transform(raw, c) * c
    back = asinh_invert(scaled, c)
    return float((np.abs(back - raw) / np.maximum(np.abs(raw), 1.0)).max())


def agree_with_exporter_inverse(
    values: np.ndarray, cofactor: float, *, atol: float = 1e-9
) -> bool:
    """Execute the cross-check against the real exporter, not a copy of it.

    The claim under test is that :func:`asinh_invert` and
    :meth:`batchdx.cyto.io.export.InversionPlan.apply` are the same function.  A test
    that compared two literals in the same module would only prove the module is
    internally consistent, so this builds the exporter's own plan and calls its own
    method.
    """

    from batchdx.cyto.io.export import InversionPlan

    plan = InversionPlan(invertible=True, scale_label="cross-check", note="", cofactor=cofactor)
    exporter = plan.apply(np.asarray(values, dtype=np.float64))
    return bool(np.allclose(exporter, asinh_invert(values, cofactor), rtol=0.0, atol=atol))
