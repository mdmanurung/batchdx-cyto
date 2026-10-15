"""Two-dimensional preliminary gating: the geometry the reference implementation uses.

WHY THIS IS NOT ``sampling.partitions.GatingRule``
---------------------------------------------------
``GatingRule`` (``batchdx.cyto.sampling.partitions:46-53``) is a **one-dimensional**
per-marker threshold with a sign, and it *labels* cells: it assigns an evaluation
population, first match by priority, and a cell that matches nothing is reported as
``ungated`` rather than dropped.  It cannot express a polygon, and it is not trying to.

The preliminary gate in the reference implementation is a different thing entirely.
``prelim_gating.R:215-216`` builds ``rectangleGate``/``polygonGate`` objects and calls
``flowCore::Subset``, which **removes** events.  The shipped gate set
(``.upstream/spiscy/config/02_prelim_gating.yaml``) is:

===================  ==========  ==========================  ==================================
gate                 type        channels                    what it is for
===================  ==========  ==========================  ==================================
``debris``           rectangle   ``FSC-A`` x ``SSC-A``      scatter window
``singlets``         polygon     ``FSC-A`` x ``FSC-H``      the diagonal doublet band
``live``             rectangle   ``FVS780`` (one axis)      viability
===================  ==========  ==========================  ==================================

Two consequences, and both are load-bearing:

* the ``singlets`` gate is a **polygon** on a scatter-vs-scatter plot.  Expressed as two
  independent one-dimensional thresholds it becomes a box, which is not the same set of
  events: the diagonal boundary is precisely what separates singlets from doublets, and
  an axis-aligned box either admits the doublet shoulder or discards real singlets.
  Shipping one-dimensional gating here would be shipping a strictly weaker gate under
  the same name;
* the gate channels (``FSC-A``, ``SSC-A``, ``FSC-H``, ``FVS780``) are **scatter and
  viability parameters, not panel markers**.  ``GatingRule`` addresses features by
  ``feature_id``, so it cannot even name them unless they have been declared as
  features.

DECISION
--------
A two-dimensional gate table, in the reference implementation's own YAML shape, so the
pinned ``02_prelim_gating.yaml`` is directly readable and any divergence is visible as a
parse error rather than as a silent reinterpretation.  ``rectangle`` and ``polygon`` are
implemented.  ``ellipsoid``, ``quad`` and ``polytope`` are **refused by name**: the
reference implementation has them commented out as TODOs, so there is no upstream
behaviour to match, and a gate this build cannot evaluate must not be read as a gate that
passed.

A one-axis ``rectangle`` is accepted, because the reference implementation's own ``live``
gate is exactly that.  A ``polygon`` on one axis is refused: a polygon needs two.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import yaml

#: Gate kinds this build evaluates.  Everything else is refused by name, with the
#: reason, rather than treated as a gate that passed.
GATE_KINDS = ("rectangle", "polygon")

#: The kinds the reference implementation leaves commented out at
#: ``prelim_gating.R:119-149``.  Listed so the refusal can say *why* rather than only
#: that it is refused.
UNIMPLEMENTED_GATE_KINDS = {
    "ellipsoid": "the reference implementation has it commented out as a TODO "
                 "(prelim_gating.R:120-128), so there is no upstream behaviour to match",
    "quad": "the reference implementation has it commented out as a TODO "
            "(prelim_gating.R:130-134)",
    "polytope": "the reference implementation has it commented out as a TODO "
                "(prelim_gating.R:136-149)",
}


class GatingError(ValueError):
    """Raised when a gate table cannot be read, or a gate cannot be evaluated."""


@dataclass(frozen=True)
class Gate:
    """One gate: a name, a kind, and the coordinates it is drawn with.

    ``bounds`` is a tuple of ``(low, high)`` intervals, one per axis, for a
    ``rectangle``.  ``vertices`` is an ``(n, 2)`` array of polygon vertices in
    ``(mk1, mk2)`` order for a ``polygon``.  Exactly one of the two is populated, and
    which one is decided by ``kind`` rather than by which happens to be non-empty.
    """

    name: str
    kind: Literal["rectangle", "polygon"]
    mk1: str
    mk2: str | None = None
    bounds: tuple[tuple[float, float], ...] = ()
    vertices: np.ndarray | None = field(default=None, repr=False)

    @property
    def channels(self) -> tuple[str, ...]:
        return (self.mk1,) if self.mk2 is None else (self.mk1, self.mk2)

    @property
    def reason(self) -> str:
        """The criterion this gate names, as it appears in a removal log."""

        axes = "/".join(self.channels)
        if self.kind == "rectangle":
            spans = ",".join(f"[{lo:g},{hi:g}]" for lo, hi in self.bounds)
            return f"{self.name}:rectangle({axes}={spans})"
        assert self.vertices is not None
        return f"{self.name}:polygon({axes},{len(self.vertices)} vertices)"

    def contains(self, x: np.ndarray, y: np.ndarray | None = None) -> np.ndarray:
        """Boolean mask of the events this gate **keeps**.

        The name is deliberate: a gate's job is to keep what is inside it, and the
        reference implementation's ``Subset`` removes everything else.  Nothing here
        ever returns the complement silently.
        """

        if self.kind == "rectangle":
            if y is None:
                low, high = self.bounds[0]
                return (x >= low) & (x <= high)
            keep = np.ones(x.shape, dtype=bool)
            for axis, (low, high) in zip((x, y), self.bounds, strict=True):
                keep &= (axis >= low) & (axis <= high)
            return keep
        if y is None:
            raise GatingError(
                f"gate {self.name!r} is a polygon, which needs two channels; it declares "
                f"only {self.mk1!r}. A one-dimensional gate is declared as "
                "'type: rectangle'"
            )
        assert self.vertices is not None
        return points_in_polygon(x, y, self.vertices)

    def validate(self) -> None:
        """Refuse a gate this build cannot evaluate, before it is applied to anything."""

        if self.kind not in GATE_KINDS:  # pragma: no cover - read_gates checks first
            raise GatingError(f"gate {self.name!r}: unknown kind {self.kind!r}")
        if self.kind == "rectangle":
            if not self.bounds:
                raise GatingError(
                    f"gate {self.name!r}: a rectangle needs one (low, high) interval per "
                    "channel"
                )
            if self.mk2 is None and len(self.bounds) != 1:
                raise GatingError(
                    f"gate {self.name!r}: declares one channel {self.mk1!r} but "
                    f"{len(self.bounds)} intervals"
                )
            if self.mk2 is not None and len(self.bounds) != 2:
                raise GatingError(
                    f"gate {self.name!r}: declares two channels "
                    f"({self.mk1}, {self.mk2}) but {len(self.bounds)} intervals"
                )
            for low, high in self.bounds:
                if not (np.isfinite(low) and np.isfinite(high)) or high < low:
                    raise GatingError(
                        f"gate {self.name!r}: interval [{low}, {high}] is not a finite, "
                        "ordered range"
                    )
            return
        if self.vertices is None or self.vertices.ndim != 2 or self.vertices.shape[1] != 2:
            raise GatingError(
                f"gate {self.name!r}: a polygon needs an (n, 2) list of (mk1, mk2) vertices"
            )
        if self.vertices.shape[0] < 3:
            raise GatingError(
                f"gate {self.name!r}: a polygon needs at least three vertices, got "
                f"{self.vertices.shape[0]}"
            )
        if not np.isfinite(self.vertices).all():
            raise GatingError(f"gate {self.name!r}: polygon vertices must be finite")
        if self.mk2 is None:
            raise GatingError(
                f"gate {self.name!r}: a polygon needs two channels; it declares only "
                f"{self.mk1!r}"
            )


def points_in_polygon(
    x: np.ndarray, y: np.ndarray, vertices: np.ndarray
) -> np.ndarray:
    """Even-odd (crossing-number) test, vectorised over events.

    Written out rather than delegated to a plotting package so that a gate has no
    dependency this package would have to declare, and so the boundary convention is
    visible: a point exactly on an edge may be counted either way by this rule, which is
    the same freedom ``flowCore::polygonGate`` does not take.  For cytometry scatter
    data, where coordinates are continuous, the boundary set has measure zero; the fact
    is recorded rather than pretended away.
    """

    v = np.asarray(vertices, dtype=np.float64)
    px = np.asarray(x, dtype=np.float64)
    py = np.asarray(y, dtype=np.float64)
    inside = np.zeros(px.shape, dtype=bool)
    n = v.shape[0]
    j = n - 1
    for i in range(n):
        xi, yi = v[i]
        xj, yj = v[j]
        straddles = (yi > py) != (yj > py)
        if not straddles.any():
            j = i
            continue
        # x coordinate of the edge at this point's y, only where the edge straddles it
        with np.errstate(divide="ignore", invalid="ignore"):
            crossing = (xj - xi) * (py - yi) / (yj - yi) + xi
        inside ^= straddles & (px < crossing)
        j = i
    return inside


@dataclass(frozen=True)
class GateTable:
    """An ordered chain of gates, plus the channels they need."""

    gates: tuple[Gate, ...]
    channels: tuple[str, ...]

    def validate(self, available: tuple[str, ...] | list[str] = ()) -> None:
        for gate in self.gates:
            gate.validate()
        missing = [c for c in self.channels if c not in set(available)]
        if missing:
            raise GatingError(
                f"the gate table needs channel(s) {missing}, which the data does not carry. "
                "Gate channels are scatter and viability parameters, not panel markers: "
                "declare them in the marker manifest or supply an FCS that has them. A "
                "missing channel is not zero, and a gate on a constant zero would remove "
                "either everything or nothing"
            )

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "name": g.name,
                    "kind": g.kind,
                    "mk1": g.mk1,
                    "mk2": g.mk2 or "",
                    "channels": "/".join(g.channels),
                    "reason": g.reason,
                }
                for g in self.gates
            ]
        )


# --------------------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------------------


def _as_coordinates(kind: str, name: str, raw) -> tuple[tuple[tuple[float, float], ...], np.ndarray | None]:
    if kind == "rectangle":
        if not isinstance(raw, dict):
            raise GatingError(
                f"gate {name!r}: a rectangle's 'coords' must be a mapping of channel -> "
                "[low, high], as in the reference implementation's config"
            )
        out: list[tuple[float, float]] = []
        for channel, span in raw.items():
            values = list(span)
            if len(values) != 2:
                raise GatingError(
                    f"gate {name!r}: coords for {channel!r} must be [low, high], got "
                    f"{span!r}"
                )
            out.append((float(values[0]), float(values[1])))
        return tuple(out), None
    values = np.asarray(raw, dtype=np.float64)
    return (), values


def read_gates(path: str | Path) -> GateTable:
    """Read a gate table in the reference implementation's own YAML shape.

    The shape is ``gate_chain`` (the order to apply) plus ``gates`` keyed by name, each
    entry carrying ``filterId``, ``mk1``, optional ``mk2``, ``type`` and ``coords``.  It
    is read rather than reinterpreted, so a table written for ``prelim_gating.R`` loads
    here unchanged and a divergence shows up as a parse error.

    The order is the chain's, and it matters: the reference implementation applies its
    gates in sequence to the surviving events, so a point outside the first gate is
    never offered to the second.  The order is preserved here for the same reason.
    """

    p = Path(path)
    if not p.exists():
        raise GatingError(f"gate table not found: {p}")
    try:
        raw = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise GatingError(f"{p}: not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise GatingError(f"{p}: top level of a gate table must be a mapping")
    declared = raw.get("gates")
    if not isinstance(declared, dict) or not declared:
        raise GatingError(f"{p}: 'gates' must be a non-empty mapping of gate name to geometry")
    chain = raw.get("gate_chain")
    if chain is None:
        names = list(declared)
    else:
        names = [str(n) for n in chain]
        unknown = [n for n in names if n not in declared]
        if unknown:
            raise GatingError(
                f"{p}: gate_chain names gate(s) {unknown} that 'gates' does not define"
            )
        unused = [n for n in declared if n not in set(names)]
        if unused:
            raise GatingError(
                f"{p}: gate(s) {unused} are defined but not in gate_chain. A gate that is "
                "declared and never applied is the accepted-and-ignored outcome this "
                "package refuses; list it in the chain or delete it"
            )

    gates: list[Gate] = []
    for name in names:
        spec = declared[name] or {}
        kind = str(spec.get("type", "")).strip().lower()
        if kind in UNIMPLEMENTED_GATE_KINDS:
            raise GatingError(
                f"gate {name!r}: type {kind!r} is not implemented by this build, because "
                f"{UNIMPLEMENTED_GATE_KINDS[kind]}. Declaring it and then not applying it "
                "would report a gate that was never evaluated as a gate that passed"
            )
        if kind not in GATE_KINDS:
            raise GatingError(
                f"gate {name!r}: type {kind!r} is not one of {list(GATE_KINDS)}"
            )
        mk1 = spec.get("mk1")
        if not mk1:
            raise GatingError(f"gate {name!r}: 'mk1' is required and is never inferred")
        mk2 = spec.get("mk2")
        bounds, vertices = _as_coordinates(kind, name, spec.get("coords"))
        gate = Gate(
            name=str(spec.get("filterId", name)),
            kind=kind,  # type: ignore[arg-type]
            mk1=str(mk1),
            mk2=str(mk2) if mk2 else None,
            bounds=bounds,
            vertices=vertices,
        )
        gate.validate()
        gates.append(gate)

    channels: list[str] = []
    for gate in gates:
        for channel in gate.channels:
            if channel not in channels:
                channels.append(channel)
    return GateTable(gates=tuple(gates), channels=tuple(channels))


def gates_from_config(config, base_dir: str | Path = ".") -> GateTable | None:
    """The gate table ``preprocessing.gating`` declares, or ``None`` if it declares none.

    This is the reader that makes ``preprocessing.gating`` a key rather than a promise.
    It reads the key and nothing else: the gate channels come from the table, the order
    comes from the chain, and no channel is inferred from the panel.
    """

    declared = config.preprocessing.gating
    if declared is None:
        return None
    return read_gates(Path(base_dir) / str(declared))


# --------------------------------------------------------------------------------------
# Applying
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class GateOutcome:
    """What one gate did, in numbers and in per-sample counts."""

    gate: Gate
    kept: int
    removed: int
    n_input: int
    removed_by_sample: dict[str, int]


@dataclass(frozen=True)
class GatingResult:
    """The keep-mask, the per-gate outcome, and the removal log rows.

    ``keep`` is a boolean array over the events *as they were handed in*.  A caller that
    gates twice must gate the surviving rows, not the original ones, or the mask indexes
    the wrong events; :func:`apply_gates` is the only supported way to chain.
    """

    keep: np.ndarray
    outcomes: tuple[GateOutcome, ...]
    removed: pd.DataFrame


def apply_gates(
    table: GateTable,
    values: np.ndarray,
    channels: list[str] | tuple[str, ...],
    *,
    sample_ids: np.ndarray,
    original_event_index: np.ndarray,
    step: str = "prelim_gating",
) -> GatingResult:
    """Apply every gate in order and return what survived plus why the rest did not.

    ``values`` is ``(n_events, n_channels)`` over ``channels``.  The gates are applied
    to the events still alive at each step, exactly as a sequence of ``Subset`` calls
    would, so a point outside the debris gate is not also counted as failing the singlets
    gate: it is counted once, at the gate that removed it.
    """

    table.validate(channels)
    index_of = {name: i for i, name in enumerate(channels)}
    missing_values = [c for c in table.channels if c not in index_of]
    if missing_values:  # pragma: no cover - validate() already refuses this
        raise GatingError(f"the data does not carry gate channel(s) {missing_values}")

    n = int(values.shape[0])
    alive = np.ones(n, dtype=bool)
    outcomes: list[GateOutcome] = []
    rows: list[pd.DataFrame] = []

    for gate in table.gates:
        x = values[:, index_of[gate.mk1]]
        y = None if gate.mk2 is None else values[:, index_of[gate.mk2]]
        keep = gate.contains(x, y) & alive
        removed_mask = alive & ~keep
        n_removed = int(removed_mask.sum())
        if n_removed:
            rows.append(
                pd.DataFrame(
                    {
                        "sample_id": np.asarray(sample_ids)[removed_mask],
                        "original_event_index": np.asarray(original_event_index)[removed_mask],
                        "step": step,
                        "reason": gate.reason,
                    }
                )
            )
        per_sample = {
            str(sid): int(c)
            for sid, c in zip(
                *np.unique(np.asarray(sample_ids)[removed_mask], return_counts=True),
                strict=True,
            )
        }
        outcomes.append(
            GateOutcome(
                gate=gate,
                kept=int(keep.sum()),
                removed=n_removed,
                n_input=int(alive.sum()),
                removed_by_sample=per_sample,
            )
        )
        alive = keep

    removed = (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=["sample_id", "original_event_index", "step", "reason"])
    )
    return GatingResult(keep=alive, outcomes=tuple(outcomes), removed=removed)
