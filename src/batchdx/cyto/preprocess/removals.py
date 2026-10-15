"""The per-event removal log, and the one writer it lacked.

``Dataset.removed_events`` (``io/dataset.py:54-56``) declares exactly the four columns
this log needs and, until now, had no writer and no reader anywhere in the repository:
the field was a promise with nothing behind it.  This module is the writer, and
:func:`read_removed_events` / :func:`write_removed_events` in
:mod:`batchdx.cyto.io.manifests` are the reader and the serialisation.

FOUR COLUMNS, AND WHAT EACH ONE IS FOR
--------------------------------------
``sample_id`` / ``original_event_index``
    Which event, in the *source* file's own numbering, before anything was removed.
    ``original_event_index`` is the identity channel carried through every write, so a
    removed event can be found in the original FCS and not merely in a derived frame.
``step``
    Which stage removed it, so a removal by the preliminary gate is never confused with
    a removal by PeacoQC.  Both are named constants here rather than free strings.
``reason``
    The criterion, in the gate's or the QC step's own terms.  Not a boolean and not a
    count: "removed by the singlets polygon over FSC-A/FSC-H" and "removed by PeacoQC
    under determine_good_cells='all'" are different claims, and a log that collapsed
    them would say only that the event is absent.

WHAT THIS DOES NOT CLAIM
------------------------
One event removed at two steps appears once, at the step that removed it.  That is not
a convenience: a chained gate is evaluated on the survivors, so an event outside the
debris gate is never offered to the singlets gate and cannot fail it.  Counting it twice
would inflate the removal total above the number of events actually lost.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

#: The columns, in order, exactly as ``Dataset.removed_events`` declares them.
REMOVED_EVENT_COLUMNS = ["sample_id", "original_event_index", "step", "reason"]

#: The step name for a removal by a preliminary gate.
STEP_GATE = "prelim_gating"

#: The step name for a removal by PeacoQC.
STEP_QC = "peacoqc"

#: Every step this build can name.  A log carrying anything else is refused rather than
#: accepted with an unrecognised label, because a step nobody can interpret is how a
#: removal log becomes decorative.
KNOWN_STEPS = (STEP_GATE, STEP_QC)


class RemovalLogError(ValueError):
    """Raised when a removal log is malformed."""


def empty_removal_log() -> pd.DataFrame:
    """A log with the right columns and no rows."""

    return pd.DataFrame(columns=REMOVED_EVENT_COLUMNS)


def removal_frame(
    *,
    sample_id: np.ndarray | list,
    original_event_index: np.ndarray | list,
    step: str,
    reason: str,
) -> pd.DataFrame:
    """Build one step's removal rows, with the dtypes the log declares.

    ``original_event_index`` is cast to a signed integer and ``sample_id`` to a string,
    because a log read back from CSV would otherwise carry them as ``float64`` and
    ``object`` respectively and every join against ``Dataset.obs`` would silently
    produce an empty frame rather than an error.
    """

    if step not in KNOWN_STEPS:
        raise RemovalLogError(
            f"removal step {step!r} is not one this build names; expected one of "
            f"{list(KNOWN_STEPS)}. An uninterpretable step label is worse than none"
        )
    if not str(reason).strip():
        raise RemovalLogError(
            "a removal without a reason is not a record: say which criterion removed the "
            "event. 'it is missing from the output' is not a reason"
        )
    sample = np.asarray(sample_id).astype(str)
    index = np.asarray(original_event_index)
    if index.dtype.kind == "f":
        if not np.isfinite(index).all():
            raise RemovalLogError("original_event_index holds a non-finite value")
        index = np.rint(index)
    index = index.astype(np.int64)
    if sample.shape != index.shape:
        raise RemovalLogError(
            f"sample_id has {sample.shape[0]} entries and original_event_index has "
            f"{index.shape[0]}; the two must be row-aligned"
        )
    return pd.DataFrame(
        {
            "sample_id": sample,
            "original_event_index": index,
            "step": np.repeat(step, sample.shape[0]),
            "reason": np.repeat(str(reason), sample.shape[0]),
        },
        columns=REMOVED_EVENT_COLUMNS,
    )


def check_removal_log(frame: pd.DataFrame, *, what: str = "removal log") -> pd.DataFrame:
    """Validate a log's columns, step names and dtypes, and return it normalised.

    Returns a copy with the dtypes the rest of the package assumes, so that a log read
    from a file and one built in memory behave identically.
    """

    if not isinstance(frame, pd.DataFrame):
        raise RemovalLogError(f"{what}: expected a DataFrame, got {type(frame).__name__}")
    missing = [c for c in REMOVED_EVENT_COLUMNS if c not in frame.columns]
    if missing:
        raise RemovalLogError(
            f"{what}: missing column(s) {missing}; expected {REMOVED_EVENT_COLUMNS}"
        )
    extra = [c for c in frame.columns if c not in REMOVED_EVENT_COLUMNS]
    if extra:
        raise RemovalLogError(
            f"{what}: unexpected column(s) {extra}. A removal log that has grown columns "
            "nothing renders is a log that has started accumulating silently"
        )
    out = frame.loc[:, REMOVED_EVENT_COLUMNS].copy()
    out["sample_id"] = out["sample_id"].astype(str)
    out["original_event_index"] = (
        out["original_event_index"].astype(float).round().astype(np.int64)
    )
    unknown = sorted(set(out["step"]) - set(KNOWN_STEPS))
    if unknown:
        raise RemovalLogError(
            f"{what}: step label(s) {unknown} are not among {list(KNOWN_STEPS)}"
        )
    if (out["reason"].astype(str).str.strip() == "").any():
        raise RemovalLogError(f"{what}: at least one row has an empty reason")
    return out.reset_index(drop=True)


def summarise_removals(frame: pd.DataFrame) -> pd.DataFrame:
    """Removed-event counts per step, per reason and per sample.

    This is the *reported* shape: a reader wants "how many events did each criterion
    remove, where", not the log itself.  It is computed from the log, never assembled
    separately, so the two cannot disagree.
    """

    checked = check_removal_log(frame)
    if checked.empty:
        return pd.DataFrame(columns=["step", "reason", "sample_id", "n_removed"])
    grouped = (
        checked.groupby(["step", "reason", "sample_id"], dropna=False)
        .size()
        .reset_index(name="n_removed")
    )
    return grouped.sort_values(["step", "reason", "sample_id"]).reset_index(drop=True)
