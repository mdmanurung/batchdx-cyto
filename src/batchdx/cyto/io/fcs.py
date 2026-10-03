"""FCS reading and writing.

Source FCS files are opened read-only and never written back to.  Export writes new
files, records the channel order and the scale of the exported values, and stamps a
keyword policy block so that a downstream tool cannot silently compensate or
transform the same values twice.

One property of the format is load-bearing and is stated rather than assumed: an FCS
event block is IEEE-754 *single* precision, so a write/read round trip is lossy at
the fourth significant decimal.  :func:`assert_fcs_round_trip` measures that loss
against the bound the format actually guarantees, so a caller that needs bit-exact
values is told rather than discovering it later.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:  # pragma: no cover - import guard
    import flowio
except ImportError as exc:  # pragma: no cover
    raise ImportError("reading and writing FCS requires the 'flowio' package") from exc

#: Keyword written into every exported file.  A reader that honours it knows the
#: values are already processed and must not be compensated or transformed again.
PROCESSING_KEYWORD = "BATCHDX_PROCESSING"
SCALE_KEYWORD = "BATCHDX_SCALE"
SOURCE_KEYWORD = "BATCHDX_SOURCE_SAMPLE"
CANDIDATE_KEYWORD = "BATCHDX_CANDIDATE"

#: FCS event values are single precision.  A round trip is therefore exact only to
#: about seven significant decimal digits, and comparing float64 against a file
#: expecting bit equality fails on values that were never wrong.
FCS_STORAGE_DTYPE = np.float32


class FCSError(ValueError):
    """Raised for an unreadable FCS file or an unusable channel layout."""


@dataclass(frozen=True)
class FCSContents:
    """One FCS file as read: events, channel identity, and the original keywords."""

    path: Path
    events: np.ndarray  # (n_events, n_channels), float64
    channel_names: list[str]  # $PnN
    channel_labels: list[str]  # $PnS, empty string where absent
    keywords: dict[str, str]

    @property
    def n_events(self) -> int:
        return int(self.events.shape[0])


def read_fcs(path: str | Path) -> FCSContents:
    """Read one FCS file.  The file is opened read-only and left untouched."""

    p = Path(path)
    if not p.exists():
        raise FCSError(f"FCS file not found: {p}")
    try:
        fd = flowio.FlowData(str(p))
    except Exception as exc:  # flowio raises several unrelated exception types
        raise FCSError(f"{p}: could not be read as FCS: {exc}") from exc

    n_channels = int(fd.channel_count)
    events = np.reshape(np.asarray(fd.events, dtype=np.float64), (-1, n_channels))

    names: list[str] = []
    labels: list[str] = []
    for i in range(1, n_channels + 1):
        # flowio keys channels by integer and lower-cases the FCS keyword names.
        meta = fd.channels.get(i) or fd.channels.get(str(i)) or {}
        names.append(str(meta.get("pnn", meta.get("PnN", f"P{i}"))))
        labels.append(str(meta.get("pns", meta.get("PnS", "")) or ""))

    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        raise FCSError(
            f"{p}: channel name(s) {dupes} appear more than once; "
            "the channel to feature mapping would be ambiguous"
        )

    keywords = {str(k): str(v) for k, v in fd.text.items()}
    return FCSContents(path=p, events=events, channel_names=names,
                       channel_labels=labels, keywords=keywords)


def write_fcs(
    path: str | Path,
    events: np.ndarray,
    channel_names: list[str],
    channel_labels: list[str] | None = None,
    *,
    scale: str,
    processing_note: str,
    candidate_id: str | None = None,
    source_sample: str | None = None,
    extra_keywords: dict[str, str] | None = None,
    overwrite: bool = False,
) -> Path:
    """Write a new FCS file with an explicit scale and processing keyword block.

    ``scale`` and ``processing_note`` are mandatory: an exported file always states
    what its numbers are, so values that could not be safely inverted back to
    instrument units are never presented as instrument measurements.
    """

    p = Path(path)
    if p.exists() and not overwrite:
        raise FCSError(f"refusing to overwrite existing FCS file {p}")
    if events.ndim != 2:
        raise FCSError(f"event array must be 2-dimensional, got shape {events.shape}")
    if events.shape[1] != len(channel_names):
        raise FCSError(
            f"event array has {events.shape[1]} columns but {len(channel_names)} channel names"
        )
    if not np.isfinite(events).all():
        raise FCSError("refusing to export non-finite values to FCS")

    metadata = {
        PROCESSING_KEYWORD: processing_note,
        SCALE_KEYWORD: scale,
    }
    if candidate_id:
        metadata[CANDIDATE_KEYWORD] = candidate_id
    if source_sample:
        metadata[SOURCE_KEYWORD] = source_sample
    if extra_keywords:
        metadata.update({str(k): str(v) for k, v in extra_keywords.items()})

    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "wb") as fh:
        flowio.create_fcs(
            fh,
            np.asarray(events, dtype=np.float64).flatten().tolist(),
            channel_names=list(channel_names),
            opt_channel_names=list(channel_labels) if channel_labels else None,
            metadata_dict=metadata,
        )
    return p


def round_trip_tolerance(events: np.ndarray) -> float:
    """The largest absolute error an FCS round trip can introduce for ``events``.

    Single precision carries about seven significant decimal digits, so the
    representable spacing at a value ``v`` is ``2**-23 * 2**floor(log2|v|)``.  The
    bound is one spacing, plus the smallest subnormal spacing so that an all-zero
    matrix has a finite, non-zero bound rather than a division by zero.
    """

    arr = np.asarray(events, dtype=np.float64)
    if arr.size == 0:
        return 0.0
    magnitude = np.max(np.abs(arr))
    if magnitude == 0.0:
        return float(np.finfo(FCS_STORAGE_DTYPE).smallest_subnormal)
    return float(np.spacing(np.float32(magnitude)))


def assert_fcs_round_trip(
    events: np.ndarray,
    contents: FCSContents,
    channel_names: list[str] | None = None,
    *,
    context: str = "",
) -> float:
    """Check that ``contents`` is the file ``events`` was written to, and by how much.

    Channel order and event count must match exactly.  Values are compared against
    the bound the format itself guarantees, so a genuine ordering mix-up fails while
    the unavoidable single-precision loss does not.  Returns the observed maximum
    absolute error, so a caller can record the real number rather than asserting
    equality it did not achieve.
    """

    expected = list(channel_names) if channel_names is not None else list(contents.channel_names)
    if list(contents.channel_names) != expected:
        raise FCSError(
            f"{contents.path}: channel order is {contents.channel_names}, expected {expected}"
        )
    original = np.asarray(events, dtype=np.float64)
    if contents.events.shape != original.shape:
        raise FCSError(
            f"{contents.path}: round trip produced {contents.events.shape}, "
            f"expected {original.shape}"
        )
    error = float(np.max(np.abs(contents.events - original))) if original.size else 0.0
    tolerance = round_trip_tolerance(original)
    if not np.isfinite(error) or error > tolerance:
        where = f" in {contents.path}" if contents.path else ""
        raise FCSError(
            f"FCS round trip{where}{(' for ' + context) if context else ''} changed values by "
            f"{error:.6g}, beyond the {tolerance:.6g} that single-precision FCS storage "
            f"allows; the events were not preserved"
        )
    return error
