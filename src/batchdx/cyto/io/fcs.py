"""FCS reading and writing.

Source FCS files are opened read-only and never written back to.  Export writes new
files, records the channel order and the scale of the exported values, and stamps a
keyword policy block so that a downstream tool cannot silently compensate or
transform the same values twice.
"""

from __future__ import annotations

import re
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

#: How the FCS Time channel is recognised.  The FCS standard does not fix the name;
#: instruments in practice use $PnN/$PnS values beginning with "Time".
TIME_CHANNEL_PATTERN = re.compile(r"^time", re.IGNORECASE)


class FCSError(ValueError):
    """Raised for an unreadable FCS file or an unusable channel layout."""


@dataclass(frozen=True)
class FCSContents:
    """One FCS file as read: events, channel identity, and the original keywords."""

    path: Path
    events: np.ndarray  # (n_events, n_channels), float64
    channel_names: list[str]  # resolved channel names (see ``channel_key``)
    channel_labels: list[str]  # $PnS, empty string where absent
    keywords: dict[str, str]
    #: Which keyword the resolved channel names came from: "$PnN", "$PnS", or
    #: "positional" when both keyword columns were duplicate or empty.
    channel_key: str = "$PnN"
    #: Name of the FCS Time channel and its per-event values, when the file has one.
    time_channel: str | None = None
    time_values: np.ndarray | None = None

    @property
    def n_events(self) -> int:
        return int(self.events.shape[0])


def _usable(values: list[str]) -> bool:
    """A keyword column is usable when every value is non-empty and unique."""
    return bool(values) and all(v.strip() for v in values) and len(set(values)) == len(values)


def _resolve_channel_names(names: list[str], labels: list[str]) -> tuple[list[str], str]:
    """Resolve channel identity: $PnN, then $PnS, then positional mapping.

    Duplicate or empty $PnN values used to be a hard error; they are now resolved by
    falling back to $PnS and then to the channel position, and the key that was used
    is recorded on :class:`FCSContents`.  Positional mapping is unambiguous by
    construction, so no genuine ambiguity remains once it is reached.
    """
    if _usable(names):
        return names, "$PnN"
    if _usable(labels):
        return labels, "$PnS"
    return [f"P{i + 1}" for i in range(len(names))], "positional"


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
    if n_channels <= 0:
        raise FCSError(f"{p}: the file declares {n_channels} channels")
    events = np.reshape(np.asarray(fd.events, dtype=np.float64), (-1, n_channels))

    names_raw: list[str] = []
    labels_raw: list[str] = []
    for i in range(1, n_channels + 1):
        # flowio keys channels by integer and lower-cases the FCS keyword names.
        meta = fd.channels.get(i) or fd.channels.get(str(i)) or {}
        names_raw.append(str(meta.get("pnn", meta.get("PnN", "")) or ""))
        labels_raw.append(str(meta.get("pns", meta.get("PnS", "")) or ""))

    names, channel_key = _resolve_channel_names(names_raw, labels_raw)

    time_idx = next(
        (i for i in range(n_channels)
         if TIME_CHANNEL_PATTERN.search(names_raw[i])
         or TIME_CHANNEL_PATTERN.search(labels_raw[i])),
        None,
    )
    time_channel = None
    time_values = None
    if time_idx is not None:
        time_channel = names_raw[time_idx] or labels_raw[time_idx] or names[time_idx]
        time_values = events[:, time_idx].copy()

    keywords = {str(k): str(v) for k, v in fd.text.items()}
    return FCSContents(path=p, events=events, channel_names=names,
                       channel_labels=labels_raw, keywords=keywords,
                       channel_key=channel_key,
                       time_channel=time_channel, time_values=time_values)


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
