"""PeacoQC event-level quality control, called from the installed Bioconductor package.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
The reference implementation runs ``PeacoQC`` from ``QC.R:102-144`` and then **discards
the removal index**: it renames the report file and moves on.  This module exists to keep
that index, so every event PeacoQC removes is attributable to a criterion rather than
merely absent from the output.

THE ENTRY POINT, READ FROM THE INSTALLED PACKAGE
------------------------------------------------
The call shape is taken from the reference call site, which is evidence of the formals
*it* passes and of nothing else.  The full signature below was read from
``formals(PeacoQC)`` of **PeacoQC 1.16.0** (the version
``.upstream/spiscy/docker/rbase.yaml:20`` pins), and every default this adapter relies on
is recorded in :data:`VERIFIED_FORMALS` so a change upstream cannot pass unnoticed:

.. code-block:: text

    PeacoQC(ff, channels, determine_good_cells = "all", plot = 20, save_fcs = TRUE,
            output_directory = ".", name_directory = "PeacoQC_results", report = TRUE,
            events_per_bin = FindEventsPerBin(...), min_cells = 150, max_bins = 500,
            step = 500, MAD = 6, IT_limit = 0.6, consecutive_bins = 5,
            remove_zeros = FALSE, suffix_fcs = "_QC", force_IT = 150,
            peak_removal = (1/3), min_nr_bins_peakdetection = 10,
            time_channel_parameter = "Time", ...)

Two of those defaults differ from what the reference call site passes, and both are
passed explicitly here so the choice is this adapter's and is recorded:
``save_fcs`` is ``TRUE`` upstream and ``FALSE`` in SPiSCy, and ``plot`` is the *number*
``20`` upstream (plot once at least 20% of measurements are removed) while SPiSCy passes
``TRUE``.

THERE IS NO ``$BadCell``
------------------------
PeacoQC 1.16.0's return value has no ``BadCell`` element.  Its ``names()`` are
``Analysis, Channels, nr_bins, EventsPerBin, WeirdChannels, <one per channel>,
ITPercentage, OutlierMads, ContributionMad, MADPercentage, ConsecutiveCells,
ConsecutiveCellsPercentage, GoodCells, PercentageRemoved, FinalFF``.  The removed set is
therefore ``!GoodCells``: a logical vector aligned to the input event order, which is
what this adapter maps back to the baseline's own event indices.  The fact that the
element is named ``GoodCells`` and not ``BadCell`` is checked at run time rather than
assumed, because a package update that renamed it would otherwise produce a silent
"nothing was removed" result.

INSUFFICIENT DATA IS A STATUS, NOT A ZERO
-----------------------------------------
``PeacoQC`` sets ``ITPercentage`` to ``NA`` when it declines to run the isolation tree
("not enough bins for a robust isolation tree analysis") and then removes on MAD alone.
So with ``determine_good_cells = "all"`` a run can complete, remove events, and have
done **no** IT analysis at all.  That is reported as ``it_available: False`` with the
reason, and the criteria that actually contributed are named, so a MAD-only pass can
never be read as the "all" that was requested.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from batchdx.cyto.adapters.r_runner import (
    RNotAvailable,
    RScriptFailed,
    r_available,
    r_package_version,
    r_version,
    run_r_script,
)
from batchdx.cyto.io.fcs import read_fcs, write_fcs
from batchdx.cyto.io.hashing import file_sha256
from batchdx.cyto.preprocess.removals import STEP_QC, removal_frame
from batchdx.cyto.schemas.contracts import SoftwareVersion

R_SCRIPT = "peacoqc_qc.R"

#: The channel carrying each event's original index through the FCS round trip, so the
#: removal index PeacoQC returns can be resolved back to the baseline.  Named so it
#: cannot collide with a real parameter by accident, and refused by the caller if it
#: does.
IDENTITY_CHANNEL = "BATCHDX_EVENT_IDX"

#: The formals of ``PeacoQC`` 1.16.0 as read from the installed package, recorded so a
#: silent upstream change is visible in a diff rather than in a wrong result.
VERIFIED_FORMALS: dict[str, Any] = {
    "ff": None,
    "channels": None,
    "determine_good_cells": "all",
    "plot": 20,
    "save_fcs": True,
    "output_directory": ".",
    "name_directory": "PeacoQC_results",
    "report": True,
    "min_cells": 150,
    "max_bins": 500,
    "step": 500,
    "MAD": 6,
    "IT_limit": 0.6,
    "consecutive_bins": 5,
    "remove_zeros": False,
    "suffix_fcs": "_QC",
    "force_IT": 150,
    "peak_removal": 1.0 / 3.0,
    "min_nr_bins_peakdetection": 10,
    "time_channel_parameter": "Time",
}

#: The values ``determine_good_cells`` accepts, per the reference call site's own
#: validation at ``QC.R:96-144``.
CRITERIA = ("all", "IT", "MAD")

#: The package version this module's constants were read from.  Recorded because a
#: formals table for a different version is a guess, not a reading.
VERIFIED_VERSION = "1.16.0"


class PeacoQCError(RuntimeError):
    """Raised when PeacoQC ran and did not produce a usable removal index."""


class PeacoQCNotAvailable(RuntimeError):
    """Raised when PeacoQC cannot run here.  A stated refusal, not a failure."""


@dataclass(frozen=True)
class PeacoQCSettings:
    """The declared PeacoQC settings, all of which are read by the R script."""

    channels: tuple[str, ...]
    criteria: str = "all"
    mad: float = 6.0
    it_limit: float = 0.6
    min_cells: int = 150
    max_bins: int = 500
    step: int = 500
    consecutive_bins: int = 5
    remove_zeros: bool = False
    force_IT: int = 150
    peak_removal: float = 1.0 / 3.0
    min_nr_bins_peakdetection: int = 10
    time_channel: str | None = None
    plot: bool | str | int | float = False
    report: bool = True

    def __post_init__(self) -> None:
        # YAML and JSON both hand back a list for a sequence and a bare string for a
        # scalar, so ``channels: CD4`` and ``channels: [CD4]`` are the same declaration.
        # Coerced here rather than at the call site, because a reader that returned a
        # one-element list would compare unequal to a one-element tuple everywhere
        # downstream and the difference is a spelling, not a setting.
        if isinstance(self.channels, str):
            coerced: tuple[str, ...] = (self.channels,)
        else:
            coerced = tuple(str(c) for c in self.channels)
        object.__setattr__(self, "channels", coerced)
        if not self.channels:
            raise PeacoQCError(
                "PeacoQC needs at least one channel to evaluate; a QC call over an empty "
                "channel list would remove nothing and report success"
            )
        if self.criteria not in CRITERIA:
            raise PeacoQCError(
                f"determine_good_cells={self.criteria!r} is not one of {list(CRITERIA)}; "
                "the reference implementation refuses anything else (QC.R:144)"
            )
        if self.mad <= 0 or self.it_limit <= 0:
            raise PeacoQCError(
                f"MAD={self.mad} and IT_limit={self.it_limit} must both be positive; a "
                "non-positive threshold would flag every measurement or none, and neither "
                "is a quality criterion"
            )
        if self.step <= 0 or self.min_cells <= 0 or self.max_bins <= 0:
            raise PeacoQCError(
                "min_cells, max_bins and step are counts and must be positive; "
                f"got min_cells={self.min_cells}, max_bins={self.max_bins}, step={self.step}"
            )
        if self.criteria in {"all", "MAD"} and self.mad is None:  # pragma: no cover
            raise PeacoQCError("criteria requires a MAD value")

    @property
    def reason(self) -> str:
        """The criterion named in the removal log, naming only what was actually asked."""

        parts = [f"determine_good_cells={self.criteria}"]
        if self.criteria in {"all", "MAD"}:
            parts.append(f"MAD={self.mad:g}")
        if self.criteria in {"all", "IT"}:
            parts.append(f"IT_limit={self.it_limit:g}")
        return f"PeacoQC{VERIFIED_VERSION}:" + ",".join(parts)

    def to_params(self, output_dir: str | Path) -> dict[str, Any]:
        return {
            "channels": list(self.channels),
            "criteria": self.criteria,
            "mad": float(self.mad),
            "it_limit": float(self.it_limit),
            "min_cells": int(self.min_cells),
            "max_bins": int(self.max_bins),
            "step": int(self.step),
            "consecutive_bins": int(self.consecutive_bins),
            "remove_zeros": bool(self.remove_zeros),
            "force_IT": int(self.force_IT),
            "peak_removal": float(self.peak_removal),
            "min_nr_bins_peakdetection": int(self.min_nr_bins_peakdetection),
            "time_channel_parameter": self.time_channel,
            "plot": self.plot,
            "report": bool(self.report),
            "output_dir": str(output_dir),
        }


#: Keys the **reference implementation's own settings file** uses, mapped onto the
#: field names here.  ``.upstream/spiscy/config/06_QC.yaml`` names its thresholds
#: ``MAD`` and ``IT`` and its channel list ``markers_to_evaluate``; none of those is a
#: formal of ``PeacoQC``'s.  They are accepted as aliases so that file loads unchanged,
#: and each one is read -- an alias is a spelling, not a key that is tolerated and
#: dropped.  Anything outside this map and the dataclass's own fields is refused.
SETTINGS_ALIASES: dict[str, str] = {
    "MAD": "mad",
    "IT": "it_limit",
    "IT_limit": "it_limit",
    "markers_to_evaluate": "channels",
    "time_channel_parameter": "time_channel",
}


def read_peacoqc_settings(path: str | Path) -> PeacoQCSettings:
    """Read a PeacoQC settings file (YAML or JSON).

    ``channels`` is required and is never inferred: PeacoQC's quality assessment is
    per-channel, and choosing the channels for someone is how a QC pass silently
    evaluates the wrong panel.  The reference implementation's own spellings
    (:data:`SETTINGS_ALIASES`) are accepted so its ``06_QC.yaml`` loads as written.
    """

    p = Path(path)
    if not p.exists():
        raise PeacoQCError(f"PeacoQC settings file not found: {p}")
    text = p.read_text(encoding="utf-8")
    raw = (
        json.loads(text) if p.suffix.lower() == ".json" else yaml_safe_load(text, p)
    )
    if not isinstance(raw, dict):
        raise PeacoQCError(f"{p}: top level must be a mapping")
    renamed: dict[str, Any] = {}
    for key, value in raw.items():
        name = SETTINGS_ALIASES.get(str(key), str(key))
        if name in renamed:
            raise PeacoQCError(
                f"{p}: setting {key!r} is a second spelling of {name!r}, which is also "
                "set. Two spellings of one setting can disagree, and the one that lost "
                "would be the one acted on"
            )
        renamed[name] = value
    if "channels" not in renamed:
        raise PeacoQCError(
            f"{p}: 'channels' is required and is never inferred. PeacoQC evaluates "
            "quality per channel, so an omitted list would mean the wrong panel"
        )
    known = set(PeacoQCSettings.__dataclass_fields__)
    unknown = sorted(set(renamed) - known)
    if unknown:
        raise PeacoQCError(
            f"{p}: unknown setting(s) {unknown}; expected a subset of {sorted(known)}. A "
            "setting this build does not read would be accepted and then ignored"
        )
    return PeacoQCSettings(**renamed)


def settings_from_config(config, base_dir: str | Path = ".") -> PeacoQCSettings | None:
    """The settings ``preprocessing.peacoqc`` declares, or ``None`` if it declares none.

    This is the reader that makes ``preprocessing.peacoqc`` a key rather than a promise.
    It reads the key and nothing else: the channels come from the settings file and are
    never taken from the panel, because PeacoQC's assessment is per-channel and choosing
    the channels for someone is how a QC pass silently evaluates the wrong panel.
    """

    declared = config.preprocessing.peacoqc
    if declared is None:
        return None
    return read_peacoqc_settings(Path(base_dir) / str(declared))


def yaml_safe_load(text: str, p: Path) -> dict:
    import yaml

    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise PeacoQCError(f"{p}: not valid YAML: {exc}") from exc


def peacoqc_availability(rscript: str = "Rscript") -> list[SoftwareVersion]:
    """R and PeacoQC as they actually are, without assuming either."""

    return [
        r_version(rscript),
        r_package_version("PeacoQC", rscript),
        r_package_version("flowCore", rscript),
    ]


# --------------------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SampleQCOutcome:
    """What PeacoQC did to one sample, including the parts it could not do."""

    sample_id: str
    n_input: int
    n_removed: int
    removed_fraction: float | None
    it_fraction: float | None
    mad_fraction: float | None
    criteria_applied: tuple[str, ...]
    it_available: bool
    it_reason: str | None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class PeacoQCResult:
    """The removal log, the per-sample outcomes, and the software that produced them.

    ``output_files`` is what PeacoQC actually wrote, counted and named by the script
    rather than predicted.  ``report = TRUE`` is not a promise PeacoQC keeps on every
    route -- it writes no report when ``determine_good_cells = "MAD"`` -- so the
    artifact is reported as it is, and a settings file asking for a report on a route
    that produces none is visible instead of silently unfulfilled.
    """

    removed: pd.DataFrame
    samples: tuple[SampleQCOutcome, ...]
    software: tuple[SoftwareVersion, ...] = ()
    warnings: tuple[str, ...] = ()
    status: str = "ok"
    reason: str | None = None
    report_paths: dict[str, str] = field(default_factory=dict)
    output_files: dict[str, list[str]] = field(default_factory=dict)

    @property
    def n_removed(self) -> int:
        return int(len(self.removed))

    def per_sample(self) -> pd.DataFrame:
        rows = [
            {
                "sample_id": s.sample_id,
                "n_input": s.n_input,
                "n_removed": s.n_removed,
                "removed_fraction": s.removed_fraction,
                "it_fraction": s.it_fraction,
                "mad_fraction": s.mad_fraction,
                "criteria_applied": "/".join(s.criteria_applied),
                "it_available": s.it_available,
                "it_reason": s.it_reason or "",
            }
            for s in self.samples
        ]
        return pd.DataFrame(rows)


def _as_list(value: Any) -> list[Any]:
    """``jsonlite`` with ``auto_unbox = TRUE`` collapses a length-one vector to a scalar."""

    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [value]
    if isinstance(value, dict):
        return [value]
    return list(value)


def _as_float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(out) else out


def run_peacoqc(
    values: np.ndarray,
    channels: list[str] | tuple[str, ...],
    *,
    sample_ids: np.ndarray,
    original_event_index: np.ndarray,
    settings: PeacoQCSettings,
    work_dir: str | Path,
    rscript: str = "Rscript",
    timeout: float | None = None,
) -> PeacoQCResult:
    """Run PeacoQC per sample and return every event it removed, with the reason.

    One FCS file per sample is written, carrying :data:`IDENTITY_CHANNEL`, because
    PeacoQC's audited interface is a ``flowFrame`` and a flowFrame is what
    ``flowCore::read.FCS`` returns.  The file is read back to resolve PeacoQC's row
    indices into the caller's own event indices, so nothing here depends on row order
    surviving by luck.

    Raises :class:`PeacoQCNotAvailable` when R or PeacoQC is missing, and
    :class:`PeacoQCError` when PeacoQC ran but the result cannot be trusted.  Neither is
    ever reported as "no events were removed".
    """

    # The call is validated BEFORE anything is looked up on PATH.  A malformed call is
    # wrong in every environment; reporting "no R on PATH" for one would hide the defect
    # behind an environment fact, and the no-R path is already covered by
    # rscript="Rscript-does-not-exist".
    if IDENTITY_CHANNEL in {str(c).upper() for c in channels}:
        raise PeacoQCError(
            f"the channel {IDENTITY_CHANNEL!r} is reserved for the event-identity check "
            "and collides with a channel PeacoQC was asked to evaluate; rename it"
        )
    if not settings.channels:
        raise PeacoQCError(
            "no channel was named for PeacoQC to evaluate; a QC call over an empty "
            "channel list removes nothing and reports success"
        )
    missing = [c for c in settings.channels if c not in set(channels)]
    if missing:
        raise PeacoQCError(
            f"PeacoQC was asked to evaluate channel(s) {missing}, which the data does not "
            "carry. A missing channel is not a quiet one: PeacoQC would evaluate whatever "
            "was closest or fail on a channel mismatch"
        )

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    in_dir = work / "fcs_in"
    samples = np.asarray(sample_ids).astype(str)
    event_index = np.asarray(original_event_index).astype(np.int64)
    if samples.shape[0] != values.shape[0] or event_index.shape[0] != values.shape[0]:
        raise PeacoQCError(
            f"values has {values.shape[0]} rows but sample_ids has {samples.shape[0]} and "
            f"original_event_index has {event_index.shape[0]}; the three must be row-aligned"
        )

    if not r_available(rscript):
        raise PeacoQCNotAvailable(
            f"{rscript!r} is not on PATH, so PeacoQC cannot run. This build does not "
            "substitute another quality-control algorithm for it"
        )
    probe = r_package_version("PeacoQC", rscript)
    if not probe.available:
        raise PeacoQCNotAvailable(
            "the PeacoQC R package is not installed in this environment "
            f"({probe.detail}). Install bioconductor-peacoqc; scripts/setup_r_env.sh "
            "lists it"
        )

    file_of_sample: dict[str, str] = {}
    for sample_id in sorted(set(samples)):
        rows = np.flatnonzero(samples == sample_id)
        if rows.size == 0:  # pragma: no cover - set() came from samples
            continue
        ordered = rows[np.argsort(event_index[rows], kind="stable")]
        index = event_index[ordered].astype(np.float64)
        block = np.column_stack([index, np.asarray(values[ordered], dtype=np.float64)])
        path = write_fcs(
            in_dir / f"{sample_id}.fcs",
            block,
            [IDENTITY_CHANNEL, *[str(c) for c in channels]],
            [IDENTITY_CHANNEL, *[str(c) for c in channels]],
            scale="raw instrument units, as handed to the QC step",
            processing_note=(
                "batchdx-cyto per-sample file for PeacoQC; values are exactly as loaded, "
                f"carrying the original event index in {IDENTITY_CHANNEL} so the removal "
                "index can be resolved back to the baseline"
            ),
            source_sample=sample_id,
            extra_keywords={"BATCHDX_PURPOSE": "peacoqc_input"},
            overwrite=True,
        )
        file_of_sample[sample_id] = str(path)

    params = settings.to_params(work / "peacoqc_out")
    ordered = sorted(file_of_sample)
    # The R script keys its results by SAMPLE ID, sent explicitly rather than derived
    # from a file name. Deriving it would make the two sides disagree by an extension on
    # the one field the removal log is built from. Two parallel character vectors rather
    # than a list of records, because jsonlite's simplifyVector collapses a list of
    # objects into a data frame.
    params["sample_ids"] = list(ordered)
    params["files"] = [file_of_sample[sid] for sid in ordered]
    try:
        r_result = run_r_script(
            R_SCRIPT, params, work, rscript=rscript,
            log_path=work / "peacoqc.log", timeout=timeout,
        )
    except RNotAvailable as exc:  # pragma: no cover - checked above
        raise PeacoQCNotAvailable(str(exc)) from exc
    except RScriptFailed as exc:
        raise PeacoQCError(
            f"the PeacoQC script failed: {exc}. Last lines of the R log:\n"
            f"{exc.log[-2000:]}"
        ) from exc

    payload = r_result.payload
    if payload.get("status") != "ok":
        raise PeacoQCError(
            f"the PeacoQC script did not report success: {payload.get('status')!r} "
            f"{payload.get('reason', '')}".strip()
        )
    per_file = payload.get("samples")
    if not isinstance(per_file, dict) or not per_file:
        raise PeacoQCError(
            "the PeacoQC script returned no per-sample result; an empty result is not a "
            "result in which nothing was removed"
        )
    absent = sorted(set(file_of_sample) - set(per_file))
    if absent:
        raise PeacoQCError(
            f"the PeacoQC script returned no entry for sample(s) {absent}; a sample that "
            "was not assessed is not a sample that passed"
        )

    reason = settings.reason
    frames: list[pd.DataFrame] = []
    outcomes: list[SampleQCOutcome] = []
    for sample_id, path in file_of_sample.items():
        entry = per_file[sample_id]
        n_input = int(entry.get("n_input", 0))
        good = _as_list(entry.get("good_cells"))
        if len(good) != n_input:
            raise PeacoQCError(
                f"PeacoQC returned {len(good)} per-cell decisions for sample "
                f"{sample_id!r}, which had {n_input} events. A decision vector of the "
                "wrong length cannot be mapped to events and is refused rather than "
                "truncated or padded"
            )
        if "GoodCells" not in entry.get("result_keys", []):
            raise PeacoQCError(
                f"PeacoQC's result for sample {sample_id!r} has no 'GoodCells' element; "
                f"it has {entry.get('result_keys')}. The removal index is derived from "
                "!GoodCells, and on a version without it nothing here can say which events "
                "were removed"
            )
        keep = np.array([bool(g) for g in good], dtype=bool)
        removed_rows = np.flatnonzero(~keep)
        # One read of the file serves both purposes: it proves the identity channel
        # survived the round trip, and it is what the removal positions are resolved
        # against. Skipping the read when nothing was removed would leave the "PeacoQC
        # looked at the file we think it looked at" claim untested in exactly the case a
        # reader would want it most (a clean sample).
        written_index = _read_identity_channel(path, sample_id)
        removed_index = (
            np.zeros(0, dtype=np.int64)
            if removed_rows.size == 0
            else written_index[removed_rows]
        )

        if removed_index.size:
            frames.append(
                removal_frame(
                    sample_id=[sample_id] * removed_index.size,
                    original_event_index=removed_index,
                    step=STEP_QC,
                    reason=reason,
                )
            )
        it_fraction = _as_float_or_none(entry.get("it_fraction"))
        mad_fraction = _as_float_or_none(entry.get("mad_fraction"))
        it_available = bool(entry.get("it_available", False))
        it_reason = entry.get("it_reason")
        applied: list[str] = []
        if it_fraction is not None and it_fraction >= 0:
            applied.append("IT")
        if mad_fraction is not None and mad_fraction >= 0:
            applied.append("MAD")
        notes: list[str] = []
        if not it_available:
            notes.append(
                "the isolation-tree analysis did not run"
                + (f" ({it_reason})" if it_reason else "")
                + f", so determine_good_cells={settings.criteria!r} was decided on "
                + ("MAD alone" if applied == ["MAD"] else "no criterion")
                + ". This is a MAD-only pass recorded as such, not the 'all' that was "
                "requested"
            )
        pct = _as_float_or_none(entry.get("percentage_removed"))
        outcomes.append(
            SampleQCOutcome(
                sample_id=sample_id,
                n_input=n_input,
                n_removed=int(removed_index.size),
                removed_fraction=pct,
                it_fraction=it_fraction,
                mad_fraction=mad_fraction,
                criteria_applied=tuple(applied),
                it_available=it_available,
                it_reason=str(it_reason) if it_reason else None,
                notes=tuple(notes),
            )
        )

    removed = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=["sample_id", "original_event_index", "step", "reason"])
    )
    if len(removed) != int(sum(o.n_removed for o in outcomes)):
        raise PeacoQCError(
            f"the removal log holds {len(removed)} rows but the per-sample outcomes "
            f"report {int(sum(o.n_removed for o in outcomes))} removals; the two are "
            "derived from one call and must agree"
        )
    return PeacoQCResult(
        removed=removed,
        samples=tuple(outcomes),
        software=tuple(peacoqc_availability(rscript)),
        warnings=tuple(str(w) for w in _as_list(payload.get("warnings"))),
        status="ok",
        report_paths={sid: str(entry.get("report_path") or "") for sid, entry in per_file.items()},
        output_files={
            sid: [str(f) for f in _as_list(entry.get("output_files"))]
            for sid, entry in per_file.items()
        },
    )


def _read_identity_channel(path: str, sample_id: str) -> np.ndarray:
    """The event indices the file was written with, read back from the file.

    Read from the file rather than assumed: if the identity channel does not come back,
    or comes back non-finite, PeacoQC's row positions cannot be resolved to events and
    the run is refused instead of reporting a set of positions that happen to look like
    event indices.
    """

    contents = read_fcs(path)
    if IDENTITY_CHANNEL not in contents.channel_names:
        raise PeacoQCError(
            f"sample {sample_id!r}: the file PeacoQC read back has no "
            f"{IDENTITY_CHANNEL!r} channel, so its removal positions cannot be resolved "
            "to events"
        )
    col = contents.channel_names.index(IDENTITY_CHANNEL)
    written = contents.events[:, col]
    if not np.isfinite(written).all():
        raise PeacoQCError(
            f"sample {sample_id!r}: the {IDENTITY_CHANNEL!r} channel came back with a "
            "non-finite value; the removal positions are not resolvable"
        )
    if written.shape[0] == 0:
        raise PeacoQCError(
            f"sample {sample_id!r}: the file PeacoQC read back holds no events; a QC "
            "call over an empty file reports success and means nothing"
        )
    return np.rint(written).astype(np.int64)


def provenance_of(result: PeacoQCResult, work_dir: str | Path) -> list:
    """Provenance refs for a PeacoQC run: the log and every input file's checksum."""

    from batchdx.cyto.schemas.contracts import ProvenanceRef

    refs = [ProvenanceRef(kind="r_log", identifier=R_SCRIPT, path=str(Path(work_dir) / "peacoqc.log"))]
    in_dir = Path(work_dir) / "fcs_in"
    if in_dir.is_dir():
        for path in sorted(in_dir.glob("*.fcs")):
            refs.append(
                ProvenanceRef(
                    kind="peacoqc_input",
                    identifier=path.stem,
                    content_hash=file_sha256(path),
                    path=str(path),
                )
            )
    del result
    return refs
