"""Readers for the sample, marker and processing manifests.

Only the metadata an operation actually needs is required.  Absent donor
information, for example, still permits descriptive diagnosis; it disables
donor-matched metrics rather than causing donors to be manufactured.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from batchdx.cyto.schemas.contracts import MarkerRecord, ProcessingManifest, SampleRecord
from batchdx.cyto.schemas.enums import MarkerRole, ReferenceRole

SAMPLE_REQUIRED = ["sample_id"]
MARKER_REQUIRED = ["feature_id"]

_TRUTHY = {"1", "true", "yes", "y", "t"}
_FALSY = {"0", "false", "no", "n", "f", ""}


class ManifestError(ValueError):
    """Raised when a manifest is missing, malformed or internally inconsistent."""


def _read_table(path: str | Path, what: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise ManifestError(f"{what} manifest not found: {p}")
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False, na_values=[""])
    df.columns = [c.strip() for c in df.columns]
    if df.empty:
        raise ManifestError(f"{what} manifest {p} has no rows")
    return df


def _as_bool(value, column: str, row: int) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    text = str(value).strip().lower()
    if text in _TRUTHY:
        return True
    if text in _FALSY:
        return False
    raise ManifestError(
        f"marker manifest row {row}: column {column!r} has value {value!r}; "
        "expected a boolean (true/false)"
    )


def _opt(row: pd.Series, key: str) -> str | None:
    if key not in row:
        return None
    value = row[key]
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


def read_sample_table(path: str | Path) -> list[SampleRecord]:
    """Read the sample manifest, keeping donor/biospecimen/aliquot identity distinct."""

    df = _read_table(path, "sample")
    missing = [c for c in SAMPLE_REQUIRED if c not in df.columns]
    if missing:
        raise ManifestError(f"sample manifest {path}: missing required column(s) {missing}")

    ids = df["sample_id"].astype(str).str.strip()
    dupes = sorted(ids[ids.duplicated()].unique())
    if dupes:
        raise ManifestError(f"sample manifest {path}: duplicate sample_id(s) {dupes}")

    known = {
        "sample_id", "donor_id", "biospecimen_id", "aliquot_id", "batch", "group",
        "stimulation", "timepoint", "reference_id", "reference_role", "n_cells",
        "source_file", "source_checksum",
    }
    technical_cols = [c for c in df.columns if c not in known]

    records: list[SampleRecord] = []
    for i, row in df.iterrows():
        role_text = _opt(row, "reference_role") or "none"
        try:
            role = ReferenceRole(role_text)
        except ValueError as exc:
            raise ManifestError(
                f"sample manifest {path} row {i}: reference_role {role_text!r} is not one of "
                f"{[r.value for r in ReferenceRole]}"
            ) from exc
        if role is not ReferenceRole.NONE and _opt(row, "reference_id") is None:
            raise ManifestError(
                f"sample manifest {path} row {i}: reference_role {role_text!r} declared without a "
                "reference_id; reference identity is never inferred from file names"
            )
        n_cells = _opt(row, "n_cells")
        records.append(
            SampleRecord(
                sample_id=str(row["sample_id"]).strip(),
                donor_id=_opt(row, "donor_id"),
                biospecimen_id=_opt(row, "biospecimen_id"),
                aliquot_id=_opt(row, "aliquot_id"),
                batch=_opt(row, "batch"),
                technical_factors={c: str(row[c]) for c in technical_cols if _opt(row, c)},
                group=_opt(row, "group"),
                stimulation=_opt(row, "stimulation"),
                timepoint=_opt(row, "timepoint"),
                reference_id=_opt(row, "reference_id"),
                reference_role=role,
                n_cells=int(n_cells) if n_cells else None,
                source_file=_opt(row, "source_file"),
                source_checksum=_opt(row, "source_checksum"),
            )
        )
    return records


def read_marker_table(path: str | Path) -> list[MarkerRecord]:
    """Read the marker manifest with independent clustering/correction/evaluation masks."""

    df = _read_table(path, "marker")
    missing = [c for c in MARKER_REQUIRED if c not in df.columns]
    if missing:
        raise ManifestError(f"marker manifest {path}: missing required column(s) {missing}")

    ids = df["feature_id"].astype(str).str.strip()
    dupes = sorted(ids[ids.duplicated()].unique())
    if dupes:
        raise ManifestError(f"marker manifest {path}: duplicate feature_id(s) {dupes}")

    if "channel" in df.columns:
        channels = df["channel"].astype(str).str.strip()
        chan_dupes = sorted(channels[channels.duplicated() & (channels != "")].unique())
        if chan_dupes:
            raise ManifestError(
                f"marker manifest {path}: channel(s) {chan_dupes} map to more than one feature_id; "
                "the channel to feature mapping must be unambiguous"
            )

    records: list[MarkerRecord] = []
    for i, row in df.iterrows():
        role_text = _opt(row, "role") or "other"
        try:
            role = MarkerRole(role_text)
        except ValueError as exc:
            raise ManifestError(
                f"marker manifest {path} row {i}: role {role_text!r} is not one of "
                f"{[r.value for r in MarkerRole]}"
            ) from exc
        records.append(
            MarkerRecord(
                feature_id=str(row["feature_id"]).strip(),
                channel=_opt(row, "channel"),
                marker=_opt(row, "marker"),
                panel_id=_opt(row, "panel_id"),
                role=role,
                use_for_clustering=_as_bool(row.get("use_for_clustering"), "use_for_clustering", i),
                use_for_correction=_as_bool(row.get("use_for_correction"), "use_for_correction", i),
                use_for_evaluation=_as_bool(row.get("use_for_evaluation"), "use_for_evaluation", i),
            )
        )
    return records


def read_processing_manifest(path: str | Path) -> ProcessingManifest:
    """Read the declared processing state.  Nothing here is inferred from the values."""

    p = Path(path)
    if not p.exists():
        raise ManifestError(f"processing manifest not found: {p}")
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ManifestError(f"processing manifest {p}: not valid JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise ManifestError(f"processing manifest {p}: top level must be an object")
    unknown = set(raw) - set(ProcessingManifest.model_fields)
    if unknown:
        raise ManifestError(
            f"processing manifest {p}: unknown key(s) {sorted(unknown)}; "
            f"expected a subset of {sorted(ProcessingManifest.model_fields)}"
        )
    if "scale" not in raw:
        raise ManifestError(
            f"processing manifest {p}: 'scale' must be declared explicitly; the package "
            "does not guess whether values are raw, transformed or corrected"
        )
    return ProcessingManifest.model_validate(raw)


# --------------------------------------------------------------------------------------
# The removal log
# --------------------------------------------------------------------------------------
#
# ``Dataset.removed_events`` has declared these four columns since the field was added
# and had no writer and no reader anywhere in the repository.  The reader and the writer
# live here, next to the other manifest readers, rather than in
# ``batchdx.cyto.preprocess.removals`` where the model is, because a file next to the
# processing manifest is this module's business.
#
# The import of the removal model is inside the functions.  ``preprocess.peacoqc`` needs
# ``batchdx.cyto.adapters``, whose package __init__ imports ``adapters.base``, which
# imports ``io.dataset``, which imports this module -- so a module-level import here is a
# cycle whenever ``io`` is imported first.  Importing it lazily costs nothing at import
# time and keeps both orders working.


def read_removed_events(path: str | Path) -> pd.DataFrame:
    """Read a per-event removal log written by :func:`write_removed_events`.

    Raises :class:`ManifestError` when the file is absent, unreadable or does not carry
    exactly the four declared columns.  It never returns an empty frame for a file it
    could not read: an absent removal log and a removal log recording nothing are
    different facts, and collapsing them is how a dropped cell becomes invisible.
    """

    from batchdx.cyto.preprocess.removals import check_removal_log

    p = Path(path)
    if not p.exists():
        raise ManifestError(
            f"removal log not found: {p}. A removal log is written by "
            "batchdx.cyto.preprocess (gating and PeacoQC); if the baseline was gated and "
            "QCed upstream, the log must come with it. No log is not a log of no removals"
        )
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    try:
        df = pd.read_csv(p, sep=sep, dtype={"sample_id": str, "step": str, "reason": str})
    except Exception as exc:  # pandas raises a family of them
        raise ManifestError(f"removal log {p} could not be read: {exc}") from exc
    if df.empty:
        raise ManifestError(
            f"removal log {p} has a header and no rows. A zero-row log is a claim that "
            "nothing was removed; if that is what happened, say so in the processing "
            "manifest's qc_history instead of shipping an empty file that reads as a "
            "record"
        )
    try:
        return check_removal_log(df, what=f"removal log {p}")
    except ValueError as exc:
        raise ManifestError(str(exc)) from exc


def write_removed_events(frame: pd.DataFrame, path: str | Path) -> Path:
    """Write a removal log in the form :func:`read_removed_events` accepts.

    The columns, their order, the integer dtype of ``original_event_index`` and the
    string dtype of ``sample_id`` are all checked before anything is written, so a log
    that cannot be read back is never written in the first place.  The extension chooses
    the separator, and a ``.csv`` suffix is required for CSV so that a reader cannot
    guess wrongly.
    """

    from batchdx.cyto.preprocess.removals import check_removal_log

    p = Path(path)
    if p.suffix.lower() not in {".csv", ".tsv", ".txt"}:
        raise ManifestError(
            f"removal log {p}: suffix must be .csv, .tsv or .txt, so the separator is "
            "declared rather than guessed from a name that could be anything"
        )
    checked = check_removal_log(frame, what=f"removal log for {p}")
    p.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    checked.to_csv(p, sep=sep, index=False)
    return p
