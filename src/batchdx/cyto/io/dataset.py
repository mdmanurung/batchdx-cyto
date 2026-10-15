"""The immutable preprocessed baseline and the per-candidate output stores.

One :class:`Dataset` is loaded once per run and is never mutated.  Correction
candidates write their own AnnData files next to it; nothing copies every candidate
into a single object, and nothing overwrites the baseline.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from batchdx.cyto.config import Config
from batchdx.cyto.io.fcs import read_fcs
from batchdx.cyto.io.hashing import file_sha256, matrix_hash
from batchdx.cyto.io.manifests import (
    ManifestError,
    read_marker_table,
    read_processing_manifest,
    read_removed_events,
    read_sample_table,
)
from batchdx.cyto.schemas.contracts import (
    DatasetManifest,
    MarkerRecord,
    ProcessingManifest,
    ProvenanceRef,
    SampleRecord,
)

CELL_ID = "cell_id"
SAMPLE_ID = "sample_id"
EVENT_INDEX = "original_event_index"
SOURCE_CHECKSUM = "source_checksum"
ACQUISITION_ORDER = "acquisition_order"
#: FCS Time channel, ingested when present; the clock-time acquisition axis.
ACQUISITION_TIME = "acquisition_time"


class DatasetValidationError(ValueError):
    """Raised when the input cannot be trusted as a preprocessed baseline."""


@dataclass
class Dataset:
    """A validated, immutable baseline plus everything derived from the manifests."""

    adata: ad.AnnData
    manifest: DatasetManifest
    expression_layer: str
    samples: dict[str, SampleRecord]
    markers: dict[str, MarkerRecord]
    processing: ProcessingManifest
    removed_events: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=["sample_id", "original_event_index", "step", "reason"]
    ))
    _baseline_hash: str = ""

    # -- accessors -------------------------------------------------------------------

    @property
    def n_cells(self) -> int:
        return int(self.adata.n_obs)

    @property
    def obs(self) -> pd.DataFrame:
        """A copy of the cell-level table.  Callers cannot mutate the baseline."""
        return self.adata.obs.copy()

    def feature_index(self, features: list[str]) -> np.ndarray:
        var_names = list(self.adata.var_names)
        missing = [f for f in features if f not in var_names]
        if missing:
            raise DatasetValidationError(
                f"feature(s) {missing} are not in the dataset; a missing marker is not "
                "zero expression and is never filled in"
            )
        pos = {name: i for i, name in enumerate(var_names)}
        return np.array([pos[f] for f in features], dtype=int)

    def expression(
        self, features: list[str] | None = None, cell_ids: list[str] | None = None
    ) -> np.ndarray:
        """Return a *copy* of the baseline expression for the requested slice."""

        mat = self.adata.layers[self.expression_layer]
        mat = np.asarray(mat, dtype=np.float64)
        if cell_ids is not None:
            pos = {cid: i for i, cid in enumerate(self.adata.obs[CELL_ID])}
            missing = [c for c in cell_ids if c not in pos]
            if missing:
                raise DatasetValidationError(f"unknown cell id(s): {missing[:5]}")
            mat = mat[np.array([pos[c] for c in cell_ids], dtype=int), :]
        if features is not None:
            mat = mat[:, self.feature_index(features)]
        return np.array(mat, copy=True)

    def baseline_hash(self) -> str:
        return self._baseline_hash

    def assert_baseline_unchanged(self) -> None:
        """Fail loudly if any read path has mutated the baseline matrix."""

        current = matrix_hash(np.asarray(self.adata.layers[self.expression_layer]))
        if current != self._baseline_hash:
            raise DatasetValidationError(
                "the preprocessed baseline was modified in memory; it is immutable by contract"
            )

    def marker_ids(self, mask: str | None = None) -> list[str]:
        return self.manifest.marker_ids(mask)

    def sample_batches(self) -> dict[str, str | None]:
        return {sid: rec.batch for sid, rec in self.samples.items()}


# --------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------


def _validate_panels(markers: list[MarkerRecord]) -> None:
    panels = {m.panel_id for m in markers if m.panel_id}
    if len(panels) > 1:
        raise DatasetValidationError(
            f"the marker manifest declares more than one panel_id ({sorted(panels)}). "
            "Compatible panels must be run separately unless an explicitly validated "
            "shared-marker mode is requested; that mode is not implemented in this build."
        )


def _validate_processing(processing: ProcessingManifest, config: Config) -> None:
    if processing.scale in {"", "unknown", None}:
        raise DatasetValidationError(
            "processing manifest: 'scale' must state the coordinate system of the values "
            "(for example 'asinh_cofactor5'); it is never inferred"
        )
    if processing.corrected:
        raise DatasetValidationError(
            "processing manifest declares corrected=true. The input to this workflow is the "
            "uncorrected preprocessed baseline; chained correction is not supported by default."
        )
    if (
        config.preprocessing.enabled
        and config.preprocessing.transformation not in (None, "none")
        and processing.transformation not in (None, "none")
    ):
        raise DatasetValidationError(
            f"the processing manifest already declares transformation "
            f"{processing.transformation!r} while preprocessing.transformation is "
            f"{config.preprocessing.transformation!r}; refusing to transform twice"
        )
    scale_match = re.fullmatch(
        r"asinh_cofactor([0-9.]+(?:[eE][+-]?[0-9]+)?)", str(processing.scale)
    )
    if scale_match and processing.transformation == "asinh":
        declared = processing.transformation_parameters.get("cofactor")
        if declared is not None and not np.isclose(
            float(scale_match.group(1)), float(declared), rtol=0, atol=1e-9
        ):
            raise DatasetValidationError(
                f"processing manifest: scale {processing.scale!r} declares cofactor "
                f"{float(scale_match.group(1))} but transformation_parameters.cofactor "
                f"is {float(declared)}; the declared scale and the declared parameter "
                "must agree"
            )


def _build_cell_index(obs: pd.DataFrame, samples: dict[str, SampleRecord]) -> pd.DataFrame:
    if SAMPLE_ID not in obs.columns:
        raise DatasetValidationError(f"the expression object has no '{SAMPLE_ID}' column in obs")
    obs = obs.copy()
    obs[SAMPLE_ID] = obs[SAMPLE_ID].astype(str)

    unknown = sorted(set(obs[SAMPLE_ID]) - set(samples))
    if unknown:
        raise DatasetValidationError(
            f"sample id(s) {unknown} appear in the data but not in the sample manifest"
        )
    absent = sorted(set(samples) - set(obs[SAMPLE_ID]))
    if absent:
        raise DatasetValidationError(
            f"sample id(s) {absent} are declared in the sample manifest but carry no events"
        )

    if EVENT_INDEX not in obs.columns:
        obs[EVENT_INDEX] = obs.groupby(SAMPLE_ID).cumcount()
    obs[EVENT_INDEX] = obs[EVENT_INDEX].astype(int)

    if SOURCE_CHECKSUM not in obs.columns:
        obs[SOURCE_CHECKSUM] = obs[SAMPLE_ID].map(
            {sid: (rec.source_checksum or "") for sid, rec in samples.items()}
        )

    if CELL_ID not in obs.columns:
        obs[CELL_ID] = obs[SAMPLE_ID].astype(str) + ":" + obs[EVENT_INDEX].astype(str)
    obs[CELL_ID] = obs[CELL_ID].astype(str)

    dup = obs[CELL_ID][obs[CELL_ID].duplicated()]
    if len(dup):
        raise DatasetValidationError(f"duplicate cell id(s): {sorted(dup.unique())[:5]}")
    return obs


def load_dataset(
    config: Config,
    base_dir: str | Path = ".",
    *,
    removed_events: str | Path | None = None,
) -> Dataset:
    """Load, validate and freeze the preprocessed baseline described by ``config``.

    ``removed_events`` is an **explicit** path to the per-event removal log the
    preprocessing step wrote.  It is a parameter and not a key in the processing manifest
    because the manifest declares *what has been done to the values* and the log records
    *which events that removed*; the two are different statements, and
    ``read_processing_manifest`` refuses a key its model does not declare.  Making the
    caller pass the path also keeps the distinction that matters here visible: when no
    path is given, ``Dataset.removed_events`` is empty because **nothing was supplied**,
    which is not the same claim as *nothing was removed*, and the caller is the only
    party that knows which of the two happened.
    """

    base = Path(base_dir)
    samples = {r.sample_id: r for r in read_sample_table(base / config.input.samples)}
    markers = read_marker_table(base / config.input.markers)
    processing = read_processing_manifest(base / config.input.processing_manifest)

    _validate_panels(markers)
    _validate_processing(processing, config)

    if config.input.format == "h5ad":
        adata, provenance = _load_h5ad(base / config.input.path, config.input.expression_layer)
    else:
        adata, provenance = _load_fcs_directory(
            base / config.input.path, config.input.fcs_glob, samples, markers,
            config.input.expression_layer,
        )

    # After the cell index exists, because the log is checked against it.
    removal_log = (
        _empty_removal_log() if removed_events is None
        else _load_removal_log(base, removed_events, adata.obs)
    )

    feature_ids = [m.feature_id for m in markers]
    var_names = list(map(str, adata.var_names))
    if var_names != feature_ids:
        missing = [f for f in feature_ids if f not in var_names]
        extra = [f for f in var_names if f not in feature_ids]
        if missing or extra:
            raise DatasetValidationError(
                "the marker manifest and the expression object disagree. "
                f"Declared but absent: {missing}. Present but undeclared: {extra}. "
                "A missing marker is not zero expression."
            )
        adata = adata[:, feature_ids].copy()

    obs = _build_cell_index(adata.obs, samples)
    adata.obs = obs
    adata.obs_names = pd.Index(obs[CELL_ID].to_numpy(), dtype=object)

    layer = config.input.expression_layer
    if layer not in adata.layers:
        raise DatasetValidationError(
            f"expression layer {layer!r} is not present; available layers: "
            f"{sorted(name for name in adata.layers if name is not None)}. "
            "The layer is declared, never guessed."
        )
    mat = np.asarray(adata.layers[layer], dtype=np.float64)
    if mat.shape != (adata.n_obs, adata.n_vars):
        raise DatasetValidationError(
            f"layer {layer!r} has shape {mat.shape}, expected {(adata.n_obs, adata.n_vars)}"
        )
    if not np.isfinite(mat).all():
        bad = np.argwhere(~np.isfinite(mat))
        raise DatasetValidationError(
            f"layer {layer!r} contains {len(bad)} non-finite value(s); the first is at "
            f"cell {adata.obs[CELL_ID].iloc[bad[0][0]]}, feature {var_names[bad[0][1]]}"
        )
    adata.layers[layer] = mat

    for sid, rec in samples.items():
        rec.n_cells = int((adata.obs[SAMPLE_ID] == sid).sum())

    manifest = DatasetManifest(
        dataset_id=config.project.id,
        n_cells=int(adata.n_obs),
        n_features=int(adata.n_vars),
        expression_layer=layer,
        samples=list(samples.values()),
        markers=markers,
        processing=processing,
        baseline_hash=matrix_hash(mat),
        provenance=provenance,
    )

    return Dataset(
        adata=adata,
        manifest=manifest,
        expression_layer=layer,
        samples=samples,
        markers={m.feature_id: m for m in markers},
        processing=processing,
        removed_events=removal_log,
        _baseline_hash=manifest.baseline_hash or "",
    )


def _empty_removal_log() -> pd.DataFrame:
    """The four declared columns and no rows.

    The model import is local for the cycle reason given in :func:`_load_removal_log`;
    this one is called on every load, so the statement is kept in one place.
    """

    from batchdx.cyto.preprocess.removals import empty_removal_log

    return empty_removal_log()


def _load_removal_log(base: Path, removed_events: str | Path, obs: pd.DataFrame) -> pd.DataFrame:
    """Read the removal log and check that it describes *this* dataset's events.

    A log naming an event the baseline does not contain is worse than no log: it would
    let a reader conclude that a cell was removed when in fact the log belongs to
    another run.  So every ``(sample_id, original_event_index)`` pair is checked against
    the loaded cell index and a mismatch is a refusal, not a warning.
    """

    # Imported here, not at module scope: batchdx.cyto.preprocess.peacoqc needs
    # batchdx.cyto.adapters, whose package __init__ imports adapters.base, which imports
    # this module. A top-level import of the removal model would be a cycle whenever io is
    # imported first. See the same note in io/manifests.py.
    from batchdx.cyto.preprocess.removals import empty_removal_log

    path = Path(removed_events)
    if not path.is_absolute():
        path = base / path
    try:
        log = read_removed_events(path)
    except ManifestError as exc:
        raise DatasetValidationError(str(exc)) from exc

    pairs = set(
        zip(
            obs[SAMPLE_ID].astype(str).to_numpy(),
            obs[EVENT_INDEX].astype(int).to_numpy(),
            strict=True,
        )
    )
    unknown = sorted(
        {
            (str(s), int(i))
            for s, i in zip(log["sample_id"], log["original_event_index"], strict=True)
        }
        - pairs
    )
    if unknown:
        raise DatasetValidationError(
            f"removal log {path} names {len(unknown)} event(s) the loaded baseline does not "
            f"contain, first {unknown[0]}. A log describing another run's events would let "
            "a reader conclude that a cell was removed when it was never here. The log and "
            "the baseline must be the same run"
        )
    return log if not log.empty else empty_removal_log()


def _load_h5ad(path: Path, layer: str) -> tuple[ad.AnnData, list[ProvenanceRef]]:
    if not path.exists():
        raise DatasetValidationError(f"input h5ad not found: {path}")
    adata = ad.read_h5ad(path)
    if layer not in adata.layers:
        # anndata >= 0.13 exposes a ``None`` key in ``.layers`` that aliases ``.X``,
        # so sorting the mapping directly raises TypeError on mixed None/str keys.
        # Report the real layer names and say that ``X`` is the default.
        available = sorted(name for name in adata.layers if name is not None)
        shown = ", ".join(available) if available else "(none)"
        raise DatasetValidationError(
            f"{path}: expression layer {layer!r} is absent; available layers "
            f"[{shown}]. Declare the layer that holds the values to analyse."
        )
    prov = [ProvenanceRef(kind="input_h5ad", identifier=path.name,
                          content_hash=file_sha256(path), path=str(path))]
    return adata, prov


def _load_fcs_directory(
    directory: Path,
    pattern: str,
    samples: dict[str, SampleRecord],
    markers: list[MarkerRecord],
    layer: str,
) -> tuple[ad.AnnData, list[ProvenanceRef]]:
    """Assemble an AnnData from per-sample FCS files named in the sample manifest."""

    if not directory.is_dir():
        raise DatasetValidationError(f"input FCS directory not found: {directory}")

    channel_to_feature: dict[str, str] = {}
    for m in markers:
        if not m.channel:
            raise DatasetValidationError(
                f"marker manifest: feature {m.feature_id!r} has no channel, which is required "
                "when reading FCS input"
            )
        channel_to_feature[m.channel] = m.feature_id

    blocks: list[np.ndarray] = []
    obs_rows: list[pd.DataFrame] = []
    prov: list[ProvenanceRef] = []

    for sid, rec in samples.items():
        if not rec.source_file:
            raise DatasetValidationError(
                f"sample {sid!r}: source_file is required for FCS input and is never "
                "inferred from file names"
            )
        path = directory / rec.source_file
        contents = read_fcs(path)
        checksum = file_sha256(path)
        rec.source_checksum = checksum
        prov.append(ProvenanceRef(kind="input_fcs", identifier=sid,
                                  content_hash=checksum, path=str(path)))

        missing = [c for c in channel_to_feature if c not in contents.channel_names]
        if missing:
            raise DatasetValidationError(
                f"{path}: declared channel(s) {missing} are absent from the file"
            )
        idx = [contents.channel_names.index(c) for c in channel_to_feature]
        blocks.append(contents.events[:, idx])

        n = contents.n_events
        row = {
            SAMPLE_ID: np.repeat(sid, n),
            EVENT_INDEX: np.arange(n, dtype=int),
            SOURCE_CHECKSUM: np.repeat(checksum, n),
            # The event index stays the acquisition-order proxy; the FCS Time channel,
            # when present, is ingested separately as clock time.
            ACQUISITION_ORDER: np.arange(n, dtype=int),
        }
        if contents.time_values is not None:
            row[ACQUISITION_TIME] = np.asarray(contents.time_values, dtype=float)
        obs_rows.append(pd.DataFrame(row))

    mat = np.vstack(blocks) if blocks else np.zeros((0, len(channel_to_feature)))
    obs = pd.concat(obs_rows, ignore_index=True)
    adata = ad.AnnData(
        X=np.zeros((mat.shape[0], mat.shape[1]), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(list(channel_to_feature.values()), dtype=object)),
    )
    adata.layers[layer] = mat
    return adata, prov


# --------------------------------------------------------------------------------------
# Candidate output store
# --------------------------------------------------------------------------------------


def write_candidate_output(
    dataset: Dataset,
    matrix: np.ndarray,
    path: str | Path,
    *,
    candidate_id: str,
    layer: str = "corrected",
    output_scale: str | None = None,
) -> Path:
    """Write one candidate's corrected matrix as its own AnnData file.

    The cell index, marker order and sample metadata are copied verbatim from the
    baseline so that event identity is verifiable downstream.
    """

    if matrix.shape != (dataset.n_cells, len(dataset.adata.var_names)):
        raise DatasetValidationError(
            f"candidate {candidate_id!r}: corrected matrix has shape {matrix.shape}, "
            f"expected {(dataset.n_cells, len(dataset.adata.var_names))}"
        )
    if not np.isfinite(matrix).all():
        raise DatasetValidationError(
            f"candidate {candidate_id!r}: corrected matrix contains non-finite values"
        )

    out = ad.AnnData(
        X=np.zeros(matrix.shape, dtype=np.float32),
        obs=dataset.adata.obs.copy(),
        var=dataset.adata.var.copy(),
    )
    out.layers[layer] = np.asarray(matrix, dtype=np.float64)
    out.uns["candidate_id"] = candidate_id
    out.uns["baseline_hash"] = dataset.baseline_hash()
    out.uns["output_scale"] = output_scale or dataset.processing.scale
    out.uns["expression_layer"] = layer

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(p)
    return p
