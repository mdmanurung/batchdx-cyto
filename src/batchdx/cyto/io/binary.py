"""Typed binary matrix exchange between Python and R.

Large delimited-text intermediates are avoided: a matrix is written as a raw
little-endian float64 block with a JSON sidecar describing shape, row ids and column
names.  The R side reads it with ``readBin`` alone, so the exchange adds no R package
dependency of its own.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

MAGIC = "batchdx-binmat-1"


def write_matrix(
    path: str | Path,
    matrix: np.ndarray,
    row_ids: list[str],
    col_names: list[str],
    extra: dict | None = None,
) -> tuple[Path, Path]:
    """Write ``matrix`` as float64 plus a ``.json`` sidecar.  Returns both paths."""

    arr = np.ascontiguousarray(np.asarray(matrix, dtype="<f8"))
    if arr.ndim != 2:
        raise ValueError(f"expected a 2-d matrix, got shape {arr.shape}")
    if arr.shape[0] != len(row_ids) or arr.shape[1] != len(col_names):
        raise ValueError(
            f"shape {arr.shape} does not match {len(row_ids)} row ids and "
            f"{len(col_names)} column names"
        )
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(arr.tobytes(order="C"))
    sidecar = p.with_suffix(p.suffix + ".json")
    sidecar.write_text(
        json.dumps({
            "magic": MAGIC,
            "dtype": "float64",
            "byte_order": "little",
            "layout": "row_major",
            "n_rows": int(arr.shape[0]),
            "n_cols": int(arr.shape[1]),
            "row_ids": list(row_ids),
            "col_names": list(col_names),
            **(extra or {}),
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    return p, sidecar


def read_matrix(path: str | Path) -> tuple[np.ndarray, list[str], list[str], dict]:
    """Read a matrix written by :func:`write_matrix`."""

    p = Path(path)
    sidecar = p.with_suffix(p.suffix + ".json")
    if not sidecar.exists():
        raise FileNotFoundError(f"sidecar not found for {p}: {sidecar}")
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    if meta.get("magic") != MAGIC:
        raise ValueError(f"{sidecar}: unexpected magic {meta.get('magic')!r}")
    arr = np.frombuffer(p.read_bytes(), dtype="<f8")
    arr = arr.reshape(int(meta["n_rows"]), int(meta["n_cols"]))
    return np.array(arr, copy=True), list(meta["row_ids"]), list(meta["col_names"]), meta
