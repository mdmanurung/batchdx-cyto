"""Content hashes used for provenance and for candidate identity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_SEP = "\x00"


def file_sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def matrix_hash(x: np.ndarray) -> str:
    """Hash of a numeric matrix, stable across byte order and C/F layout."""

    arr = np.ascontiguousarray(np.asarray(x, dtype=np.float64))
    h = hashlib.sha256()
    h.update(str(arr.shape).encode())
    h.update(arr.tobytes())
    return h.hexdigest()


def frame_hash(df: pd.DataFrame) -> str:
    h = hashlib.sha256()
    h.update(json.dumps(list(map(str, df.columns)), sort_keys=True).encode())
    for col in df.columns:
        s = df[col]
        if pd.api.types.is_numeric_dtype(s):
            h.update(np.ascontiguousarray(s.to_numpy(dtype=np.float64)).tobytes())
        else:
            h.update(_SEP.join(map(str, s.tolist())).encode())
    return h.hexdigest()


def canonical_hash(obj: Any) -> str:
    """Hash of a JSON-serialisable object under a canonical serialisation."""

    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()
