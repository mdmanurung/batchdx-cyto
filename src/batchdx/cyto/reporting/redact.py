"""Redaction of personal identifiers and absolute private paths.

Redaction is on by default.  The mapping from pseudonym back to the original value is
written to a separate file so that it can be kept, or not kept, independently of the
report that gets circulated.

Two kinds of value are redacted, and they are not interchangeable.  An *identifier*
becomes a stable pseudonym, so two rows that share a sample still visibly share one.
A *path* becomes a location relative to the run's own base, which is not reversible at
all: the point is that the report cannot carry the directory layout of the machine that
produced it, not that a reader can look the file up afterwards.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


class Redactor:
    """Stable pseudonyms for identifiers, plus private-path stripping."""

    def __init__(self, enabled: bool = True, base_dir: str | Path = ".", salt: str = "batchdx"):
        self.enabled = enabled
        self.base = Path(base_dir).resolve()
        self.salt = salt
        self._map: dict[str, str] = {}

    def identifier(self, value: str | None, kind: str = "ID") -> str:
        if value is None:
            return ""
        if not self.enabled:
            return str(value)
        text = str(value)
        if text not in self._map:
            digest = hashlib.sha256((self.salt + "|" + text).encode()).hexdigest()[:8]
            self._map[text] = f"{kind}-{digest}"
        return self._map[text]

    def path(self, value: str | None) -> str:
        """Rewrite one path as ``<run>/…`` when it lies inside the run's own base.

        A relative path is resolved against ``base``, not against the process's working
        directory: a configuration is written relative to the project it describes, so
        ``data/preprocessed.h5ad`` means ``<project>/data/preprocessed.h5ad`` and would
        otherwise be resolved against whatever directory the report happened to be
        rendered from, which makes the same report redact differently on two machines.

        A path outside the base keeps its file name and loses everything above it, so
        the reader still learns which file was read and not which machine it came from.
        """

        if value is None:
            return ""
        text = str(value)
        if not text.strip():
            return text
        if not self.enabled:
            return text
        candidate = Path(text).expanduser()
        if not candidate.is_absolute():
            candidate = self.base / candidate
        try:
            resolved = candidate.resolve()
        except OSError:  # pragma: no cover - defensive
            return f"<redacted>/{candidate.name}"
        try:
            return str(Path("<run>") / resolved.relative_to(self.base))
        except ValueError:
            return f"<redacted>/{resolved.name}"

    def write_map(self, path: str | Path) -> Path | None:
        """Write the pseudonym mapping, or nothing when redaction is off."""

        if not self.enabled or not self._map:
            return None
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(
            json.dumps({"note": "keep this file separate from any circulated report",
                        "mapping": self._map}, indent=2) + "\n",
            encoding="utf-8",
        )
        return p
