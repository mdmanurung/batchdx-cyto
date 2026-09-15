"""Redaction of personal identifiers and absolute private paths.

Redaction is on by default.  The mapping from pseudonym back to the original value is
written to a separate file so that it can be kept, or not kept, independently of the
report that gets circulated.
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
        if value is None:
            return ""
        text = str(value)
        if not self.enabled:
            return text
        try:
            resolved = Path(text).resolve()
        except OSError:  # pragma: no cover - defensive
            return text
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
