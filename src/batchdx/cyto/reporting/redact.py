"""Redaction of personal identifiers and absolute private paths.

Redaction is on by default.  The mapping from pseudonym back to the original value is
written to a separate file so that it can be kept, or not kept, independently of the
report that gets circulated.
"""

from __future__ import annotations

import hashlib
import json
import secrets
from pathlib import Path


def run_salt(run_dir: str | Path) -> str:
    """A per-run random salt, persisted beside the identifier map.

    Pseudonyms derived from a shipped, public salt are deterministic across every
    deployment and trivially reversible for low-entropy identifiers, so each run gets
    its own secret.  The salt is stored only in the run directory, next to the
    pseudonym mapping that is kept separately from any circulated report.
    """

    p = Path(run_dir) / "identifier_salt.json"
    if p.exists():
        try:
            salt = json.loads(p.read_text(encoding="utf-8")).get("salt", "")
            if salt:
                return str(salt)
        except (OSError, ValueError):  # pragma: no cover - defensive
            pass
    salt = secrets.token_hex(16)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"note": "secret salt for this run's pseudonyms; treat it "
                                 "with the same care as the identifier map",
                             "salt": salt}, indent=2) + "\n", encoding="utf-8")
    return salt


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
            digest = hashlib.sha256((self.salt + "|" + text).encode()).hexdigest()[:12]
            self._map[text] = f"{kind}-{digest}"
        return self._map[text]

    def strip_paths(self, obj):
        """Recursively replace path-like strings in a JSON-ready structure.

        Used on the resolved configuration before it is embedded in the report: input
        and checkout locations are private filesystem details and must not leak into a
        circulated document, whatever the 'identifiers redacted' pill promises.
        """

        if not self.enabled:
            return obj
        if isinstance(obj, dict):
            return {k: self.strip_paths(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [self.strip_paths(v) for v in obj]
        if isinstance(obj, str):
            if obj.startswith(("/", "~/")) or (len(obj) > 2 and obj[1] == ":"):
                return self.path(obj)
        return obj

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
