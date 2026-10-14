"""The version is declared twice (pyproject and _version); they must agree."""

from __future__ import annotations

import tomllib
from pathlib import Path

from batchdx.cyto._version import __version__


def test_pyproject_and_package_version_agree(repo_root: Path) -> None:
    pyproject = tomllib.loads((repo_root / "pyproject.toml").read_text())
    assert pyproject["project"]["version"] == __version__
