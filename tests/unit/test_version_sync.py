"""The version is declared twice (pyproject and _version); they must agree."""

from __future__ import annotations

from pathlib import Path

import tomllib

from batchdx.cyto._version import __version__


def test_pyproject_and_package_version_agree(repo_root: Path) -> None:
    pyproject = tomllib.loads((repo_root / "pyproject.toml").read_text())
    assert pyproject["project"]["version"] == __version__
