"""Locating packaged assets from an installed wheel or from a source checkout.

An installed package must find its Snakefile, R scripts, report templates and example
configurations without the user standing in the source tree.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_ASSET_DIR = _PACKAGE_DIR / "workflow_assets"


@cache
def _repo_root() -> Path | None:
    """The source checkout root, when the package is being used from one."""

    for parent in _PACKAGE_DIR.parents:
        if (parent / "workflow" / "Snakefile").exists() and (parent / "pyproject.toml").exists():
            return parent
    return None


def workflow_root() -> Path:
    """Directory holding ``Snakefile``, ``rules/``, ``scripts/``, ``envs/``, ``profiles/``."""

    installed = _ASSET_DIR / "workflow"
    if (installed / "Snakefile").exists():
        return installed
    root = _repo_root()
    if root is not None:
        return root / "workflow"
    raise FileNotFoundError(
        "the packaged workflow could not be located; expected it at "
        f"{installed} (installed) or alongside pyproject.toml (checkout)"
    )


def snakefile() -> Path:
    return workflow_root() / "Snakefile"


def config_root() -> Path:
    """Directory holding the shipped example configurations."""

    installed = _ASSET_DIR / "config"
    if installed.exists():
        return installed
    root = _repo_root()
    if root is not None:
        return root / "config"
    raise FileNotFoundError("the packaged config directory could not be located")


def examples_root() -> Path:
    installed = _ASSET_DIR / "examples"
    if installed.exists():
        return installed
    root = _repo_root()
    if root is not None:
        return root / "examples"
    raise FileNotFoundError("the packaged examples directory could not be located")


def templates_root() -> Path:
    """Report templates, which always ship inside the package itself."""

    return _PACKAGE_DIR / "reporting" / "templates"
