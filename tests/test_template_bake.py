"""Tests that the Cookiecutter template bakes successfully (without running hooks)."""

from __future__ import annotations

from pathlib import Path

import pytest
import toml
from cookiecutter.main import cookiecutter

REPO_ROOT = Path(__file__).resolve().parent.parent

EXTRA_CONTEXT: dict[str, str] = {
    "project_name": "Cookiecutter Test Project",
    "short_description": "A project for bake tests.",
    "github_username": "testuser",
    "full_name": "Test User",
    "email": "test@example.com",
}

EXPECTED_REPO_NAME = "cookiecutter-test-project"
EXPECTED_MODULE_NAME = "cookiecutter_test_project"


def _assert_no_jinja_leaks(project_root: Path) -> None:
    """Fail if any generated Python or TOML still contains unrendered ``{{``."""
    for path in project_root.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        assert "{{" not in text, f"Unrendered Jinja in {path}"

    for path in project_root.rglob("*.toml"):
        text = path.read_text(encoding="utf-8")
        assert "{{" not in text, f"Unrendered Jinja in {path}"


@pytest.fixture
def baked_project(tmp_path: Path) -> Path:
    """Generate the template into a temporary directory (hooks disabled)."""
    result = cookiecutter(
        str(REPO_ROOT),
        no_input=True,
        extra_context=EXTRA_CONTEXT,
        output_dir=str(tmp_path),
        overwrite_if_exists=True,
        accept_hooks=False,
    )
    project_dir = tmp_path / EXPECTED_REPO_NAME
    assert Path(result).resolve() == project_dir.resolve()
    return project_dir


def test_cookiecutter_bake_creates_expected_layout(baked_project: Path) -> None:
    """Baked project contains core files and package layout."""
    assert baked_project.is_dir()

    assert (baked_project / "pyproject.toml").is_file()
    assert (baked_project / "README.md").is_file()
    assert (baked_project / ".pre-commit-config.yaml").is_file()

    pkg = baked_project / EXPECTED_MODULE_NAME
    assert pkg.is_dir()
    assert (pkg / "__init__.py").is_file()
    assert (pkg / "cli.py").is_file()

    tests_dir = baked_project / "tests"
    assert tests_dir.is_dir()
    assert (tests_dir / "test_cli.py").is_file()
    assert (tests_dir / "test_models.py").is_file()
    assert (tests_dir / "test___init__.py").is_file()


def test_baked_pyproject_metadata(baked_project: Path) -> None:
    """pyproject.toml parses and matches Hatch + project identity."""
    data = toml.load(baked_project / "pyproject.toml")
    assert data["build-system"]["build-backend"] == "hatchling.build"
    assert "hatchling" in data["build-system"]["requires"][0]
    assert data["project"]["name"] == EXPECTED_REPO_NAME
    assert data["tool"]["hatch"]["build"]["targets"]["wheel"]["packages"] == [
        EXPECTED_MODULE_NAME
    ]


def test_baked_project_has_no_jinja_leaks(baked_project: Path) -> None:
    """Rendered sources do not contain raw Jinja markers."""
    _assert_no_jinja_leaks(baked_project)
