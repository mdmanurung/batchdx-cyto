"""Shared fixtures.

Everything here is small, deterministic and synthetic: the unit and scientific suites
never need a downloaded dataset.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

TESTS_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TESTS_ROOT.parent
sys.path.insert(0, str(TESTS_ROOT))

from fixtures.backends import register_fixture_backends  # noqa: E402

register_fixture_backends()


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def upstream_dir() -> Path:
    """Where the pinned upstream checkouts live, if they were fetched."""

    return REPO_ROOT / ".upstream"


@pytest.fixture
def example_project(tmp_path: Path) -> Path:
    """A materialised paired-reference example, small enough to run in seconds."""

    from batchdx.cyto.examples import materialise_example

    materialise_example("paired_reference", tmp_path, cells_per_sample=400)
    return tmp_path


@pytest.fixture
def confounded_project(tmp_path: Path) -> Path:
    from batchdx.cyto.examples import materialise_example

    materialise_example("confounded", tmp_path, cells_per_sample=400)
    return tmp_path


@pytest.fixture
def loaded(example_project: Path):
    """``(config, dataset)`` for the paired-reference example."""

    from batchdx.cyto.config import load_config
    from batchdx.cyto.io import load_dataset

    config = load_config(example_project / "config" / "config.yaml")
    return config, load_dataset(config, example_project)


@pytest.fixture
def frozen(loaded):
    """``(config, dataset, design, sampling)`` with the common partition frozen."""

    from batchdx.cyto.design import audit_design
    from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels

    config, dataset = loaded
    design = audit_design(
        list(dataset.samples.values()),
        dataset_id=dataset.manifest.dataset_id,
        correction_factor=config.design.correction_batch,
        protected=config.design.protect,
        donor_factor=config.design.donor,
        minimum_batches=config.design.minimum_batches,
        require_shared_reference=config.design.require_shared_reference,
    )
    base = Path(dataset.manifest.provenance[0].path).parents[1]
    labels, tier, confidence = resolve_labels(dataset, config, base)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    return config, dataset, design, sampling
