"""Release metadata has to agree with itself, and has to agree with the code.

These are cheap assertions about files that ship, and they exist because each of them
was once wrong in a way nothing caught:

* ``pyproject.toml`` and ``CITATION.cff`` both declared MIT while **no licence text was
  shipped at all** -- a hard blocker for distribution, and invisible to every test in
  the suite;
* three files carried the version independently and had drifted apart before;
* ``THIRD_PARTY_NOTICES.md`` asserted that no OpenProblems implementation had been
  extracted, which was false: seven components are registered ``compatibility`` /
  ``upstream_verified=True``;
* ``workflow/envs/README.md`` said one environment file lived in the directory it
  describes, and two did;
* ``config/upstream_lock.yaml`` described CytoBatchFlagR's dependency risk as six
  packages when the enforced set is 29;
* ``containers/Dockerfile.r-backends`` installed none of them, and declared an
  environment variable that nothing in the repository reads.

The theme is the one the whole project cares about: a claim written into a document has
to be one that can be checked. Each test below therefore either measures something or
refuses to pass. Where a check needs something that may be absent (a pinned upstream
checkout, the isolated R library, ``tomllib``), it skips with a reason instead of
guessing.

Two of these tests are stronger than a string comparison: the container ones are
re-derived from the upstream ``library()`` block itself, so they fail if upstream's needs
change or if a package is dropped from the image, rather than only if this file's
prose changes.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest
import yaml

from batchdx.cyto import __version__ as runtime_version
from batchdx.cyto._version import CONTRACT_SCHEMA_VERSION
from batchdx.cyto._version import __version__ as module_version

try:  # Python 3.11+; 3.10 has no tomllib and no declared toml backport.
    import tomllib
except ImportError:  # pragma: no cover - exercised only on 3.10
    tomllib = None  # type: ignore[assignment]

requires_toml = pytest.mark.skipif(
    tomllib is None, reason="tomllib is stdlib from 3.11; this project also supports 3.10"
)

#: The R packages that ship with R itself, so no image has to install them.
R_BUNDLED = frozenset({"MASS", "stats"})

#: The anchored, falsifiable facts the two corrected documents must keep stating. These
#: are asserted as text because the failure they guard against *is* a text claim; each
#: one names a file, or a measurement, rather than restating an opinion.
UPSTREAM_LIBRARY_BLOCK = ".upstream/CytoBatchFlagR/R/CytoBatchFlagR_functions.R"
UPSTREAM_LIBRARY_LINES = "6-34"
ENFORCED_LIBRARY_COUNT = 29

#: What counts as executable or configured, for a scan that must ignore prose.
_CODE_SUFFIXES = frozenset({".py", ".R", ".r", ".sh", ".bash", ".yaml", ".yml", ".toml", ".json"})


def _pyproject(repo_root: Path) -> dict:
    assert tomllib is not None
    return tomllib.loads((repo_root / "pyproject.toml").read_text())


def _library_block_packages(upstream_dir: Path) -> list[str]:
    """The package names upstream loads unconditionally, read from the pinned checkout."""
    source = upstream_dir / "CytoBatchFlagR" / "R" / "CytoBatchFlagR_functions.R"
    if not source.exists():
        pytest.skip("run scripts/fetch_upstream.sh; the pinned CytoBatchFlagR is absent")
    return re.findall(r"library\((\w+)\)", source.read_text())


def _r_library_names(repo_root: Path) -> set[str]:
    # CI builds the same package set in a temporary micromamba prefix and exports
    # its library through R_LIBS. Resolve the executing environment first.
    library = Path(os.environ.get("R_LIBS") or repo_root / "renv" / "lib" / "R" / "library")
    if not library.is_dir():
        pytest.skip(f"the active isolated R library is absent: {library}")
    return {entry.name for entry in library.iterdir() if entry.is_dir()}


def _install_lists(repo_root: Path) -> tuple[set[str], set[str]]:
    """Every package name the container asks CRAN and Bioconductor for, in that order.

    Unioned over all the install calls, not just the largest, so adding a small extra
    CRAN step cannot silently hide a missing package from the check below.
    """
    text = (repo_root / "containers" / "Dockerfile.r-backends").read_text()
    bioc = re.findall(r"BiocManager::install\(c\(([^)]*)\)", text)
    cran = re.findall(r"install\.packages\(c\(([^)]*)\)", text)
    assert bioc, "the container no longer installs anything from Bioconductor"
    assert cran, "the container no longer installs anything from CRAN"
    return (
        {name for group in cran for name in re.findall(r"'([\w.]+)'", group)},
        {name for group in bioc for name in re.findall(r"'([\w.]+)'", group)},
    )


# --------------------------------------------------------------------------------------
# 1.9  The licence
# --------------------------------------------------------------------------------------


def test_the_mit_licence_text_ships(repo_root: Path) -> None:
    """Two files declared MIT; until now neither shipped the text it was promising."""
    licence = repo_root / "LICENSE"
    assert licence.is_file(), "LICENSE is declared in two places but no licence text ships"
    text = licence.read_text()
    assert text.startswith("MIT License"), "the licence file must be the MIT text"
    # The two clauses that make it MIT rather than a paraphrase of it.
    assert "Permission is hereby granted, free of charge" in text
    assert 'THE SOFTWARE IS PROVIDED "AS IS"' in text
    assert "Copyright (c)" in text, "an MIT licence with no copyright line is incomplete"


@requires_toml
def test_the_declared_licence_is_mit_in_all_three_places(repo_root: Path) -> None:
    project = _pyproject(repo_root)["project"]
    # PEP 639: an SPDX expression string, not the deprecated {text = ...} table.
    assert project["license"] == "MIT"
    citation = yaml.safe_load((repo_root / "CITATION.cff").read_text())
    assert citation["license"] == "MIT"


@requires_toml
def test_every_declared_licence_file_exists(repo_root: Path) -> None:
    """`license-files` is a promise about the filesystem; check it."""
    for relative in _pyproject(repo_root)["project"]["license-files"]:
        assert (repo_root / relative).is_file(), f"license-files names a missing file: {relative}"


@requires_toml
def test_the_build_backend_floor_is_pinned(repo_root: Path) -> None:
    """hatchling below 1.27 silently emits the deprecated License: field instead.

    Verified by building both ways: 1.27.0 emits ``License-Expression`` and
    ``License-File``; 1.26.3 succeeds but drops ``License-File``. An unpinned
    ``requires`` would make the published metadata depend on the resolver.
    """
    requires = _pyproject(repo_root)["build-system"]["requires"]
    hatchling = next(r for r in requires if r.startswith("hatchling"))
    assert ">=1.27" in hatchling, f"the PEP 639 floor is not pinned: {hatchling!r}"


# --------------------------------------------------------------------------------------
# 1.10  Version and classifiers
# --------------------------------------------------------------------------------------


@requires_toml
def test_the_version_is_the_same_in_every_place_it_is_written(repo_root: Path) -> None:
    """Four copies of one number: the build, the module, the package, the citation."""
    declared = _pyproject(repo_root)["project"]["version"]
    citation = yaml.safe_load((repo_root / "CITATION.cff").read_text())["version"]
    assert declared == citation, f"pyproject says {declared}, CITATION.cff says {citation}"
    assert module_version == declared, f"_version.py says {module_version}"
    assert runtime_version == declared, f"the imported package says {runtime_version}"


@requires_toml
def test_the_version_is_not_a_development_only_sentinel(repo_root: Path) -> None:
    """`0.1.0.dev0` reads as 'not releasable'. Alpha gates are real; a dev sentinel is not."""
    assert not _pyproject(repo_root)["project"]["version"].endswith(".dev0")
    assert not module_version.endswith(".dev0")


@requires_toml
def test_every_python_version_requires_python_allows_has_a_classifier(repo_root: Path) -> None:
    """The floor is a promise; the classifiers are what PyPI renders from it.

    Compared as integer pairs, not floats: `float("3.10") == 3.1`, which would make this
    test silently agree with a mismatch.
    """
    project = _pyproject(repo_root)
    floor = re.search(r">=\s*(\d+)\.(\d+)", project["project"]["requires-python"])
    assert floor is not None, "requires-python no longer states a version floor"
    lowest = (int(floor.group(1)), int(floor.group(2)))
    classified = {
        (int(m.group(1)), int(m.group(2)))
        for c in project["project"]["classifiers"]
        if (m := re.fullmatch(r"Programming Language :: Python :: (\d+)\.(\d+)", c))
    }
    assert classified, "no Python version classifiers are declared"
    assert min(classified) == lowest, (
        f"requires-python allows {lowest[0]}.{lowest[1]} but the lowest classifier is "
        f"{min(classified)[0]}.{min(classified)[1]}"
    )
    # And nothing classified may be below the floor, which `min` alone would not catch.
    assert all(version >= lowest for version in classified), (
        f"a classifier sits below the requires-python floor: {sorted(classified)}"
    )


@requires_toml
def test_r_is_declared_as_a_language_this_package_actually_uses(repo_root: Path) -> None:
    """`Programming Language :: R` is a claim that R scripts do real work here."""
    classifiers = _pyproject(repo_root)["project"]["classifiers"]
    assert "Programming Language :: R" in classifiers
    scripts = sorted((repo_root / "workflow" / "scripts").rglob("*.R"))
    scripts += sorted((repo_root / "scripts").rglob("*.R"))
    assert scripts, "the R classifier is declared but no R script is tracked"


@requires_toml
def test_the_development_status_matches_the_gate_states(repo_root: Path) -> None:
    """Several Stage gates are PARTIAL and backends record not_available.

    Alpha (3) is the honest trove value for that. Beta or later would be a claim the
    gate table in docs/implementation_status.md does not support.
    """
    classifiers = _pyproject(repo_root)["project"]["classifiers"]
    status = [c for c in classifiers if c.startswith("Development Status ::")]
    assert status == ["Development Status :: 3 - Alpha"], (
        f"the development status drifted away from the gate states: {status}"
    )


def test_the_contract_and_config_schema_versions_are_still_declared(repo_root: Path) -> None:
    """They live in the same module as the version, so a careless edit lands here too."""
    assert CONTRACT_SCHEMA_VERSION == 1
    assert (repo_root / "src" / "batchdx" / "cyto" / "_version.py").is_file()


# --------------------------------------------------------------------------------------
# 1.11  The corrected documents
# --------------------------------------------------------------------------------------


def test_the_notices_no_longer_deny_the_extracted_components(repo_root: Path) -> None:
    """Seven components are registered `compatibility` / `upstream_verified=True`."""
    notices = (repo_root / "THIRD_PARTY_NOTICES.md").read_text()
    assert "no implementation was extracted" not in notices, (
        "the notices still deny an extraction that happened; see docs/blocked.md B-007, "
        "which is RESOLVED"
    )
    for component in (
        "emd_mean_ct_horiz",
        "emd_mean_ct_vert",
        "ratio_consistent_peaks",
        "average_batch_r2_ct",
        "iLisi",
        "cLisi",
        "flowsom_mean_mapping_similarity",
    ):
        assert component in notices, f"the notices omit the verified component {component}"


def test_the_notices_list_the_metrics_extra_dependencies(repo_root: Path) -> None:
    """The Python dependency list omitted scikit-learn and scib-metrics entirely."""
    notices = (repo_root / "THIRD_PARTY_NOTICES.md").read_text()
    for package in ("scikit-learn", "scib-metrics"):
        assert package in notices, f"the notices omit the metrics-extra dependency {package}"


@requires_toml
def test_every_declared_python_dependency_is_named_in_the_notices(repo_root: Path) -> None:
    """The notices are the file a distributor reads; the list must not drift."""
    notices = (repo_root / "THIRD_PARTY_NOTICES.md").read_text().lower()
    core = _pyproject(repo_root)["project"]["dependencies"]
    for requirement in core:
        name = re.split(r"[<>=!~\[ ]", requirement, maxsplit=1)[0].strip().lower()
        assert name in notices, f"the notices omit the core dependency {name!r}"


def test_the_envs_readme_names_every_environment_file_next_to_it(repo_root: Path) -> None:
    """The README claimed one environment file while two sat beside it.

    Written as a walk over the directory rather than as a count, so adding a third file
    fails this test instead of quietly making the README stale again.
    """
    envs = repo_root / "workflow" / "envs"
    files = sorted(p.name for p in envs.glob("*.yaml"))
    assert files, "workflow/envs/ holds no environment file to describe"
    readme = (envs / "README.md").read_text()
    for name in files:
        assert name in readme, f"workflow/envs/README.md does not mention {name}"


def test_the_cytobatchflagr_dependency_risk_names_the_enforced_library_block(
    repo_root: Path,
) -> None:
    """The enforced set is the unconditional `library()` block, not a summary of it."""
    lock = yaml.safe_load((repo_root / "config" / "upstream_lock.yaml").read_text())
    risk = lock["resources"]["cytobatchflagr"]["dependency_risk"]
    assert "CytoBatchFlagR_functions.R" in risk
    assert UPSTREAM_LIBRARY_LINES in risk, (
        f"the risk must cite the anchor {UPSTREAM_LIBRARY_LINES} of the library block"
    )
    assert str(ENFORCED_LIBRARY_COUNT) in risk


def test_the_cytobatchflagr_dependency_risk_does_not_claim_test_coverage(
    repo_root: Path,
) -> None:
    """The risk must match the coverage that actually exists, in either direction.

    This test originally required the text "no cytobatchflagr r code executes", because
    when it was written both CytoBatchFlagR tests stopped at an eligibility guard before
    reaching `run_r_script`, and a document that implied otherwise would be over-claiming.

    CytoBatchFlagR R code now executes for real, so that assertion is inverted.  The
    check is kept in the direction that still matters: the document must not claim
    coverage it does not have.  It now requires the honest positive statement, and
    requires that the specific test file which provides the coverage is named, so the
    claim is traceable rather than bare.
    """
    risk = yaml.safe_load(
        (repo_root / "config" / "upstream_lock.yaml").read_text()
    )["resources"]["cytobatchflagr"]["dependency_risk"]
    assert "tests/integration/test_cytobatchflagr_backend.py" in risk, (
        "the risk text claims CytoBatchFlagR R code executes, so it must name the test "
        "file that executes it; a coverage claim with no traceable source is a claim"
    )
    # Matched case-insensitively: the risk text writes "DOES now execute", and a test
    # that only passes for one capitalisation of an emphatic word is a test that will
    # break the next time somebody tidies the prose.
    assert "does now execute" in risk.lower()
    # And it must not simultaneously claim the opposite.
    assert "NO CytoBatchFlagR R code executes" not in risk, (
        "the risk text now claims the R code executes and simultaneously says it does "
        "not; one of the two is stale"
    )


def test_the_cytobatchflagr_dependency_risk_enumerates_the_whole_enforced_set(
    repo_root: Path, upstream_dir: Path
) -> None:
    """Every enforced package must be named, not just the missing ones.

    A count in prose drifts from the list next to it: this file once said "the other
    twelve are CRAN" beside a list of twenty-two.  Comparing the prose count to the
    names it just enumerated catches that class of error, which a spot-check cannot.
    """
    enforced = _library_block_packages(upstream_dir)
    assert len(enforced) == ENFORCED_LIBRARY_COUNT
    risk = yaml.safe_load(
        (repo_root / "config" / "upstream_lock.yaml").read_text()
    )["resources"]["cytobatchflagr"]["dependency_risk"]
    for package in enforced:
        assert package in risk, f"the risk omits an unconditionally loaded package: {package}"

    named_bioc = set(re.search(r"five are Bioconductor \(([^)]*)\)", risk).group(1).replace(" ", "").split(","))
    assert named_bioc == (set(enforced) & {"flowCore", "FlowSOM", "limma", "ComplexHeatmap",
                                           "ConsensusClusterPlus"}), (
        f"the Bioconductor subset named in the risk does not match the enforced set: {named_bioc}"
    )
    cran_count = int(re.search(r"(\d+) are CRAN", risk).group(1))
    assert cran_count == len(set(enforced) - named_bioc - R_BUNDLED), (
        f"the risk claims {cran_count} CRAN packages but the enforced set has "
        f"{len(set(enforced) - named_bioc - R_BUNDLED)}"
    )

def test_the_cytobatchflagr_dependency_risk_names_what_is_actually_missing(
    repo_root: Path, upstream_dir: Path
) -> None:
    """The strongest form: measure the gap and require the document to account for it.

    This test originally asserted that ``missing`` was NON-EMPTY, so that the lock's
    ``dependency_risk`` would have to name every package absent from the R library.  That
    was the right check while six of the 29 enforced packages were unavailable.  It is now
    obsolete in that form: all 29 are installed, so a non-empty ``missing`` would only be
    satisfiable by breaking the environment.

    The check is inverted instead, and it is stricter rather than weaker.  If a package
    IS missing, the document must still name it -- so the check is kept alive for the
    environments where it matters.  And when nothing is missing, the document must say so
    and must not carry a stale claim that something is absent, which is the failure mode
    that would otherwise go unnoticed once the environment is complete.
    """

    installed = _r_library_names(repo_root)
    missing = sorted(set(_library_block_packages(upstream_dir)) - installed - R_BUNDLED)
    risk = yaml.safe_load(
        (repo_root / "config" / "upstream_lock.yaml").read_text()
    )["resources"]["cytobatchflagr"]["dependency_risk"]
    for package in missing:
        assert package in risk, (
            f"{package} is loaded unconditionally by upstream and is absent from the "
            "installed R library, but the dependency_risk does not name it"
        )
    if not missing:
        assert "ALL 29 are now installed" in risk, (
            "no enforced package is missing, so the dependency_risk must record that the "
            f"gap is closed. Actual risk text starts: {risk[:200]!r}"
        )
        assert "not installed" not in risk, (
            "the risk text still claims a package is not installed, but all 29 are"
        )


# --------------------------------------------------------------------------------------
# 1.12  The Snakefile
# --------------------------------------------------------------------------------------


def test_the_snakefile_imports_the_error_it_raises(repo_root: Path) -> None:
    """Snakemake injects `WorkflowError` into the Snakefile namespace, but that is an
    implementation detail of `snakemake.workflow.globals()`, not a documented API.

    Resolved empirically either way: a Snakefile containing nothing but
    `raise WorkflowError(...)` works under both `snakemake` and `snakemake --lint`. So
    this is not a fix for an observed crash -- it is an explicit import in place of an
    undocumented injection, and the test keeps it that way.
    """
    snakefile = (repo_root / "workflow" / "Snakefile").read_text()
    assert "from snakemake.exceptions import WorkflowError" in snakefile
    assert "raise WorkflowError(" in snakefile, (
        "the import would be dead weight if nothing raises it"
    )


# --------------------------------------------------------------------------------------
# 1.13  The container, which has never been built
# --------------------------------------------------------------------------------------


def test_the_container_base_image_is_digest_pinned(repo_root: Path) -> None:
    """A tag is a mutable pointer; a digest is the image."""
    dockerfile = (repo_root / "containers" / "Dockerfile.r-backends").read_text()
    from_lines = [line for line in dockerfile.splitlines() if line.startswith("FROM ")]
    assert from_lines, "the container defines no FROM"
    for line in from_lines:
        assert re.search(r"@\s*sha256:[0-9a-f]{64}\s*$", line), (
            f"the base image is not digest-pinned: {line!r}"
        )


def test_the_container_banner_does_not_claim_it_has_been_built(repo_root: Path) -> None:
    """No container runtime exists on the node this was written on."""
    dockerfile = (repo_root / "containers" / "Dockerfile.r-backends").read_text()
    banner = dockerfile.split("FROM ", 1)[0]
    assert "NEVER BUILT" in banner, "the banner must say the image has never been built"
    for overclaim in ("verified image", "tested image", "reproducible build"):
        assert overclaim not in banner.lower(), f"the banner overclaims: {overclaim!r}"


def test_the_container_installs_everything_the_upstream_library_block_requires(
    repo_root: Path, upstream_dir: Path
) -> None:
    """`library()` is unconditional, so one missing package aborts the whole source().

    Re-derived from the pinned checkout rather than from this repository's prose, so
    the image cannot quietly fall behind what upstream actually loads.
    """
    required = set(_library_block_packages(upstream_dir))
    assert len(required) == ENFORCED_LIBRARY_COUNT
    cran, bioc = _install_lists(repo_root)
    installed = cran | bioc | R_BUNDLED
    missing = sorted(required - installed)
    assert not missing, (
        "the container does not install these packages that the upstream library block "
        f"loads unconditionally: {missing}"
    )


def _scannable_files(repo_root: Path) -> list[Path]:
    """Code and configuration, from the roots where an environment read could live.

    Bounded on purpose. A whole-tree walk here is both slow and wrong: it reaches the
    cookiecutter template and ``pixi.lock``, and a two-letter variable name such as
    ``PATH`` then matches in dozens of unrelated places.
    """
    roots = ("src", "workflow", "scripts", "hooks", "config", "tests", ".github")
    this_file = Path(__file__).resolve()
    found: list[Path] = []
    for root in roots:
        base = repo_root / root
        if not base.is_dir():
            continue
        found += [
            path
            for path in base.rglob("*")
            if path.is_file()
            and path.resolve() != this_file
            and "__pycache__" not in path.parts
            and (path.suffix in _CODE_SUFFIXES or path.name == "Snakefile")
        ]
    return sorted(found)


def test_the_container_declares_no_environment_variable_the_repository_ignores(
    repo_root: Path,
) -> None:
    """`BATCHDX_UPSTREAM_DIR` was declared and read by nothing but its own line.

    Prose is excluded, and so is this file: a test that mentions the variable it is
    hunting for would otherwise satisfy itself, which is how a check quietly stops
    being able to fail.
    """
    dockerfile = repo_root / "containers" / "Dockerfile.r-backends"
    declared = re.findall(
        r"^ENV\s+([A-Za-z_][A-Za-z0-9_]*)", dockerfile.read_text(), re.MULTILINE
    )
    scanned = _scannable_files(repo_root)
    assert scanned, "the scan found no candidate files, so it could not have failed"
    for name in declared:
        hits = [
            path.relative_to(repo_root).as_posix()
            for path in scanned
            if name in path.read_text(errors="ignore")
        ]
        assert hits, (
            f"the container sets {name} and nothing reads it. An environment variable "
            "that no code or configuration consults is a promise the repository does not "
            "keep; wire it or delete it."
        )


# --------------------------------------------------------------------------------------
# 1.14  The compatibility matrix
# --------------------------------------------------------------------------------------


def test_the_compatibility_matrix_exists_and_states_its_rule(repo_root: Path) -> None:
    """A matrix without a rule for what `verified` means is an opinion."""
    matrix = repo_root / "docs" / "compatibility_matrix.md"
    assert matrix.is_file()
    text = matrix.read_text()
    for status in ("verified", "declared-only", "untested"):
        assert status in text, f"the matrix never uses the status {status!r}"


def test_the_compatibility_matrix_does_not_claim_an_unrun_python_version_was_run(
    repo_root: Path,
) -> None:
    """3.10, 3.11 and 3.13 are declared here and have never been executed from this
    checkout. Only 3.12 has. The matrix must keep saying so, per row."""
    text = (repo_root / "docs" / "compatibility_matrix.md").read_text()
    for version in ("3.10", "3.11", "3.13"):
        row = next(
            (
                line
                for line in text.splitlines()
                if line.startswith(f"| {version} ") and "|" in line
            ),
            None,
        )
        assert row is not None, f"the matrix has no row for Python {version}"
        assert "`verified`" not in row, (
            f"the matrix claims Python {version} was verified; it has never been executed "
            f"from this checkout -- row reads: {row}"
        )
