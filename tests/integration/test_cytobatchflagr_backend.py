"""The CytoBatchFlagR happy path: upstream R code that actually executes.

Until this file existed, no CytoBatchFlagR R code had ever run in this repository. Both
pre-existing tests of the adapter stopped at a guard before ``run_r_script``: one had no
checkout, the other had no thresholds, and the shared ``frozen`` fixture has three
batches, so the adapter returned ``ineligible`` at ``min_batches = 4``. The
``requires_cytobatchflagr`` marker was therefore declared and used by zero tests.

What is exercised here, and why each part is load-bearing:

* **The design has to be its own.** Four batches with a measured technical replicate in
  every one, because the IQR rule is degenerate below four. ``simulate.py`` has no
  four-batch scenario, so the spec is built here from ``SampleSpec``/``SimulationSpec``
  on the pattern of ``tests/scientific/test_stage3_design.py:150-154``.
* **A real ``threshold_path``.** The adapter refuses to pick positivity thresholds itself
  (five analysts chose five different thresholds for the same markers), so without one it
  never reaches R at all. The panel is CD3 and CD4 at 1.75, which are the shipped
  example's own declared gating thresholds.
* **Why the panel is two markers and not six.** Upstream's helpers drop any control with
  fewer than 100 cells on one side of a threshold, silently, and ``min_cells`` is
  disabled by the wrapper so the shortfall can be reported instead. At 400 cells per
  sample only CD3 and CD4 clear 100 cells on *both* sides of 1.75 in every batch: CD14,
  CD19, CD56 and CD8 leave 30-92 cells above the threshold, which the adapter now refuses
  to evaluate. The cell count was not raised to dodge this, and the counts are asserted
  in ``test_cell_counts_justify_the_panel`` so the choice cannot rot silently.
* **A negative control, and a candid one.** The same design is run with and without an
  injected shift, on the same seed, so every cell is identical except the shift. The
  clean run produces no flags on this fixture, but that is *not* evidence that the rule
  is specific: ``detect_outlier`` applies a 1.5 IQR Tukey fence to four values, and the
  measured probability that at least one of four normally distributed values falls outside
  that fence is 20.7% (5.15% per value, 400,000 draws). The 1.5 IQR rule assumes a large
  sample; at four batches it has roughly a one-in-five chance of firing per test, so a
  panel of tests is expected to flag something on clean data and the number does not
  improve with more batches (27% at six, 27% at eight, 30% at twenty). Six tests here
  means a clean run flags nothing about one time in four. The test therefore asserts the
  *differential* -- the marker under test is flagged in the shifted run and not in the
  clean run -- which no amount of fence noise can satisfy by accident, and records the
  rate so the clean result is not read as specificity.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from batchdx.cyto.adapters.r_runner import r_package_version
from batchdx.cyto.config import load_config
from batchdx.cyto.design import audit_design
from batchdx.cyto.diagnostics import DetectorContext, get_detector
from batchdx.cyto.diagnostics.cytobatchflagr import MIN_CELLS_PER_CONTROL_MARKER
from batchdx.cyto.examples import materialise_example
from batchdx.cyto.io import load_dataset
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation

pytestmark = pytest.mark.requires_cytobatchflagr

#: config/upstream_lock.yaml:94. Asserted against the lock as well, so moving the pin
#: without moving this test is a failure rather than a silently weaker assertion.
PINNED_COMMIT = "ad21e3b1"

BATCHES = ("B1", "B2", "B3", "B4")
SHIFTED_BATCH = "B4"
SHIFTED_MARKER = "CD3"
UNSHIFTED_MARKER = "CD4"
REFERENCE_ID = "REF1"
CELLS_PER_SAMPLE = 400
PANEL: tuple[tuple[str, float], ...] = (("CD3", 1.75), ("CD4", 1.75))
METRICS = ("-MFI", "+MFI", "%pos")

#: Ship with R itself, so they are present in any R installation rather than installed.
_BASE_R_PACKAGES = frozenset({"stats", "MASS", "graphics", "grDevices", "utils", "methods"})

#: Measured on 400,000 draws, not assumed. See the module docstring.
FENCE_FALSE_POSITIVE_RATE_AT_FOUR_BATCHES = 0.207


def _require_checkout(upstream_dir: Path) -> Path:
    functions = upstream_dir / "CytoBatchFlagR" / "R" / "CytoBatchFlagR_functions.R"
    if not functions.exists():
        pytest.skip(
            "no CytoBatchFlagR checkout at the pinned commit; run "
            "scripts/fetch_upstream.sh. This is the only condition under which these tests "
            "skip: an R library that cannot satisfy upstream's unconditional library() "
            "block must fail these tests, not skip them, or the R layer goes unexercised "
            "again in exactly the way it had been."
        )
    return functions


def _spec(shift: dict[str, dict[str, float]] | None) -> SimulationSpec:
    """Four batches, one measured reference aliquot in each, balanced case/control."""

    samples: list[SampleSpec] = []
    for b_i, batch in enumerate(BATCHES):
        samples.append(SampleSpec(
            sample_id=f"ref_{batch}", donor_id="REFDONOR", batch=batch, group="reference",
            reference_id=REFERENCE_ID, reference_role="technical_replicate",
            n_cells=CELLS_PER_SAMPLE,
        ))
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                donor = f"D{b_i}{k}{tag}"
                samples.append(SampleSpec(
                    sample_id=f"{donor}_{batch}", donor_id=donor, batch=batch, group=group,
                    n_cells=CELLS_PER_SAMPLE,
                ))
    return SimulationSpec(samples=samples, seed=473, batch_shift=shift or {})


def _build(tmp_path: Path, shift: dict[str, dict[str, float]] | None,
           panel: tuple[tuple[str, float], ...] = PANEL):
    """A materialised four-batch project, loaded, with a real detector threshold table."""

    project = tmp_path / "project"
    materialise_example("paired_reference", project, cells_per_sample=CELLS_PER_SAMPLE)
    write_simulation(_spec(shift), project)

    # A separate file: config/thresholds.tsv is the label source and needs the full
    # gating schema (population/sign/priority), so it cannot double as the detector's
    # fixed positivity table.
    thresholds = project / "config" / "iqr_thresholds.tsv"
    thresholds.write_text(
        pd.DataFrame(panel, columns=["feature_id", "threshold"]).to_csv(sep="\t", index=False),
        encoding="utf-8",
    )

    config = load_config(project / "config" / "config.yaml")
    dataset = load_dataset(config, project)
    design = audit_design(
        list(dataset.samples.values()), dataset_id="cbf_happy",
        correction_factor="batch", protected=["group"], require_shared_reference=True,
        minimum_batches=3,
    )
    labels, tier, confidence = resolve_labels(dataset, config, project)
    sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
    return project, dataset, design, sampling, thresholds


def _run(dataset, design, sampling, thresholds, out_dir, upstream_dir, *, panel=PANEL):
    context = DetectorContext(
        run_id="t", design=design, sampling=sampling, out_dir=str(out_dir), seed=473,
        upstream_dir=str(upstream_dir), threshold_path=str(thresholds),
    )
    return get_detector("cytobatchflagr").run(dataset, context)


def _flagged(result, *, marker: str | None = None) -> list:
    return [a for a in result.assessments
            if a.flag and (marker is None or a.marker == marker)]


def _cell_counts(result) -> pd.DataFrame:
    path = result.tables.get("cell_counts")
    assert path, "the wrapper must publish the counts that justify the panel it evaluated"
    return pd.read_csv(path)


# --------------------------------------------------------------------------------------
# The happy path
# --------------------------------------------------------------------------------------


def test_cytobatchflagr_flags_the_injected_shift_and_not_merely_a_batch(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """Four batches, a control in each, one batch shifted: only that batch is flagged."""

    _require_checkout(upstream_dir)
    _, dataset, design, sampling, thresholds = _build(
        tmp_path, {SHIFTED_BATCH: {SHIFTED_MARKER: 1.2}}
    )
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", upstream_dir)

    assert result.status == "ok", f"expected upstream R to run and answer; got {result.status}: {result.reason}"
    assert result.eligibility is not None and result.eligibility.eligible
    # every metric x marker x batch, so the fence really was given four values per group
    assert len(result.assessments) == len(METRICS) * len(PANEL) * len(BATCHES)
    assert {a.statistic for a in result.assessments} == set(METRICS)
    assert all(a.evidence_class == "technical_reference" for a in result.assessments)
    assert all(a.threshold == 1.5 for a in result.assessments)

    flags = _flagged(result)
    assert flags, "an injected +1.2 shift on CD3 in one of four batches must trip the fence"
    # The shifted batch, not merely A batch: every flag is on B4, and B4 is flagged for
    # the shifted marker on all three metrics, and for no other marker.
    assert {a.batch for a in flags} == {SHIFTED_BATCH}
    assert {(a.marker, a.statistic) for a in flags} == {
        (SHIFTED_MARKER, m) for m in METRICS
    }
    assert not _flagged(result, marker=UNSHIFTED_MARKER)
    # The unshifted batches of the shifted marker are present and unflagged, which is what
    # makes "B4 and only B4" a statement about the shift rather than about the panel.
    cd3 = {(a.batch, a.flag) for a in result.assessments if a.marker == SHIFTED_MARKER}
    assert cd3 == {(b, b == SHIFTED_BATCH) for b in BATCHES}

    # The shift moved the value, and the value moved in the same direction for all three
    # metrics: not a single row that happens to sit on the other side of its own fence.
    by_batch = pd.DataFrame(
        [(a.statistic, a.batch, a.value) for a in result.assessments
         if a.marker == SHIFTED_MARKER],
        columns=["statistic", "batch", "value"],
    )
    pivot = by_batch.pivot(index="statistic", columns="batch", values="value")
    for statistic, row in pivot.iterrows():
        others = row.drop(SHIFTED_BATCH)
        assert row[SHIFTED_BATCH] > others.max(), (
            f"{statistic}: the shifted batch should be the extreme one, got "
            f"{row[SHIFTED_BATCH]} against an unshifted maximum of {others.max()}"
        )
        # and the three unshifted batches agree with each other to within the sampling
        # noise of a median over ~200 cells, so the shift is not the only thing that moved
        assert others.max() - others.min() < row[SHIFTED_BATCH] - others.max(), (
            f"{statistic}: the spread among unshifted batches "
            f"({others.max() - others.min()}) rivals the shift "
            f"({row[SHIFTED_BATCH] - others.max()}), so the fence would be flagging noise"
        )


def test_cytobatchflagr_clean_design_does_not_flag_the_marker_under_test(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """The negative control: the same design, the same seed, no injected shift.

    Paired with ``test_cytobatchflagr_flags_the_injected_shift_and_not_merely_a_batch``:
    the only difference between the two datasets is the +1.2 on CD3 in B4, so every cell
    is identical elsewhere. The differential -- flagged with the shift, not flagged
    without it, on the same marker -- is the claim; see the module docstring for why
    "zero flags" is not by itself evidence of specificity at four batches.
    """

    _require_checkout(upstream_dir)
    _, dataset, design, sampling, thresholds = _build(tmp_path, None)
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", upstream_dir)

    assert result.status == "ok", f"expected upstream R to run and answer; got {result.status}: {result.reason}"
    assert len(result.assessments) == len(METRICS) * len(PANEL) * len(BATCHES)
    assert not _flagged(result, marker=SHIFTED_MARKER), (
        "the marker that is shifted in the paired run must not be flagged here, or the "
        "shifted run's flag proves nothing about the shift"
    )
    assert not _flagged(result), (
        f"a clean design flagged {[a.assessment_id for a in _flagged(result)]}. The 1.5 IQR "
        f"fence fires on ~{FENCE_FALSE_POSITIVE_RATE_AT_FOUR_BATCHES:.0%} of four-value tests "
        f"even with no batch effect, so this is the expected failure mode of the rule rather "
        f"than a detector bug -- but the differential assertion above is the one that carries "
        f"the result, and it still holds if this does not."
    )


def test_cytobatchflagr_pins_the_upstream_commit_and_publishes_its_table(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """The arm's-length pin has to be observable, or it is not a pin."""

    _require_checkout(upstream_dir)
    _, dataset, design, sampling, thresholds = _build(
        tmp_path, {SHIFTED_BATCH: {SHIFTED_MARKER: 1.2}}
    )
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", upstream_dir)
    assert result.status == "ok"

    cytobatchflagr = [s for s in result.software if s.name == "CytoBatchFlagR"]
    assert len(cytobatchflagr) == 1
    recorded = cytobatchflagr[0].commit
    assert recorded, "the commit the checkout was actually read at must be recorded"
    assert recorded.startswith(PINNED_COMMIT), (
        f"expected the checkout pinned in config/upstream_lock.yaml ({PINNED_COMMIT}), "
        f"got {recorded!r}"
    )

    # The lock is the authority; asserting only the literal would let the pin drift.
    lock = upstream_dir.parent / "config" / "upstream_lock.yaml"
    assert lock.exists()
    assert PINNED_COMMIT in lock.read_text(encoding="utf-8")

    iqr = result.tables.get("iqr")
    assert iqr, "result.tables['iqr'] must name the table the flags were read from"
    assert Path(iqr).exists()
    frame = pd.read_csv(iqr)
    assert len(frame) == len(result.assessments)
    assert set(frame["batch"].astype(str)) == set(BATCHES)
    assert frame.loc[frame["batch"].astype(str) == SHIFTED_BATCH, "flagged_batch"].notna().all()
    assert frame["flagged_batch"].astype(str).eq("none").sum() == (
        len(frame) - len(_flagged(result))
    )


def test_cell_counts_justify_the_panel_that_was_evaluated(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """Per-batch, per-control cell counts, so a future silent drop cannot go unnoticed.

    The point of asserting the counts rather than the outcome is that a regression which
    reintroduces a dropped control would leave the *shifted* test vacuously green: fewer
    than four values means the fence cannot fire, and a fence that cannot fire produces
    no flag and no error. Asserting the inputs closes that door.
    """

    _require_checkout(upstream_dir)
    _, dataset, design, sampling, thresholds = _build(
        tmp_path, {SHIFTED_BATCH: {SHIFTED_MARKER: 1.2}}
    )
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", upstream_dir)
    assert result.status == "ok"

    counts = _cell_counts(result)
    assert result.requirements["min_cells_per_control_marker"] == MIN_CELLS_PER_CONTROL_MARKER
    assert set(counts["control"].astype(str)) == {REFERENCE_ID}
    assert set(counts["batch"].astype(str)) == set(BATCHES)
    assert set(counts["marker"].astype(str)) == {m for m, _ in PANEL}
    assert len(counts) == len(BATCHES) * len(PANEL)

    assert (counts["n_cells"] == CELLS_PER_SAMPLE).all()
    for column in ("n_negative", "n_positive"):
        worst = int(counts[column].min())
        assert worst >= MIN_CELLS_PER_CONTROL_MARKER, (
            f"{column} falls to {worst} cells, under the "
            f"{MIN_CELLS_PER_CONTROL_MARKER} the adapter requires; with that control dropped "
            f"the fence would be evaluated on three values and could not fire"
        )


def test_a_control_too_small_for_its_threshold_is_ineligible_not_clean(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """The vacuity guard, exercised on the failure it exists to prevent.

    Same 400-cell data, but CD3's threshold moved to 4.0: almost no cell in a control
    lands above it, so the positive-MFI median would be built on a handful of cells and
    upstream's own helpers would drop the control outright (functions.R:1084). The result
    must be a stated ``ineligible`` naming the metric, the marker and the batches. Were it
    ``ok`` with no flags, a batch effect would have been reported as an absence of one --
    and the test suite would have called that a pass.
    """

    _require_checkout(upstream_dir)
    panel = (("CD3", 4.0), ("CD4", 1.75))
    _, dataset, design, sampling, thresholds = _build(tmp_path, None, panel=panel)
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", upstream_dir)

    assert result.status == "ineligible", (
        f"a control that cannot support the metric must be ineligible, not clean; got "
        f"{result.status}: {result.reason}"
    )
    assert not result.assessments, "an ineligible detector produces no flags at all"
    assert result.requirements["min_cells_per_control_marker"] == MIN_CELLS_PER_CONTROL_MARKER
    assert "+MFI" in result.reason
    assert SHIFTED_MARKER in result.reason
    assert REFERENCE_ID in result.reason
    assert "min_cells" in result.reason
    for batch in BATCHES:
        assert batch in result.reason, "the reason must name the batches, not just the count"
    # And the shortfall is real, not asserted: the counts behind it are published.
    counts = _cell_counts(result)
    cd3_positive = counts.loc[counts["marker"] == SHIFTED_MARKER, "n_positive"]
    assert (cd3_positive < MIN_CELLS_PER_CONTROL_MARKER).all()


# --------------------------------------------------------------------------------------
# The preflight over upstream's unconditional library() block
# --------------------------------------------------------------------------------------


def test_every_package_in_the_upstream_library_block_is_installed(upstream_dir: Path) -> None:
    """The binding dependency set, checked from the file rather than from a list.

    Upstream ships no DESCRIPTION and no NAMESPACE, so the ``library()`` calls are the
    only machine-readable statement of what it needs, and they are unconditional. A
    package that is not installed aborts the whole ``source()``. Parsing them here means
    an added or renamed ``library()`` call fails this test until the environment is
    extended, instead of aborting a detection run at runtime.
    """

    _require_checkout(upstream_dir)
    text = (upstream_dir / "CytoBatchFlagR/R/CytoBatchFlagR_functions.R").read_text(
        encoding="utf-8"
    )
    packages = sorted(set(re.findall(r"library\(\s*([A-Za-z][A-Za-z0-9._]*)", text)))
    assert len(packages) == 29, f"the library() block changed: {packages}"
    assert "ComplexHeatmap" in packages and "patchwork" in packages

    missing = [
        p for p in packages
        if not r_package_version(p).available and p not in _BASE_R_PACKAGES
    ]
    assert not missing, (
        f"upstream's unconditional library() block needs {missing}, which are not installed. "
        f"library() is not wrapped in require(), so one absent package aborts the entire "
        f"source(); build them with scripts/setup_r_env.sh."
    )


def test_a_missing_r_package_is_not_available_with_a_stated_remedy(
    tmp_path: Path, upstream_dir: Path
) -> None:
    """A missing package must be an unavailability, not a failed detection.

    The environment this test runs in has all 29 packages installed, so the negative
    branch of the preflight is reached with a stand-in checkout: a file that is not
    upstream and contains no upstream code, naming one package that genuinely is not
    installed. The wrapper returns before it would ever source it, so nothing upstream is
    involved. Without this, the preflight's failure path is only reachable on a machine
    that is missing a package, which is what left it unexercised in the first place.
    """

    _require_checkout(upstream_dir)
    stand_in = tmp_path / "standin" / "CytoBatchFlagR" / "R"
    stand_in.mkdir(parents=True)
    absent = "a.package.that.is.not.installed"
    assert not r_package_version(absent).available
    (stand_in / "CytoBatchFlagR_functions.R").write_text(
        "## NOT upstream source. A test double for the preflight's negative branch only.\n"
        f"library(dplyr)\nlibrary({absent})\n",
        encoding="utf-8",
    )

    _, dataset, design, sampling, thresholds = _build(tmp_path, None)
    result = _run(dataset, design, sampling, thresholds, tmp_path / "out", tmp_path / "standin")

    assert result.status == "not_available", (
        f"an absent R package must be reported as unavailability, not as a failure; got "
        f"{result.status}: {result.reason}"
    )
    assert not result.assessments
    assert absent in result.reason
    assert "setup_r_env.sh" in result.reason
    assert "never reported as an absence of batch effects" in result.reason
    assert result.requirements["missing_packages"] == [absent]


def test_no_r_interpreter_is_not_available_never_an_absence_of_batch_effects(
    tmp_path: Path, upstream_dir: Path
) -> None:
    _, dataset, design, sampling, thresholds = _build(
        tmp_path, {SHIFTED_BATCH: {SHIFTED_MARKER: 1.2}}
    )
    context = DetectorContext(
        run_id="t", design=design, sampling=sampling, out_dir=str(tmp_path / "out"),
        seed=473, upstream_dir=str(upstream_dir), threshold_path=str(thresholds),
        rscript="Rscript-does-not-exist",
    )
    result = get_detector("cytobatchflagr").run(dataset, context)
    assert result.status == "not_available"
    assert not result.assessments
    assert "not on PATH" in result.reason
