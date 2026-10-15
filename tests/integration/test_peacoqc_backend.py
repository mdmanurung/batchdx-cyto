"""The PeacoQC backend: the real call, the removal index, and the unavailable path.

Three groups, in increasing order of what they need from the environment.

**No R.**  The refusals.  A channel PeacoQC was not asked to evaluate, an
``rscript`` that does not exist, an identity-channel name collision, misaligned rows.
None of these can be checked by *having* PeacoQC, which is why they are here rather than
skipped along with the rest: a path that only exists when the dependency is installed is
a path that is never executed on the machines that need it most.

**R, no PeacoQC.**  The unavailable path, driven with
``rscript="Rscript-does-not-exist"`` -- the pattern that made the Stage 5 unavailable
path testable in every environment.  It needs no R at all.

**PeacoQC installed** (``requires_r``).  Real calls, against a synthetic fixture with a
*known* contaminated block, asserting that the block is removed and that the removal is
attributed.  The environment-dependent tests skip *visibly* and say why.  They are not
marked xfail: a skipped test is a recorded fact, whereas an xfail would let a broken
implementation look like an expected outcome.

The FlowVS test in the last group is the load-bearing one for the transform convention.
It does not read ``flowVS``'s documentation and compare it to a comment; it runs
``flowVS::transFlowVS`` on the same numbers and compares element-wise.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml

from batchdx.cyto.adapters.r_runner import r_package_version, script_path
from batchdx.cyto.io.fcs import read_fcs, write_fcs
from batchdx.cyto.preprocess import (
    IDENTITY_CHANNEL,
    VERIFIED_FORMALS,
    VERIFIED_VERSION,
    PeacoQCError,
    PeacoQCNotAvailable,
    PeacoQCSettings,
    asinh_invert,
    asinh_transform,
    read_gates,
    read_peacoqc_settings,
    run_peacoqc,
)

#: The name used to make "there is no R here" testable without uninstalling R.
NO_SUCH_RSCRIPT = "Rscript-does-not-exist"

#: The channels a fixture carries.  Four are scatter/viability parameters of the kind a
#: preliminary gate uses and two are panel markers, so a test that accidentally gated
#: instead of QCing would be visible.
CHANNELS = ["FSC-A", "SSC-A", "FSC-H", "FVS780", "CD4", "CD3"]

requires_r = pytest.mark.requires_r


def _fixture(n: int, block: tuple[int, int], seed: int = 473) -> np.ndarray:
    """A synthetic sample with a contiguous, *known* contaminated acquisition block.

    The contamination is a contiguous range rather than scattered indices because that
    is what a clog or a flow-rate change looks like, and because it gives the test a
    ground truth it can assert against: every event in ``block`` is removed, and the
    count is known exactly.
    """

    rng = np.random.default_rng(seed)
    values = np.column_stack(
        [
            np.abs(rng.normal(150_000, 18_000, n)),
            np.abs(rng.normal(60_000, 11_000, n)),
            np.abs(rng.normal(140_000, 17_000, n)),
            np.abs(rng.normal(1_800, 700, n)),
            np.abs(rng.normal(2_500, 700, n)),
            np.abs(rng.normal(1_800, 600, n)),
        ]
    )
    lo, hi = block
    values[lo:hi, 4] = values[lo:hi, 4] * 25 + 30_000
    values[lo:hi, 5] = values[lo:hi, 5] * 20 + 25_000
    debris = rng.choice(n, max(1, n // 12), replace=False)
    values[debris, 4] *= 12
    return values


def _call(tmp_path: Path, settings: PeacoQCSettings, n: int, block: tuple[int, int], **kwargs):
    values = _fixture(n, block)
    return run_peacoqc(
        values, CHANNELS,
        sample_ids=np.repeat("s1", n),
        original_event_index=np.arange(n),
        settings=settings,
        work_dir=tmp_path / "work",
        **kwargs,
    )


def _need_peacoqc() -> None:
    version = r_package_version("PeacoQC", "Rscript")
    if not version.available:
        pytest.skip(
            "the PeacoQC R package is not installed in this environment; run "
            f"scripts/setup_r_env.sh, which installs bioconductor-peacoqc ({version.detail})"
        )


# --------------------------------------------------------------------------------------
# No R: the refusals
# --------------------------------------------------------------------------------------


def test_a_missing_rscript_is_a_stated_unavailability_not_a_silent_zero(tmp_path: Path) -> None:
    """The pattern that makes an unavailable path testable in every environment.

    A run that could not execute must not come back saying "no events were removed",
    which is indistinguishable from a run that executed and found a clean sample.
    """

    with pytest.raises(PeacoQCNotAvailable) as excinfo:
        _call(tmp_path, PeacoQCSettings(channels=("CD4",)), n=200, block=(50, 80),
              rscript=NO_SUCH_RSCRIPT)
    message = str(excinfo.value)
    assert NO_SUCH_RSCRIPT in message
    assert "does not substitute another quality-control algorithm" in message


def test_row_misalignment_is_refused_rather_than_silently_truncated(tmp_path: Path) -> None:
    """A sample vector shorter than the matrix, which every downstream index would hide."""

    with pytest.raises(PeacoQCError, match="row-aligned"):
        run_peacoqc(
            _fixture(200, (50, 80)), CHANNELS,
            sample_ids=np.repeat("s1", 199),
            original_event_index=np.arange(200),
            settings=PeacoQCSettings(channels=("CD4",)),
            work_dir=tmp_path / "work",
            rscript=NO_SUCH_RSCRIPT,
        )


def test_the_reserved_identity_channel_cannot_collide_with_an_evaluated_channel(
    tmp_path: Path,
) -> None:
    """The identity channel is what makes a removal index resolvable, so it is reserved."""

    channels = [*CHANNELS, IDENTITY_CHANNEL]
    values = np.column_stack([_fixture(200, (50, 80)), np.arange(200, dtype=float)])
    with pytest.raises(PeacoQCError, match="reserved for the event-identity check"):
        run_peacoqc(
            values, channels,
            sample_ids=np.repeat("s1", 200),
            original_event_index=np.arange(200),
            settings=PeacoQCSettings(channels=(IDENTITY_CHANNEL,)),
            work_dir=tmp_path / "work",
            rscript=NO_SUCH_RSCRIPT,
        )


def test_the_r_script_exists_where_run_r_script_looks_for_it() -> None:
    """A packaged script that is not where the runner looks for it is a silent failure."""

    assert script_path("peacoqc_qc.R").exists()
    assert script_path("peacoqc_qc.R").name == "peacoqc_qc.R"


def test_the_verified_formals_table_records_the_differences_from_the_reference_call_site() -> None:
    """The recorded signature is a reading, so the assumptions built on it are pinned.

    A table that had drifted from the package it was read off would make every claim in
    ``preprocess/peacoqc.py``'s docstring stale, and nothing else in this build would
    notice.  Checked here against Python so it holds with or without R.
    """

    assert VERIFIED_FORMALS["determine_good_cells"] == "all"
    assert VERIFIED_FORMALS["save_fcs"] is True, (
        "upstream's default is TRUE and this adapter passes FALSE; if the package ever "
        "changed that default the difference would stop being this adapter's choice"
    )
    assert VERIFIED_FORMALS["plot"] == 20, (
        "upstream's default is the number 20 (plot once >=20% of measurements removed), "
        "not TRUE; the reference call site passes TRUE and this adapter passes False"
    )
    assert VERIFIED_FORMALS["MAD"] == 6 and VERIFIED_FORMALS["IT_limit"] == 0.6
    assert VERIFIED_FORMALS["min_cells"] == 150 and VERIFIED_FORMALS["step"] == 500
    assert VERIFIED_FORMALS["time_channel_parameter"] == "Time"
    assert "BadCell" not in VERIFIED_FORMALS, (
        "there is no BadCell argument and no BadCell result element; the removed set is "
        "!GoodCells"
    )
    assert "GoodCells" not in VERIFIED_FORMALS, "GoodCells is a RESULT element, not an argument"


# --------------------------------------------------------------------------------------
# PeacoQC installed
# --------------------------------------------------------------------------------------


@requires_r
def test_peacoqc_is_installed_at_the_version_whose_formals_were_read() -> None:
    _need_peacoqc()
    version = r_package_version("PeacoQC", "Rscript")
    assert version.version == VERIFIED_VERSION, (
        f"PeacoQC {version.version} is installed but VERIFIED_FORMALS was read from "
        f"{VERIFIED_VERSION}. Re-read the formals before trusting the table"
    )


@requires_r
def test_flowvs_transflowvs_agrees_with_this_build_elementwise(tmp_path: Path) -> None:
    """Item 3.2, executed rather than argued.

    ``flowVS::transFlowVS`` is run on a synthetic flowFrame and its output compared with
    ``preprocess.transform.asinh_transform`` element by element, and ``raw`` recovered
    with ``preprocess.transform.asinh_invert``.  If the two conventions ever diverged --
    if flowVS multiplied by the cofactor on the way out, or this build did -- this fails
    with the offending numbers rather than leaving the disagreement in a docstring.

    The scaled convention is computed here too and asserted to be *different*, because
    the exporter's exact-inverse claim at ``io/export.py:47-52`` is true only for the
    unscaled form and nothing else in the suite would notice if that changed.
    """

    if not r_package_version("flowVS", "Rscript").available:
        pytest.skip(
            "the flowVS R package is not installed. No script calls it; it is present only "
            "so this cross-check can run. scripts/setup_r_env.sh installs it"
        )
    raw = [0.0, 100.0, 2500.0, 50000.0, 250000.0, 1000000.0]
    path = tmp_path / "probe.fcs"
    write_fcs(
        path, np.array(raw, dtype=np.float64).reshape(1, 6), CHANNELS,
        scale="raw instrument units",
        processing_note="transform-convention cross-check fixture", overwrite=True,
    )
    script = tmp_path / "probe.R"
    script.write_text(
        "suppressPackageStartupMessages({library(flowCore); library(flowVS)})\n"
        "a <- commandArgs(trailingOnly = TRUE)\n"
        "ff <- read.FCS(a[1], transformation = FALSE, truncate_max_range = FALSE)\n"
        "fs <- flowSet(list(ff))\n"
        "out <- transFlowVS(fs, channels = colnames(fs), cofactors = rep(5, 6))\n"
        "cat(jsonlite::toJSON(as.numeric(exprs(out[[1]])[1, ]), digits = 17), '\\n')\n",
        encoding="utf-8",
    )
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(script), str(path)], capture_output=True, text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    theirs = np.array(json.loads(proc.stdout.strip().splitlines()[-1]), dtype=np.float64)
    ours = asinh_transform(np.array(raw, dtype=np.float64), 5.0)
    assert theirs.shape == ours.shape
    # The bound is the medium's, not the formula's. flowCore holds expression values as
    # single-precision doubles, so flowVS's output is float32(ours) and the largest
    # relative difference is a few times 1e-8. Anything larger would be a different
    # formula, and anything at the 1e-12 level would mean the two agreed bit for bit --
    # which they cannot through a float32 store.
    assert np.allclose(theirs, ours, rtol=1e-6, atol=1e-6), (
        "flowVS::transFlowVS and batchdx.cyto.preprocess.asinh_transform disagree by more "
        f"than single precision allows. Theirs: {theirs[:3]}. Ours: {ours[:3]}. flowVS "
        "divides by the cofactor and does NOT multiply by it on the way out, so "
        "io.export's sinh(x) * cofactor IS its exact inverse"
    )
    assert np.allclose(theirs, ours.astype(np.float32).astype(np.float64), rtol=0, atol=0), (
        "flowVS's output is not float32(asinh(raw / cofactor)). It is a different function "
        "of the input, and io.export's inverse would then be wrong"
    )
    assert not np.allclose(theirs, ours * 5.0, rtol=1e-3), (
        "flowVS's output is the SCALED asinh. Then this build's convention is wrong, "
        "io.export's inverse is wrong, and both must change together"
    )
    assert asinh_invert(theirs, 5.0) == pytest.approx(raw, rel=1e-5, abs=1.0), (
        "raw was not recovered from flowVS's own output by sinh(x) * cofactor"
    )


@requires_r
def test_peacoqc_removes_the_known_contaminated_block_and_attributes_it(
    tmp_path: Path,
) -> None:
    """A real PeacoQC call, and a ground truth it has to agree with.

    20 000 events, a contiguous 3 000-event contaminated block starting at 7 000.  With
    the default ``step`` there are not enough bins for the isolation tree, so PeacoQC
    removes on MAD alone -- and the test asserts that this is *reported* rather than
    presented as the ``determine_good_cells="all"`` that was asked for.
    """

    _need_peacoqc()
    settings = PeacoQCSettings(
        channels=("CD4", "CD3"), criteria="all", mad=7.0, it_limit=0.6, report=False
    )
    n, lo, hi = 20_000, 7_000, 10_000
    result = _call(tmp_path, settings, n, block=(lo, hi))
    assert result.status == "ok", f"{result.status}: {result.reason}"
    assert result.n_removed > 0, (
        "PeacoQC removed nothing from a fixture with a 15% contaminated block. A run that "
        "reports no removals here is indistinguishable from a broken adapter, so this fails "
        "rather than passing"
    )
    removed = result.removed
    assert removed.columns.tolist() == [
        "sample_id", "original_event_index", "step", "reason"
    ]
    assert set(removed["sample_id"]) == {"s1"}
    assert set(removed["step"]) == {"peacoqc"}
    assert removed["reason"].str.startswith(f"PeacoQC{VERIFIED_VERSION}:").all()
    assert "determine_good_cells=all" in removed["reason"].iloc[0]

    # The ground truth: every event in the block, identified by its ORIGINAL index.
    block = set(range(lo, hi))
    survivors = block - set(removed["original_event_index"])
    assert not survivors, (
        f"{len(survivors)} event(s) of the known contaminated block survived, first "
        f"{sorted(survivors)[:5]}; the fixture's ground truth and PeacoQC's answer "
        "disagree"
    )
    outcome = result.samples[0]
    assert outcome.n_input == n
    assert outcome.n_removed == len(removed)
    assert outcome.removed_fraction == pytest.approx(100.0 * len(removed) / n, abs=1e-6)

    # The isolation tree did NOT run, and the result says so.
    assert outcome.it_fraction is None
    assert outcome.it_available is False
    assert outcome.it_reason and "isolation tree" in outcome.it_reason
    assert outcome.criteria_applied == ("MAD",), (
        "with the isolation tree unavailable, determine_good_cells='all' was decided on "
        "MAD alone; reporting it as anything else misstates what was done"
    )
    assert any("MAD-only pass recorded as such" in n for n in outcome.notes), outcome.notes
    assert bool(result.per_sample()["it_available"].iloc[0]) is False


@requires_r
def test_peacoqc_runs_the_isolation_tree_when_there_are_enough_bins(tmp_path: Path) -> None:
    """The other regime, so the MAD-only path above is not the only one exercised.

    ``step = 250`` over 100 000 events gives 400 bins, enough for the isolation tree.  It
    runs and removes on IT; MAD removes nothing.  The criteria that actually contributed
    are named, so a reader can tell this run from the MAD-only one.
    """

    _need_peacoqc()
    settings = PeacoQCSettings(
        channels=("CD4", "CD3"), criteria="all", mad=7.0, it_limit=0.6,
        step=250, report=False,
    )
    n, lo, hi = 100_000, 35_000, 50_000
    result = _call(tmp_path, settings, n, block=(lo, hi))
    outcome = result.samples[0]
    assert outcome.it_available is True, (
        f"the isolation tree did not run here ({outcome.it_reason}); the test meant to "
        "exercise it would then be asserting the MAD path a second time"
    )
    assert outcome.it_fraction is not None and outcome.it_fraction > 0
    assert outcome.criteria_applied[0] == "IT"
    assert not outcome.notes
    assert set(range(lo, hi)) <= set(result.removed["original_event_index"])


@requires_r
def test_peacoqc_removal_indices_resolve_to_original_event_indices_not_row_positions(
    tmp_path: Path,
) -> None:
    """The identity check, with the event indices deliberately *not* 0..n-1.

    PeacoQC returns row positions in the file it read.  Those are event indices only if
    the file was written in event order starting at zero, and this build makes no such
    assumption: it writes a dedicated identity channel and resolves through it.  Here the
    indices are reversed and offset, so a row-position answer would be wrong for every
    row.
    """

    _need_peacoqc()
    n = 20_000
    values = _fixture(n, (7_000, 10_000))
    # Reversed, doubled and offset, so the event indices are *disjoint from the row
    # positions*: a row-position answer would be below 20 000 and even, this is above
    # 1 000 000. Any value under the floor is therefore decisive evidence that the
    # identity channel was not used.
    event_index = (np.arange(n)[::-1] * 2 + 1_000_000).astype(np.int64)
    result = run_peacoqc(
        values, CHANNELS,
        sample_ids=np.repeat("s1", n),
        original_event_index=event_index,
        settings=PeacoQCSettings(channels=("CD4", "CD3"), criteria="MAD", mad=7.0,
                                 report=False),
        work_dir=tmp_path / "work",
    )
    removed = set(result.removed["original_event_index"])
    assert removed, "nothing was removed from a fixture with a contaminated block"
    assert min(removed) >= 1_000_000, (
        f"a removal index came back at {min(removed)}, below the identity channel's "
        "offset floor, so it is a row position and the identity channel was not used"
    )
    # The contaminated block is a range of ROWS; the log must name the event indices those
    # rows carry. A row-position answer would have produced a contiguous 7000..9999.
    assert set(event_index[7_000:10_000]) <= removed, (
        "the known contaminated block was not removed by its original event index"
    )
    assert len(removed) == result.samples[0].n_removed
    # The file PeacoQC read is written in EVENT order, not in the caller's row order, so
    # its identity channel is the sorted event index. That permutation is the whole reason
    # a dedicated channel is needed: resolving a row position against the caller's own
    # ordering would have been wrong here.
    written = read_fcs(tmp_path / "work" / "fcs_in" / "s1.fcs")
    assert IDENTITY_CHANNEL in written.channel_names
    assert np.array_equal(
        np.rint(
            written.events[:, written.channel_names.index(IDENTITY_CHANNEL)]
        ).astype(np.int64),
        np.sort(event_index),
    ), "the identity channel did not survive the FCS round trip, or was written in the "\
       "caller's row order rather than in event order"


@requires_r
def test_two_samples_are_qced_separately_and_the_log_names_the_sample(tmp_path: Path) -> None:
    """Per-sample attribution, and no cross-contamination of one sample's removals."""

    _need_peacoqc()
    n = 6_000
    values = np.vstack([_fixture(n, (2_000, 2_600), seed=1), _fixture(n, (1_000, 1_200), seed=2)])
    result = run_peacoqc(
        values, CHANNELS,
        sample_ids=np.array(["s1"] * n + ["s2"] * n),
        original_event_index=np.concatenate([np.arange(n), np.arange(n)]),
        settings=PeacoQCSettings(channels=("CD4", "CD3"), criteria="MAD", mad=7.0,
                                 report=False),
        work_dir=tmp_path / "work",
    )
    assert {s.sample_id for s in result.samples} == {"s1", "s2"}
    for outcome in result.samples:
        mine = result.removed[result.removed["sample_id"] == outcome.sample_id]
        assert len(mine) == outcome.n_removed
    summary = result.per_sample()
    assert {"n_input", "n_removed", "criteria_applied", "it_available"} <= set(summary.columns)
    assert (summary["n_input"] == n).all()
    assert (tmp_path / "work" / "fcs_in" / "s1.fcs").exists()
    assert (tmp_path / "work" / "fcs_in" / "s2.fcs").exists()


@requires_r
def test_a_clean_sample_reports_zero_removals_and_the_log_is_empty_not_absent(
    tmp_path: Path,
) -> None:
    """The legitimate zero, and it looks different from the unavailable case.

    A sample with no contamination and a MAD of 6 removes nothing.  That is a real result
    and it is reported as ``n_removed = 0`` with the IT status stated -- which is exactly
    what an unrunnable call must never produce, and the reason the two are told apart by
    status and reason rather than by the removal count.
    """

    _need_peacoqc()
    rng = np.random.default_rng(7)
    n = 6_000
    clean = np.column_stack(
        [
            np.abs(rng.normal(150_000, 18_000, n)), np.abs(rng.normal(60_000, 11_000, n)),
            np.abs(rng.normal(140_000, 17_000, n)), np.abs(rng.normal(1_800, 700, n)),
            np.abs(rng.normal(2_500, 700, n)), np.abs(rng.normal(1_800, 600, n)),
        ]
    )
    result = run_peacoqc(
        clean, CHANNELS,
        sample_ids=np.repeat("s1", n),
        original_event_index=np.arange(n),
        settings=PeacoQCSettings(channels=("CD4", "CD3"), criteria="MAD", mad=6.0,
                                 report=False),
        work_dir=tmp_path / "work",
    )
    assert result.status == "ok"
    assert result.samples[0].n_input == n
    assert result.samples[0].n_removed == 0
    assert result.removed.empty, (
        "a zero-removal run must leave an empty log carrying the four declared columns, "
        "not rows for events that were not removed"
    )
    assert result.removed.columns.tolist() == [
        "sample_id", "original_event_index", "step", "reason"
    ]
    # And the IT status is stated even on a clean sample: a clean sample and a sample
    # that was never assessed must not look alike.
    assert result.samples[0].it_available is False
    assert result.samples[0].it_reason


@requires_r
def test_the_script_result_payload_names_the_versions_that_produced_it(tmp_path: Path) -> None:
    """A result file written without a PeacoQC call would report a QC pass that never
    happened, so the payload must carry the version this module's formals were read from."""

    _need_peacoqc()
    _call(tmp_path, PeacoQCSettings(channels=("CD4",), criteria="MAD", mad=6.0, report=False),
          n=2_000, block=(500, 600))
    payload = json.loads((tmp_path / "work" / "result.json").read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["peacoqc_version"] == VERIFIED_VERSION
    assert payload["channels"] == ["CD4"]
    assert payload["criteria"] == "MAD"
    assert payload["r_version"].startswith("4.")
    assert "s1" in payload["samples"], (
        "results are keyed by sample_id, not by file name: the Python side knows samples "
        "by the manifest's sample_id and the two must not disagree by an extension on the "
        "field the removal log is built from"
    )
    assert payload["samples"]["s1"]["n_input"] == 2_000
    assert "GoodCells" in payload["samples"]["s1"]["result_keys"], (
        "the payload must name the result element the removal index is derived from, so a "
        "package update that renamed it fails here instead of reporting no removals"
    )
    assert len(payload["samples"]["s1"]["good_cells"]) == 2_000


@requires_r
def test_the_report_setting_is_honoured_and_reported_as_counted(tmp_path: Path) -> None:
    """``report`` is a declared setting, and the artifact is counted, not promised.

    Measured against PeacoQC 1.16.0 rather than read off its documentation: a report file
    is written on both the ``all`` and the ``MAD`` route, so the earlier assumption that
    ``MAD`` produces none is wrong and is not asserted.  What the test does check is that
    the setting changes what is on disk, and that the count reported by the script
    matches the files that are there -- so a settings file asking for a report can never
    be silently unfulfilled.
    """

    _need_peacoqc()
    with_report = _call(
        tmp_path, PeacoQCSettings(channels=("CD4", "CD3"), criteria="all", mad=7.0,
                                  it_limit=0.6, report=True),
        n=6_000, block=(2_000, 2_600),
    )
    directory = Path(with_report.report_paths["s1"])
    assert directory.is_dir(), f"PeacoQC wrote no output directory at {directory}"
    files = with_report.output_files["s1"]
    assert files, "report = TRUE produced no file; the setting was accepted and did nothing"
    assert all(Path(f).exists() for f in files), files
    assert any(Path(f).name.startswith("PeacoQC_report") for f in files), files

    without = tmp_path / "off"
    no_report = _call(
        without, PeacoQCSettings(channels=("CD4", "CD3"), criteria="MAD", mad=7.0, report=False),
        n=6_000, block=(2_000, 2_600),
    )
    assert not no_report.output_files["s1"], (
        "report = FALSE still produced output files; the setting is not being read"
    )


# --------------------------------------------------------------------------------------
# Item 3.6, self-checked
# --------------------------------------------------------------------------------------


def test_every_r_function_this_build_calls_is_declared_in_the_upstream_lock(
    repo_root: Path,
) -> None:
    """The Stage 0 gate would not catch this, so this test does.

    ``PeacoQC::PeacoQC`` is the upstream R function this build newly calls.
    ``test_every_adapter_referenced_r_function_appears_in_the_lock`` only scans the
    cyCombine and CytoBatchFlagR sources, and the ``spiscy`` lock entry has no
    ``adapter_calls`` list at all, so nothing in the suite would notice a missing
    declaration.  This test fails, specifically, until the lock entry gains one -- which
    is the intended state: a missing declaration is a real defect and is reported as one
    rather than skipped.
    """

    script = (repo_root / "workflow" / "scripts" / "peacoqc_qc.R").read_text(encoding="utf-8")
    # R namespaces are separated by "::" and packages by "."; the lock records bare
    # function names (as the cycombine and cytobatchflagr entries do), so both
    # separators have to be stripped.  Splitting on "." alone left
    # "PeacoQC::PeacoQC" intact, so this check could never match its own declaration
    # and would have failed for a correctly declared lock.
    used = {name.split("::")[-1].split(".")[-1] for name in ("PeacoQC::PeacoQC",)}
    assert "PeacoQC::PeacoQC" in script
    lock = yaml.safe_load((repo_root / "config" / "upstream_lock.yaml").read_text(encoding="utf-8"))
    declared = set(lock["resources"]["spiscy"].get("adapter_calls") or [])
    assert not (used - declared), (
        "config/upstream_lock.yaml: the spiscy entry declares "
        f"{sorted(declared) or 'no adapter_calls'}, so the R function this build now "
        f"calls ({sorted(used)}) is recorded nowhere the Stage 0 audit reads. "
        "Expected YAML block, to be added under the spiscy entry:\n\n"
        "    adapter_calls:\n"
        "      - PeacoQC\n"
        "      - flowVS.transFlowVS\n"
    )


# --------------------------------------------------------------------------------------
# The reference implementation's own configuration
# --------------------------------------------------------------------------------------


def test_the_reference_qc_settings_file_loads_here_unchanged(repo_root: Path) -> None:
    """``06_QC.yaml`` is read, not reinterpreted.

    The reference implementation names its thresholds ``MAD`` and ``IT`` and its channel
    list ``markers_to_evaluate``, none of which is a formal of ``PeacoQC``.  They are
    accepted as aliases so that file loads as written; if it stops loading, the divergence
    is a parse error rather than a quiet reinterpretation.
    """

    path = repo_root / ".upstream" / "spiscy" / "config" / "06_QC.yaml"
    if not path.exists():
        pytest.skip("the pinned SPiSCy checkout is absent; run scripts/fetch_upstream.sh")
    settings = read_peacoqc_settings(path)
    assert settings.criteria == "all", "the reference implementation's own default criteria"
    assert settings.mad == 7.0 and settings.it_limit == 0.6
    assert "CD4" in settings.channels and "ST2" in settings.channels
    assert len(settings.channels) == 15


def test_the_reference_gate_table_loads_here_unchanged(repo_root: Path) -> None:
    """``02_prelim_gating.yaml`` is read by :func:`read_gates` as written."""

    path = repo_root / ".upstream" / "spiscy" / "config" / "02_prelim_gating.yaml"
    if not path.exists():
        pytest.skip("the pinned SPiSCy checkout is absent; run scripts/fetch_upstream.sh")
    table = read_gates(path)
    assert [g.name for g in table.gates] == ["debris", "singlets", "live"]
    assert table.channels == ("FSC-A", "SSC-A", "FSC-H", "FVS780")
    singlets = next(g for g in table.gates if g.name == "singlets")
    assert singlets.kind == "polygon" and singlets.vertices is not None
    assert singlets.vertices.shape == (5, 2), "the file closes the polygon explicitly"
    assert int(table.to_frame()["reason"].str.contains("polygon").sum()) == 1
