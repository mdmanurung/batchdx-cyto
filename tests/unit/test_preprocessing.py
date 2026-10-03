"""Preprocessing: the gate geometry, the transform convention, and the removal log.

Three groups, and the grouping is the point:

* **No R, no data.**  Gate reading and evaluation, including the point-in-polygon test
  against a hand-computed answer and against the reference implementation's own shipped
  gate table; the arcsinh convention, cross-checked against the real
  ``io.export.InversionPlan.apply`` rather than against a copy of it; the removal log's
  writer, reader, dtypes and refusals.
* **The reference implementation's shapes.**  The pinned SPiSCy
  ``config/02_prelim_gating.yaml`` and ``config/06_QC.yaml`` are read here, so a
  divergence between this build and the reference shows up as a parse failure instead of
  as a silent reinterpretation.
* **The unavailable path.**  What the configuration says when a step cannot run, and
  what the schema refuses.

Every test in this file runs with no R interpreter and no network.  The tests that need
an installed PeacoQC are in ``tests/integration/test_peacoqc_backend.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from batchdx.cyto.io.export import INVERTIBLE, plan_inversion
from batchdx.cyto.io.manifests import ManifestError, read_removed_events, write_removed_events
from batchdx.cyto.preprocess import (
    ASINH_FORM,
    ASINH_INVERSE_FORM,
    SCALED_ASINH_FORM,
    STEP_GATE,
    STEP_QC,
    GatingError,
    PeacoQCError,
    PeacoQCSettings,
    RemovalLogError,
    TransformError,
    agree_with_exporter_inverse,
    apply_gates,
    asinh_invert,
    asinh_round_trip_error,
    asinh_transform,
    check_removal_log,
    empty_removal_log,
    gates_from_config,
    points_in_polygon,
    read_gates,
    read_peacoqc_settings,
    removal_frame,
    resolve_convention,
    scaled_asinh_mismatch,
    settings_from_config,
    summarise_removals,
)
from batchdx.cyto.schemas.contracts import ProcessingManifest

#: The reference implementation's three gates, transcribed from
#: ``.upstream/spiscy/config/02_prelim_gating.yaml``.  The integration test reads the
#: file itself; this copy is here so the unit suite is independent of the checkout.
SPISCY_GATES = """
gate_chain: [debris, singlets, live]
gates:
  debris:
    filterId: "debris"
    mk1: "FSC-A"
    mk2: "SSC-A"
    type: "rectangle"
    coords:
      mk1: [100000, 250000]
      mk2: [25000, 150000]
  singlets:
    filterId: "singlets"
    mk1: "FSC-A"
    mk2: "FSC-H"
    type: "polygon"
    coords:
      - [50000, 25000]
      - [400000, 230000]
      - [350000, 275000]
      - [50000, 110000]
      - [50000, 25000]
  live:
    filterId: "live"
    mk1: "FVS780"
    type: "rectangle"
    coords:
      mk1: [-2000, 4000]
"""


# --------------------------------------------------------------------------------------
# 3.1  The gating model
# --------------------------------------------------------------------------------------


def test_the_reference_implementation_ships_a_polygon_and_this_build_evaluates_it() -> None:
    """The reason the gate model is two-dimensional, asserted rather than described.

    Two of the three gates the reference implementation actually configures are
    rectangles; the third, ``singlets``, is a **polygon** on FSC-A against FSC-H.  A
    one-dimensional threshold table cannot express it, and this test fails if someone
    reduces the gate model to one dimension and calls the coverage equivalent.
    """

    table = read_gates_from_text(SPISCY_GATES)
    kinds = {g.name: g.kind for g in table.gates}
    assert kinds == {"debris": "rectangle", "singlets": "polygon", "live": "rectangle"}
    assert table.channels == ("FSC-A", "SSC-A", "FSC-H", "FVS780"), (
        "the gate channels are scatter and viability parameters, not panel markers; the "
        "order is the chain's, and it is what the application walks"
    )


def test_a_polygon_is_evaluated_against_a_hand_computed_inside_and_outside_set(
    tmp_path: Path,
) -> None:
    """The crossing-number test, against points whose answer is known by hand.

    The square is (0,0), (4,0), (4,4), (0,4).  Inside: (2,2).  Outside on each of the
    four sides.  A gate that got the boundary convention wrong, or transposed the axes,
    or tested ``>=``/``<=`` the wrong way round, would move at least one of these.
    """

    square = np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]])
    x = np.array([2.0, -1.0, 5.0, 2.0, 2.0, 0.5, 3.5])
    y = np.array([2.0, 2.0, 2.0, -1.0, 5.0, 0.5, 3.5])
    inside = points_in_polygon(x, y, square)
    assert inside.tolist() == [True, False, False, False, False, True, True]


def test_a_non_convex_polygon_excludes_the_notch() -> None:
    """A concave polygon is where a bounding-box shortcut would give the wrong answer.

    An L-shape: the notch at (2,2)-(3,3) is outside the polygon and inside its bounding
    box.  A gate implemented as "inside the box" would admit those points.
    """

    ell = np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [1.0, 1.0], [1.0, 4.0], [0.0, 4.0]])
    inside = points_in_polygon(
        np.array([0.5, 2.5, 0.5, 3.5]), np.array([3.5, 2.5, 0.5, 0.5]), ell
    )
    assert inside.tolist() == [True, False, True, True], (
        "the point in the L's notch is inside the bounding box and outside the polygon; "
        "a bounding-box gate would have admitted it"
    )


def test_the_two_dimensional_gates_change_the_outcome_and_the_one_dimensional_one_would_not(
    tmp_path: Path,
) -> None:
    """The finding behind the decision, measured.

    A ``singlets`` gate written as two independent one-dimensional thresholds is a box.
    This builds the diagonal polygon and the box that a 1-D table would have to become,
    and shows they select different events.  That is the whole argument for a 2-D gate
    table: the weaker gate is not a simplification of the stronger one, it is a different
    gate that lets doublet-shoulder events through.
    """

    table = read_gates_from_text(
        """
        gate_chain: [singlets]
        gates:
          singlets:
            filterId: singlets
            mk1: FSC-A
            mk2: FSC-H
            type: polygon
            coords:
              - [0, 0]
              - [10, 10]
              - [10, 4]
              - [4, 10]
              - [0, 0]
        """
    )
    gate = table.gates[0]
    # Two events with the *same* FSC-A, one on the diagonal and one off it.  FSC-A alone
    # cannot separate them; the polygon can.
    x = np.array([8.0, 8.0])
    y = np.array([8.0, 1.0])
    assert gate.contains(x, y).tolist() == [True, False], (
        "the diagonal polygon keeps the on-diagonal event and rejects the shoulder"
    )
    # The shape a one-dimensional table would have to take: each axis thresholded
    # independently over its own range, i.e. the bounding box.  For a diagonal band that
    # box is not the band -- here it does not even keep the event the polygon keeps.
    box = (x >= 0) & (x <= 10) & (y >= 0) & (y <= 4)
    assert box.tolist() == [False, True], (
        "the axis-aligned box disagrees with the diagonal polygon on BOTH events. One-"
        "dimensional gating is not a weaker version of two-dimensional gating; it is a "
        "different gate, and here it is not even a consistent one"
    )
    assert x[0] == x[1], (
        "the two events share FSC-A exactly, so no threshold on FSC-A -- at any value -- "
        "can separate them. That is the whole argument for a 2-D gate table"
    )


def test_a_gate_type_this_build_does_not_implement_is_refused_by_name() -> None:
    """A gate that was never evaluated must not be reported as a gate that passed."""

    with pytest.raises(GatingError) as excinfo:
        read_gates_from_text(
            """
            gate_chain: [el]
            gates:
              el:
                filterId: el
                mk1: FSC-A
                mk2: SSC-A
                type: ellipsoid
                coords: {}
            """
        )
    message = str(excinfo.value)
    assert "ellipsoid" in message
    assert "not implemented" in message
    assert "never evaluated" in message


def test_a_gate_defined_but_not_in_the_chain_is_refused() -> None:
    """A declared gate that is never applied is the accepted-and-ignored outcome."""

    with pytest.raises(GatingError, match="not in gate_chain"):
        read_gates_from_text(
            """
            gate_chain: [debris]
            gates:
              debris:
                filterId: debris
                mk1: FSC-A
                type: rectangle
                coords: {mk1: [0, 1]}
              singlets:
                filterId: singlets
                mk1: FSC-A
                type: rectangle
                coords: {mk1: [0, 2]}
            """
        )


def test_a_gate_channel_the_data_lacks_is_refused_rather_than_gated_on_zero() -> None:
    """A missing channel is not a constant zero, and a gate on a constant removes all."""

    table = read_gates_from_text(SPISCY_GATES)
    with pytest.raises(GatingError) as excinfo:
        table.validate(["FSC-A", "SSC-A"])  # no FSC-H, no FVS780
    message = str(excinfo.value)
    assert "FSC-H" in message and "FVS780" in message
    assert "not zero" in message


def test_gates_are_chained_so_an_event_is_counted_once_at_the_gate_that_removed_it(
    tmp_path: Path,
) -> None:
    """Two gates, and the accounting: no double counting, and the reason names the gate."""

    table = read_gates_from_text(
        """
        gate_chain: [debris, singlets]
        gates:
          debris:
            filterId: debris
            mk1: FSC-A
            type: rectangle
            coords: {mk1: [0, 10]}
          singlets:
            filterId: singlets
            mk1: FSC-A
            mk2: FSC-H
            type: rectangle
            coords: {mk1: [0, 10], mk2: [0, 10]}
        """
    )
    # 8 events.  Four fail the debris gate (FSC-A above 10); of the four that survive,
    # one fails the singlets gate (FSC-H above 10).
    fsc_a = np.array([1.0, 2.0, 20.0, 30.0, 40.0, 5.0, 6.0, 50.0])
    fsc_h = np.array([1.0, 2.0, 20.0, 30.0, 40.0, 1.0, 50.0, 1.0])
    values = np.column_stack([fsc_a, fsc_h])
    result = apply_gates(
        table, values, ["FSC-A", "FSC-H"],
        sample_ids=np.array(["s1"] * 8),
        original_event_index=np.arange(8),
    )
    assert result.keep.tolist() == [True, True, False, False, False, True, False, False]
    assert [o.removed for o in result.outcomes] == [4, 1]
    assert [o.n_input for o in result.outcomes] == [8, 4], (
        "the second gate is evaluated on the survivors, so the removed-by-gate counts sum "
        "to the number of events actually lost and no event is charged to two gates"
    )
    assert list(result.removed["original_event_index"]) == [2, 3, 4, 7, 6]
    assert set(result.removed["step"]) == {STEP_GATE}
    assert result.removed["reason"].str.startswith("debris:rectangle").sum() == 4
    assert result.removed["reason"].str.startswith("singlets:rectangle").sum() == 1
    assert result.outcomes[1].removed_by_sample == {"s1": 1}


def test_a_one_axis_rectangle_is_a_half_line_and_is_accepted() -> None:
    """The reference implementation's own ``live`` gate is one axis; it must load."""

    table = read_gates_from_text(SPISCY_GATES)
    live = next(g for g in table.gates if g.name == "live")
    assert live.channels == ("FVS780",)
    # A negative lower bound is legal: viability can be negative before compensation.
    assert live.contains(np.array([-2000.0, 0.0, 4000.0, 4001.0])).tolist() == [
        True, True, True, False
    ]


def test_a_polygon_on_one_axis_is_refused() -> None:
    with pytest.raises(GatingError, match="needs two channels"):
        read_gates_from_text(
            """
            gate_chain: [p]
            gates:
              p:
                filterId: p
                mk1: FSC-A
                type: polygon
                coords: [[0, 0], [1, 0], [0, 1]]
            """
        )


def test_gates_from_config_reads_the_declared_key_and_nothing_is_inferred(
    repo_root: Path, tmp_path: Path
) -> None:
    """``preprocessing.gating`` is a path, and it is read."""

    from batchdx.cyto.config import parse_config

    (tmp_path / "gates.yaml").write_text(SPISCY_GATES, encoding="utf-8")
    base = yaml.safe_load(
        (repo_root / "config" / "example_paired_reference.yaml").read_text(encoding="utf-8")
    )
    base["preprocessing"] = dict(base.get("preprocessing") or {}, gating="gates.yaml")
    config = parse_config(base, source="test.yaml")
    table = gates_from_config(config, tmp_path)
    assert table is not None
    assert [g.name for g in table.gates] == ["debris", "singlets", "live"]

    # And with no key declared it declines, rather than inventing a gate table.
    base["preprocessing"] = {"enabled": False}
    assert gates_from_config(parse_config(base, source="test.yaml"), tmp_path) is None


# --------------------------------------------------------------------------------------
# 3.2  The transform convention
# --------------------------------------------------------------------------------------


def test_the_transform_is_the_unscaled_asinh_that_flowvs_computes() -> None:
    """The forward form, as a number, for the values the gates above would pass on.

    ``flowVS::transFlowVS``'s body is one line, ``asinh(mat[, channels[j]] /
    cofactors[j])``, and its own ``Details`` section says the same in words.  These are
    the numbers it returns for a cofactor of 5, computed here independently.  The
    integration test ``test_flowvs_transflowvs_agrees_with_this_build_elementwise``
    executes the R function and compares; this one pins the arithmetic so the test fails
    on *which* value is wrong rather than only that the two disagree.
    """

    raw = np.array([0.0, 100.0, 2500.0, 50000.0, 250000.0, 1e6])
    got = asinh_transform(raw, 5.0)
    assert got == pytest.approx(
        [0.0, 3.689503869, 6.907756279, 9.903487555, 11.512925465, 12.899219826], abs=1e-9
    ), (
        "flowVS divides by the cofactor and does NOT multiply by it on the way out. If "
        "this assertion fails, this build has adopted the scaled convention and "
        "io.export's inverse is no longer exact"
    )


def test_the_exporter_inverse_is_exact_for_this_convention_and_is_the_same_function() -> None:
    """Item 3.2, executed against the real exporter rather than a copy of it.

    ``io.export.InversionPlan.apply`` computes ``sinh(x) * cofactor``.  That is the
    identity map's exact inverse of ``asinh(raw / cofactor)``, so on this build's
    convention the exporter's exact-inverse claim is **true**.  The cross-check calls the
    exporter's own method; a test comparing two literals would only prove this module is
    internally consistent.
    """

    raw = np.array([0.0, 1.0, 100.0, 2500.0, 50000.0, 250000.0, 1e6])
    for cofactor in (5.0, 1.0, 150.0, 0.5):
        absolute, relative = asinh_round_trip_error(raw, cofactor)
        assert absolute < 1e-6, (cofactor, absolute)
        assert relative < 1e-9, (cofactor, relative)
        assert agree_with_exporter_inverse(raw, cofactor), (
            f"cofactor {cofactor}: batchdx.cyto.preprocess.asinh_invert and "
            "io.export.InversionPlan.apply have diverged"
        )


def test_the_scaled_convention_is_not_what_this_build_inverts_and_breaks_badly() -> None:
    """The other reading of the exporter's formula, and what it would cost.

    If the baseline were on ``asinh(raw / cofactor) * cofactor`` -- the reading the
    work-package plan recorded -- then ``sinh(x) * cofactor`` would not recover ``raw``
    at all, and the exporter's note would be false.  It is not on that convention, so the
    note is true.  This test puts a number on the difference so the distinction cannot be
    quietly forgotten.
    """

    raw = np.array([0.0, 100.0, 2500.0, 50000.0, 250000.0, 1e6])
    mismatch = scaled_asinh_mismatch(raw, 5.0)
    assert mismatch > 1e10, (
        f"the scaled convention is only {mismatch:g} away from being inverted exactly, "
        "which would mean the two conventions had converged and the distinction in "
        "preprocess/transform.py was describing a distinction that no longer exists"
    )
    # And the scaled convention is a real, different function of the input.
    scaled = asinh_transform(raw, 5.0) * 5.0
    assert not np.allclose(scaled, asinh_transform(raw, 5.0))
    assert ASINH_FORM == "asinh(raw / cofactor)"
    assert ASINH_INVERSE_FORM == "sinh(x) * cofactor"
    assert SCALED_ASINH_FORM == "asinh(raw / cofactor) * cofactor"


def test_the_exporter_declares_asinh_the_only_invertible_transformation() -> None:
    """Cross-checked from the other side: the exporter's own table and its own plan."""

    assert INVERTIBLE == {"asinh"}
    plan = plan_inversion(
        ProcessingManifest(
            scale="asinh_cofactor5", transformation="asinh", transformation_parameters={"cofactor": 5.0}
        )
    )
    assert plan.invertible and plan.cofactor == 5.0
    values = asinh_transform(np.array([2500.0, 1e6]), 5.0)
    assert plan.apply(values) == pytest.approx([2500.0, 1e6], rel=1e-9), (
        "plan_inversion claims an exact inverse; on this build's convention it is one"
    )


def test_a_non_positive_or_non_finite_cofactor_is_refused() -> None:
    """There is no universal cofactor, and a zero one has no asinh interpretation."""

    for bad in (0.0, -5.0, float("nan"), float("inf")):
        with pytest.raises(TransformError, match="strictly positive"):
            asinh_transform(np.array([1.0]), bad)
        with pytest.raises(TransformError, match="strictly positive"):
            asinh_invert(np.array([1.0]), bad)


def test_logicle_is_declared_unimplemented_and_says_why() -> None:
    """Item 3.2's logicle question, answered in the code and not only in the report."""

    convention = resolve_convention("logicle")
    assert not convention.invertible
    assert "no closed-form inverse" in convention.note
    assert "no logicle" in convention.note
    assert convention.forward != ASINH_FORM
    with pytest.raises(TransformError, match="no transformation was declared"):
        resolve_convention("none")
    with pytest.raises(TransformError, match="is not implemented"):
        resolve_convention("biexponential")


# --------------------------------------------------------------------------------------
# 3.4  The removal log
# --------------------------------------------------------------------------------------


def test_the_removal_log_declares_exactly_the_four_columns_dataset_removed_events_names() -> None:
    """The writer and the field agree, so a log written here populates the field.

    ``Dataset.removed_events`` is declared in ``io/dataset.py:54-56`` with these four
    columns.  The test reads the declaration out of the dataclass rather than from this
    module, so a change to the field without a change to the writer fails here.
    """

    from batchdx.cyto.io.dataset import Dataset

    declared = list(Dataset.__dataclass_fields__["removed_events"].default_factory().columns)
    assert declared == ["sample_id", "original_event_index", "step", "reason"]
    assert empty_removal_log().columns.tolist() == declared


def test_a_removal_row_carries_typed_indices_and_a_reason_naming_the_criterion() -> None:
    frame = removal_frame(
        sample_id=["s1", "s1", "s2"],
        original_event_index=np.array([4.0, 9.0, 0.0]),  # as a CSV would read it back
        step=STEP_GATE,
        reason="debris:rectangle(FSC-A/SSC-A=[1e+05,2.5e+05],[25000,150000])",
    )
    assert frame["original_event_index"].dtype == np.int64
    assert pd.api.types.is_string_dtype(frame["sample_id"])
    assert frame["original_event_index"].tolist() == [4, 9, 0]
    assert frame["sample_id"].tolist() == ["s1", "s1", "s2"]


def test_a_removal_without_a_reason_or_with_an_unknown_step_is_refused() -> None:
    """An uninterpretable log is worse than no log."""

    with pytest.raises(RemovalLogError, match="not a record"):
        removal_frame(sample_id=["s1"], original_event_index=[0], step=STEP_GATE, reason="  ")
    with pytest.raises(RemovalLogError, match="not one this build names"):
        removal_frame(
            sample_id=["s1"], original_event_index=[0], step="clipping", reason="whatever"
        )


def test_a_log_that_grew_columns_is_refused() -> None:
    """Columns nothing renders are how a log starts accumulating silently."""

    frame = removal_frame(
        sample_id=["s1"], original_event_index=[0], step=STEP_QC, reason="PeacoQC"
    )
    frame["confidence"] = 0.9
    with pytest.raises(RemovalLogError, match="unexpected column"):
        check_removal_log(frame)


def test_a_zero_row_log_file_is_refused_rather_than_read_as_nothing_was_removed() -> None:
    """A missing value is never a zero.  An empty file is a claim, and a wrong one."""

    path = Path("/tmp") / "never-written-by-this-test.tsv"
    path.write_text("sample_id\toriginal_event_index\tstep\treason\n", encoding="utf-8")
    try:
        with pytest.raises(ManifestError, match="no rows"):
            read_removed_events(path)
    finally:
        path.unlink(missing_ok=True)
    with pytest.raises(ManifestError, match="not a log of no removals"):
        read_removed_events(path)


def test_the_removal_log_round_trips_through_a_file_with_its_dtypes_intact(
    tmp_path: Path,
) -> None:
    frame = pd.concat(
        [
            removal_frame(
                sample_id=["s1", "s1"], original_event_index=[4, 9], step=STEP_GATE,
                reason="debris:rectangle(FSC-A/SSC-A)",
            ),
            removal_frame(
                sample_id=["s2"], original_event_index=[0], step=STEP_QC,
                reason="PeacoQC1.16.0:determine_good_cells=all,MAD=6,IT_limit=0.6",
            ),
        ],
        ignore_index=True,
    )
    path = write_removed_events(frame, tmp_path / "removed.tsv")
    back = read_removed_events(path)
    pd.testing.assert_frame_equal(back, frame, check_dtype=False)
    assert back["original_event_index"].dtype == np.int64
    assert back["sample_id"].tolist() == ["s1", "s1", "s2"]
    summary = summarise_removals(back)
    assert summary["n_removed"].sum() == 3
    assert set(summary["step"]) == {STEP_GATE, STEP_QC}
    # And CSV works too, so the separator is declared rather than guessed.
    csv = write_removed_events(frame, tmp_path / "removed.csv")
    pd.testing.assert_frame_equal(read_removed_events(csv), frame, check_dtype=False)


def test_a_removal_log_suffix_that_declares_no_separator_is_refused(tmp_path: Path) -> None:
    frame = removal_frame(
        sample_id=["s1"], original_event_index=[0], step=STEP_QC, reason="PeacoQC"
    )
    with pytest.raises(ManifestError, match="separator is declared"):
        write_removed_events(frame, tmp_path / "removed.log")


def test_summarising_an_empty_log_returns_no_rows_and_no_zero() -> None:
    summary = summarise_removals(empty_removal_log())
    assert summary.empty
    assert "n_removed" in summary.columns
    assert len(summary) == 0, "an empty log must not manufacture a row of zeroes"


# --------------------------------------------------------------------------------------
# 3.3 / 3.5  What the configuration says when a step cannot run
# --------------------------------------------------------------------------------------


def test_peacoqc_settings_are_read_from_the_file_and_the_channels_are_required(
    tmp_path: Path,
) -> None:
    """The channels are never inferred: PeacoQC assesses quality per channel."""

    from batchdx.cyto.config import parse_config

    path = tmp_path / "qc.yaml"
    # The reference implementation's own spellings: 06_QC.yaml says MAD, IT and
    # markers_to_evaluate, and that file must load here unchanged.
    path.write_text(
        "markers_to_evaluate: [CD4, CD3]\ncriteria: MAD\nMAD: 7\nIT: 0.6\n", encoding="utf-8"
    )
    settings = read_peacoqc_settings(path)
    assert settings.channels == ("CD4", "CD3")
    assert settings.criteria == "MAD"
    assert settings.mad == 7.0
    assert settings.it_limit == 0.6
    # criteria=MAD must not claim an IT threshold it did not apply.
    assert "IT_limit" not in settings.reason
    assert "MAD=7" in settings.reason

    (tmp_path / "bad.yaml").write_text("criteria: MAD\n", encoding="utf-8")
    with pytest.raises(PeacoQCError, match="never inferred"):
        read_peacoqc_settings(tmp_path / "bad.yaml")

    (tmp_path / "unknown.yaml").write_text("channels: [A]\nludicrous: 3\n", encoding="utf-8")
    with pytest.raises(PeacoQCError, match="accepted and then ignored"):
        read_peacoqc_settings(tmp_path / "unknown.yaml")

    # Two spellings of one setting is a conflict, not a convenience.
    (tmp_path / "twice.yaml").write_text("channels: [A]\nmad: 3\nMAD: 7\n", encoding="utf-8")
    with pytest.raises(PeacoQCError, match="second spelling"):
        read_peacoqc_settings(tmp_path / "twice.yaml")

    base = yaml.safe_load(
        (Path(__file__).resolve().parents[2] / "config" / "example_paired_reference.yaml").read_text()
    )
    base["preprocessing"] = dict(base.get("preprocessing") or {}, peacoqc="qc.yaml")
    config = parse_config(base, source="test.yaml")
    assert settings_from_config(config, tmp_path).channels == ("CD4", "CD3")
    base["preprocessing"] = {"enabled": False}
    assert settings_from_config(parse_config(base, source="test.yaml"), tmp_path) is None


def test_a_peacoqc_call_with_no_criteria_this_build_knows_is_refused() -> None:
    with pytest.raises(PeacoQCError, match="determine_good_cells"):
        PeacoQCSettings(channels=("CD4",), criteria="whatever")
    with pytest.raises(PeacoQCError, match="at least one channel"):
        PeacoQCSettings(channels=())


def test_a_blank_declared_path_is_refused_rather_than_resolved_to_the_cwd(
    repo_root: Path,
) -> None:
    """``gating: ""`` is a path to the working directory, not an absence."""

    from batchdx.cyto.config import ConfigError, parse_config

    base = yaml.safe_load(
        (repo_root / "config" / "example_paired_reference.yaml").read_text(encoding="utf-8")
    )
    for key in ("gating", "peacoqc"):
        raw = {**base, "preprocessing": {**(base.get("preprocessing") or {}), key: "  "}}
        with pytest.raises(ConfigError) as excinfo:
            parse_config(raw, source="test.yaml")
        assert f"preprocessing.{key}" in str(excinfo.value), str(excinfo.value)
        assert "leave the key out" in str(excinfo.value)


def test_preprocessing_time_channel_is_still_refused_and_says_why(repo_root: Path) -> None:
    """The one key of the three that stays unsupported, and it says why.

    ``preprocessing.gating`` and ``preprocessing.peacoqc`` left
    ``config.UNSUPPORTED_KEYS`` on the same change that gave them readers.
    ``preprocessing.time_channel`` did not, because nothing in Tranche 3 makes the
    acquisition-drift detector read a configured channel, and resurrecting it would
    recreate the exact failure the table exists to prevent.
    """

    from batchdx.cyto.config import UNSUPPORTED_KEYS, ConfigError, parse_config

    assert set(UNSUPPORTED_KEYS) == {"preprocessing.time_channel"}
    base = yaml.safe_load(
        (repo_root / "config" / "example_paired_reference.yaml").read_text(encoding="utf-8")
    )
    raw = {
        **base,
        "preprocessing": {**(base.get("preprocessing") or {}), "time_channel": "Time"},
    }
    with pytest.raises(ConfigError) as excinfo:
        parse_config(raw, source="test.yaml")
    message = str(excinfo.value)
    assert "preprocessing.time_channel" in message
    assert "acquisition order recorded at ingestion" in message


def test_the_run_still_never_transforms_events_even_with_the_new_keys(
    repo_root: Path,
) -> None:
    """The contract the new keys must not have broken.

    ``_validate_processing`` refuses a configuration that asks the run to transform a
    baseline the processing manifest says is already transformed.  Declaring a gate
    table and a PeacoQC pass changes nothing about that: they describe the upstream step
    that produced the baseline, not an instruction to the run.
    """

    from batchdx.cyto.config import parse_config
    from batchdx.cyto.io.dataset import DatasetValidationError, _validate_processing

    base = yaml.safe_load(
        (repo_root / "config" / "example_paired_reference.yaml").read_text(encoding="utf-8")
    )
    base["preprocessing"] = {
        "enabled": True,
        "transformation": "asinh",
        "cofactor_by_panel": {"panel_A": 5.0},
        "gating": "gates.yaml",
        "peacoqc": "qc.yaml",
    }
    config = parse_config(base, source="test.yaml")
    with pytest.raises(DatasetValidationError, match="refusing to transform twice"):
        _validate_processing(
            ProcessingManifest(scale="asinh_cofactor5", transformation="asinh"), config
        )


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------


def read_gates_from_text(text: str):
    import tempfile

    handle = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False, encoding="utf-8")
    try:
        handle.write(text)
        handle.close()
        return read_gates(handle.name)
    finally:
        Path(handle.name).unlink(missing_ok=True)
