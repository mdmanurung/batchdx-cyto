"""CytoNorm backend: the FCS round trip, the split, and the event-identity check.

Three groups of tests, in increasing order of what they need from the environment.

**No CytoNorm required.**  The capability declarations, the parameter refusals, the
FCS round trip this backend is built on, the reason a run reports
``not_available``, the reasons a run reports itself *ineligible* rather than
producing a correction, and — most importantly — the check that a corrected file
really is the baseline event for event.  That last one is tested by feeding the
verifier a deliberately reordered file, so the guarantee does not depend on CytoNorm
behaving.

**CytoNorm required** (``requires_cytonorm``).  One real correction: fit on the
training split's measured reference aliquots, apply to every sample, read the
written files back, and confirm that the technical offset between batches shrank.

The environment-dependent test skips *visibly* and says why.  It is not marked
xfail and the skip is not hidden: a skipped test is a recorded fact, whereas an
xfail would let a broken implementation look like an expected outcome.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from batchdx.cyto.adapters import (
    BackendContext,
    BackendError,
    BackendNotAvailable,
    available_backends,
    get_backend,
    resolve_parameters,
)
from batchdx.cyto.adapters.cytonorm import (
    FIT_APPLY_VERIFIED,
    IDENTITY_CHANNEL,
    TRANSFORM_LIST,
    EventIdentityError,
    verify_event_identity,
)
from batchdx.cyto.adapters.r_runner import RScriptFailed, r_package_version, run_r_script
from batchdx.cyto.config import load_config
from batchdx.cyto.design import audit_design
from batchdx.cyto.io import load_dataset
from batchdx.cyto.io.fcs import (
    FCSError,
    assert_fcs_round_trip,
    read_fcs,
    round_trip_tolerance,
    write_fcs,
)
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.sampling.splits import (
    SPLIT_NAMES,
    missing_arms,
    split_members,
    split_unit_field,
)
from batchdx.cyto.schemas.contracts import SamplingManifest
from batchdx.cyto.schemas.enums import CorrectionStatus, ReferenceRole
from batchdx.cyto.simulate import SampleSpec, SimulationSpec, write_simulation

#: Small on purpose.  The FlowSOM quantile fit needs enough cells per cluster per
#: batch; a few hundred cells per sample with nQ=99 would put most clusters on the
#: identity function, which is a legitimate upstream outcome but a useless smoke
#: test.  These values are stated rather than tuned silently.
SMOKE_N_CELLS = 600
SMOKE_PARAMS = {"nCells": 10_000, "xdim": 6, "ydim": 6, "nClus": 3, "nQ": 20}


def cytonorm_available() -> bool:
    return r_package_version("CytoNorm").available


def _requires_cytonorm() -> None:
    if not cytonorm_available():
        pytest.skip(
            "CytoNorm is not installed in this environment (see docs/blocked.md B-005). "
            "The Stage 6 gate for the cytonorm backend is recorded as not_available, not "
            "passed, until this test executes."
        )


# --------------------------------------------------------------------------------------
# The fixture: enough measured references in enough batches to fit on
# --------------------------------------------------------------------------------------


def _reference_rich_spec(n_cells: int = SMOKE_N_CELLS) -> SimulationSpec:
    """Three batches, two measured reference *donors* each, plus balanced biology.

    Two reference donors per batch rather than one shared donor is the point: a split
    keyed on the donor cannot give one donor's aliquots to more than one split, so a
    single shared reference donor would leave the training split without a control in
    some batches and the fit would be ineligible for a structural reason that has
    nothing to do with CytoNorm.
    """

    samples: list[SampleSpec] = []
    for b_i, batch in enumerate(["B1", "B2", "B3"]):
        for r in range(2):
            samples.append(SampleSpec(
                sample_id=f"ref{b_i}{r}_{batch}", donor_id=f"REF{b_i}{r}", batch=batch,
                group="reference", reference_id="REF1",
                reference_role="technical_replicate", n_cells=n_cells,
            ))
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                samples.append(SampleSpec(
                    sample_id=f"D{b_i}{k}{tag}_{batch}", donor_id=f"D{b_i}{k}{tag}",
                    batch=batch, group=group, n_cells=n_cells,
                ))
    return SimulationSpec(
        samples=samples,
        seed=473,
        batch_shift={"B1": {"CD4": 0.0}, "B2": {"CD4": 0.9}, "B3": {"CD4": -0.7}},
        group_shift={"case": {"CD69": 0.9}},
    )


def _materialise(target: Path) -> dict:
    """A complete project, loaded and frozen the way a run would."""

    paths = write_simulation(_reference_rich_spec(), target)
    thresholds = target / "config" / "thresholds.tsv"
    thresholds.write_text(
        "population\tfeature_id\tsign\tthreshold\tpriority\n"
        "B\tCD19\t+\t1.75\t0\n"
        "Mono\tCD14\t+\t1.75\t1\n"
        "NK\tCD3\t-\t1.75\t2\n"
        "NK\tCD56\t+\t1.65\t2\n"
        "T_CD4\tCD3\t+\t1.75\t3\n"
        "T_CD4\tCD4\t+\t1.75\t3\n"
        "T_CD8\tCD3\t+\t1.75\t4\n"
        "T_CD8\tCD8\t+\t1.65\t4\n",
        encoding="utf-8",
    )
    config = {
        "schema_version": 2,
        "project": {"id": "cytonorm_smoke", "seed": 473},
        "input": {
            "format": "h5ad", "path": "data/preprocessed.h5ad",
            "expression_layer": "preprocessed", "samples": "config/samples.tsv",
            "markers": "config/markers.tsv",
            "processing_manifest": "config/processing.json",
        },
        "preprocessing": {"enabled": False},
        "design": {
            "correction_batch": "batch", "protect": ["group"], "donor": "donor_id",
            "reference_identity": "reference_id", "reference_role": "reference_role",
            "allow_exploratory_override": False, "minimum_batches": 3,
            "require_shared_reference": True,
        },
        "evaluation": {
            "label_source": "gating", "gating_thresholds": "config/thresholds.tsv",
            "representation": {
                "space": "marker", "marker_mask": "evaluation", "scaling": "none",
                "distance": "euclidean", "k": 30,
            },
            "sampling": {
                "unit": "sample_id", "max_cells_per_sample": 1200,
                "min_cells_per_population": 50,
            },
            "partitions": "frozen_common", "uncertainty_unit": "donor_id",
            "policy": "conservative",
        },
        "correction": {
            "control_policy": "measured_only",
            "candidates": [
                {"id": "uncorrected", "backend": "none", "parameters": {}},
                {"id": "cytonorm_reference", "backend": "cytonorm", "parameters": {}},
            ],
        },
        "detectors": {
            "acquisition_drift": True, "composition": True,
            "cycombine_express": False, "cycombine_cluster": False, "cytobatchflagr": False,
        },
        "upstream": {"lock": "config/upstream_lock.yaml", "checkout_dir": None,
                     "rscript": "Rscript"},
        "report": {"offline": True, "redact_identifiers": True,
                   "title": "batchdx-cyto cytonorm smoke"},
        "exports": {"h5ad": True, "fcs": False},
    }
    config_path = target / "config" / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    _ = paths

    loaded = load_config(config_path)
    dataset = load_dataset(loaded, target)
    design = audit_design(
        list(dataset.samples.values()), dataset_id=dataset.manifest.dataset_id,
        correction_factor="batch", protected=["group"], donor_factor="donor_id",
        minimum_batches=3, require_shared_reference=True,
    )
    labels, tier, confidence = resolve_labels(dataset, loaded, target)
    sampling = build_sampling_manifest(dataset, loaded, design, labels, tier, confidence)
    return {"config": loaded, "dataset": dataset, "design": design, "sampling": sampling,
            "path": target}


def _context(project: dict, tmp_path: Path, **kwargs) -> BackendContext:
    return BackendContext(
        design=project["design"], sampling=project["sampling"],
        work_dir=str(tmp_path / "work"), seed=project["config"].project.seed, **kwargs,
    )


# --------------------------------------------------------------------------------------
# Capabilities and refusals — no CytoNorm required
# --------------------------------------------------------------------------------------


def test_cytonorm_is_registered_with_the_capabilities_it_now_has() -> None:
    caps = available_backends()["cytonorm"]
    assert caps.supports_fit_apply is True
    assert "measured_reference_in_every_batch" in caps.reference_layouts
    assert caps.supports_covariates is False
    assert set(caps.accepted_parameters) == {
        "nCells", "xdim", "ydim", "nClus", "nQ", "scale"
    }
    # The transform list is pinned by the adapter, not chosen by a caller.
    assert "transformList" in caps.rejected_parameters
    assert "aggregate_control" in caps.rejected_parameters
    assert "transform the values a second time" in caps.rejected_parameters["transformList"]


def test_the_verified_flag_cannot_be_claimed_without_the_package() -> None:
    """A capability flag is a claim about an executed test, so it is tied to reality.

    If CytoNorm is not installed here, the smoke test below cannot run, and a
    ``fit_apply_verified=True`` would then be an unbacked claim in the capabilities of
    every emitted record.
    """

    if not cytonorm_available() and FIT_APPLY_VERIFIED:
        pytest.fail(
            "adapters/cytonorm.py declares fit_apply_verified=True but CytoNorm is not "
            "installed in this environment, so the smoke test that earns the flag cannot "
            "have run. Set FIT_APPLY_VERIFIED to False."
        )
    if cytonorm_available():
        assert FIT_APPLY_VERIFIED is True, (
            "CytoNorm is installed but the backend does not claim a verified fit/apply; "
            "run the smoke test and record the result rather than leaving the flag unset"
        )


def test_cytonorm_refuses_a_declared_transform_list() -> None:
    caps = get_backend("cytonorm").capabilities
    with pytest.raises(BackendError) as excinfo:
        resolve_parameters(caps, {"transformList": "cytofTransform"})
    assert "transformList" in str(excinfo.value)
    assert "second time" in str(excinfo.value)


def test_the_r_script_pins_the_transform_list_to_the_identity() -> None:
    """The Python-side refusal is mirrored in R, because R is where it would matter."""

    script = (Path(__file__).resolve().parents[2] / "workflow" / "scripts"
              / "cytonorm_correct.R").read_text(encoding="utf-8")
    assert 'identical(transform, "identity")' in script
    assert "refusing transform_list" in script
    assert "transform the values twice" in script
    # And it must actually call the two separable entry points the lock names.
    assert "CytoNorm.train(" in script
    assert "CytoNorm.normalize(" in script


def test_cytonorm_reports_not_available_with_an_actionable_reason(frozen, tmp_path) -> None:
    """The absent-interpreter case is reproduced, not waited for.

    Gating on whether CytoNorm happens to be installed makes this a no-op in exactly
    the environment where the backend works, which is how the unavailable path stays
    untested.
    """

    _, dataset, _, sampling = frozen
    context = BackendContext(
        design=frozen[2], sampling=sampling, work_dir=str(tmp_path / "w"),
        seed=frozen[0].project.seed, rscript="Rscript-does-not-exist",
    )
    with pytest.raises(BackendNotAvailable) as excinfo:
        get_backend("cytonorm").run(dataset, "c", {}, context)
    message = str(excinfo.value)
    assert "not installed" in message
    assert "never substitutes another algorithm" in message
    assert "GPL" in message, "the licence constraint is part of the install instructions"
    _ = dataset


def test_cytonorm_needs_a_split_it_can_fit_on(frozen, tmp_path) -> None:
    """Fit and apply cannot be separated without a partition to fit on."""

    config, dataset, design, sampling = frozen
    for empty in (None, SamplingManifest.model_validate(
        {**json.loads(sampling.to_json()), "splits": {}}
    )):
        context = BackendContext(
            design=design, sampling=empty, work_dir=str(tmp_path / "w"),
            seed=config.project.seed,
        )
        with pytest.raises(BackendError) as excinfo:
            get_backend("cytonorm").run(dataset, "c", {}, context)
        assert "train/tune/validate" in str(excinfo.value)


def test_cytonorm_is_ineligible_when_the_training_split_holds_no_measured_reference(
    frozen, tmp_path
) -> None:
    """An honest refusal, not a fit on data the split holds out.

    The shipped paired-reference example has one reference donor whose three aliquots
    are indivisible, so the training split can end up without a control in any batch.
    CytoNorm cannot fit there, and the correct outcome is to say so.
    """

    config, dataset, design, sampling = frozen
    members = split_members(sampling.splits)
    train_controls = [
        sid for sid in members["train"]
        if dataset.samples[sid].reference_role is not ReferenceRole.NONE
    ]
    if train_controls:
        pytest.skip(
            f"this seed puts reference aliquot(s) {train_controls} in the training split, so "
            "the ineligible path is not exercised here; the eligibility contract is asserted "
            "directly on CytoNormBackend._training_controls below"
        )

    context = BackendContext(
        design=design, sampling=sampling, work_dir=str(tmp_path / "w"),
        seed=config.project.seed,
    )
    with pytest.raises(BackendError) as excinfo:
        get_backend("cytonorm").run(dataset, "c", {}, context)
    message = str(excinfo.value)
    assert "no measured technical reference" in message
    assert "ineligible" in message
    assert "pseudo control is forbidden" in message


def test_the_training_controls_are_the_training_split_measured_references(tmp_path) -> None:
    """The eligibility rule itself, checked directly rather than through a run.

    A run can only exercise it when the package happens to be installed, and the
    package's presence must not decide whether the rule is tested.
    """

    from batchdx.cyto.adapters.cytonorm import CytoNormBackend

    project = _materialise(tmp_path / "proj")
    dataset, sampling = project["dataset"], project["sampling"]
    members = split_members(sampling.splits)
    assert split_ok(members, dataset), "the fixture must produce a usable split"

    backend = CytoNormBackend()
    controls, batches = backend._training_controls(dataset, members, "donor_id")
    assert controls, "the fixture must leave a measured reference in the training split"
    assert set(controls) <= set(members["train"])
    assert not (set(controls) & set(members["validate"]))
    for sid in controls:
        assert dataset.samples[sid].reference_role in {
            ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR
        }
    assert len({batches[s] for s in controls}) >= 2, (
        "CytoNorm needs a per-batch quantile goal, so the fixture must leave at least two "
        f"batches in the training split; it left {sorted({batches[s] for s in controls})}"
    )


def split_ok(members: dict[str, list[str]], dataset) -> bool:
    return all(members[name] for name in SPLIT_NAMES) and set(
        members["train"]
    ).isdisjoint(members["validate"])


# --------------------------------------------------------------------------------------
# The FCS round trip and the event-identity check — no CytoNorm required
# --------------------------------------------------------------------------------------


def test_the_fcs_round_trip_this_backend_depends_on_is_exact_to_the_format_bound(
    tmp_path,
) -> None:
    """CytoNorm is handed files, so the write/read pair has to be sound and stated."""

    rng = np.random.default_rng(11)
    events = np.column_stack([
        np.arange(500, dtype=np.float64),          # the identity channel: exact integers
        rng.normal(size=(500, 6)) * 40.0 - 5.0,   # marker values, a realistic range
    ])
    names = [IDENTITY_CHANNEL, "CD3", "CD4", "CD8", "CD19", "CD14", "CD56"]
    path = write_fcs(
        tmp_path / "S1.fcs", events, names, names, scale="asinh_cofactor5",
        processing_note="round trip", source_sample="S1",
    )
    contents = read_fcs(path)
    error = assert_fcs_round_trip(events, contents, names)
    assert error <= round_trip_tolerance(events)
    # The identity channel is the verification channel, so it must be exact: values
    # under 2**24 are representable in single precision without loss.
    assert np.array_equal(contents.events[:, 0], events[:, 0])
    # flowio lowers keyword names; the scale must still be recoverable.
    assert contents.keywords["batchdx_scale"] == "asinh_cofactor5"


def test_the_round_trip_check_catches_a_channel_mix_up(tmp_path) -> None:
    events = np.column_stack([np.arange(10, dtype=np.float64), np.ones((10, 2))])
    names = [IDENTITY_CHANNEL, "CD3", "CD4"]
    path = write_fcs(
        tmp_path / "S1.fcs", events, names, names, scale="asinh_cofactor5",
        processing_note="round trip",
    )
    contents = read_fcs(path)

    with pytest.raises(FCSError, match="channel order"):
        assert_fcs_round_trip(events, contents, ["CD3", IDENTITY_CHANNEL, "CD4"])
    with pytest.raises(FCSError, match="did not change values|beyond the"):
        assert_fcs_round_trip(events + 0.5, contents, names)
    with pytest.raises(FCSError, match="produced"):
        assert_fcs_round_trip(events[:5], contents, names)


def test_a_reordered_corrected_file_is_refused_rather_than_matched(tmp_path) -> None:
    """The event-identity guarantee, tested by breaking the thing it guarantees.

    If this only ran on a well-behaved CytoNorm it would pass whether or not the check
    exists, so the input is corrupted here instead.
    """

    rng = np.random.default_rng(3)
    events = np.column_stack([
        np.arange(200, dtype=np.float64), rng.normal(size=(200, 3)),
    ])
    names = [IDENTITY_CHANNEL, "CD3", "CD4", "CD8"]
    path = write_fcs(
        tmp_path / "S1.fcs", events, names, names, scale="asinh_cofactor5",
        processing_note="identity",
    )
    expected = events[:, 0].astype(np.int64)
    assert np.array_equal(verify_event_identity(read_fcs(path), expected, "S1"), expected)

    # A permutation: the same events, in the wrong order.  This is the failure that
    # would be invisible if values were matched by position rather than by identity.
    order = rng.permutation(200)
    shuffled = events[order]
    shuffled_path = write_fcs(
        tmp_path / "shuffled.fcs", shuffled, names, names, scale="asinh_cofactor5",
        processing_note="identity",
    )
    with pytest.raises(EventIdentityError) as excinfo:
        verify_event_identity(read_fcs(shuffled_path), expected, "S1")
    assert "does not correspond event by event" in str(excinfo.value)
    assert "first diverges at row" in str(excinfo.value)

    # A dropped event.
    with pytest.raises(EventIdentityError, match="could not be verified"):
        verify_event_identity(read_fcs(path), expected[:100], "S1")

    # A duplicated event: the right count, the wrong events.  Row 100 is overwritten
    # with a copy of row 50, so event 100 disappears and event 50 appears twice.
    duplicated = events.copy()
    duplicated[100] = duplicated[50]
    dup_path = write_fcs(
        tmp_path / "dup.fcs", duplicated, names, names, scale="asinh_cofactor5",
        processing_note="identity",
    )
    with pytest.raises(EventIdentityError) as excinfo:
        verify_event_identity(read_fcs(dup_path), expected, "S1")
    assert "first diverges at row 100" in str(excinfo.value)
    assert issubclass(EventIdentityError, BackendError), (
        "the identity failure must be typed as a backend failure, so the pipeline records "
        "FAILED rather than NOT_AVAILABLE"
    )


# --------------------------------------------------------------------------------------
# The real correction — CytoNorm required
# --------------------------------------------------------------------------------------


@pytest.mark.requires_cytonorm
def test_the_cytonorm_backend_corrects_a_held_out_split_and_verifies_event_identity(
    tmp_path,
) -> None:
    _requires_cytonorm()
    project = _materialise(tmp_path / "proj")
    dataset, sampling = project["dataset"], project["sampling"]
    members = split_members(sampling.splits)
    assert split_ok(members, dataset)

    baseline = dataset.expression()
    marker = "CD4"
    pos = dataset.feature_index([marker])[0]
    batches = {s: dataset.samples[s].batch for s in dataset.samples}
    obs = dataset.obs
    per_sample_batch = obs["sample_id"].map(batches).to_numpy()

    def batch_spread(matrix: np.ndarray) -> float:
        return float(np.std([
            matrix[per_sample_batch == b, pos].mean() for b in sorted(set(batches.values()))
        ]))

    before = batch_spread(baseline)
    matrix, result = get_backend("cytonorm").run(
        dataset, "cytonorm_reference", dict(SMOKE_PARAMS), _context(project, tmp_path)
    )

    assert result.status is CorrectionStatus.SUCCESS
    assert matrix.shape == baseline.shape
    assert not np.array_equal(matrix, baseline), (
        "CytoNorm reported success but returned the input unchanged; that is the one "
        "outcome this adapter must not present as a correction"
    )
    assert result.output_scale == dataset.processing.scale
    assert result.validation["excluded_markers_unchanged"]
    assert result.n_cells == dataset.n_cells

    # The correction did what it claims: the technical offset between batches shrank.
    after = batch_spread(matrix)
    assert after < before, f"batch spread on {marker} did not shrink: {before:.4g} -> {after:.4g}"

    # Markers outside the correction scope are byte-identical, which is only true if
    # the scattered values landed on the right cells.
    out_of_scope = [m for m in dataset.marker_ids() if m not in dataset.marker_ids("correction")]
    for m in out_of_scope:
        col = dataset.feature_index([m])[0]
        assert np.array_equal(matrix[:, col], baseline[:, col]), (
            f"marker {m!r} is outside the correction scope but changed"
        )

    # The fit is recorded, and it is recorded as having used the training split only.
    joined = " ".join(result.warnings)
    assert "training split" in joined
    train_controls = [
        s for s in members["train"]
        if dataset.samples[s].reference_role is not ReferenceRole.NONE
    ]
    assert train_controls, "the fixture must leave a reference in the training split"
    for sid in train_controls:
        assert f"{sid} ({batches[sid]})" in joined, (
            f"the fitted-on control {sid!r} must be named in the record, so a reader can "
            "check for themselves that no validation sample contributed to the fit"
        )
    for sid in members["validate"]:
        assert f"{sid} ({batches[sid]})" not in joined.split("applied to every sample")[0], (
            f"validation sample {sid!r} appears among the fitted-on controls"
        )

    # Upstream's own degradation is surfaced rather than swallowed: a cluster left on
    # the identity function is a cluster that came back unchanged, and the record says
    # so.  This fixture is sized so that it does not happen, which is asserted rather
    # than assumed.
    assert result.resources["n_metaclusters"] > 0
    degenerate = [w for w in result.warnings if "identity function" in w]
    assert degenerate == [], (
        f"CytoNorm fell back to the identity function on this fixture: {degenerate}. Raise "
        "SMOKE_N_CELLS or lower nQ so the quantile goal is estimable; a partial correction "
        "reported as a complete one is the failure this backend exists to avoid"
    )
    assert result.resources["fcs_write_max_abs_error"] >= 0.0
    assert result.resources["max_abs_change"] > 0.0
    assert (Path(_context(project, tmp_path).work_dir) / "cytonorm.log").exists()


@pytest.mark.requires_cytonorm
def test_a_batch_the_split_never_fitted_is_reported_as_uncorrected_not_corrected(
    tmp_path,
) -> None:
    """The one outcome this adapter must not produce silently: unchanged data as success.

    ``CytoNorm.normalize`` skips a file it has no fitted spline for, and its final
    combine step then writes that file back with no cluster replaced — byte identical
    to the input.  So "it wrote a file" and "it corrected the cells" are different
    claims, and only the second may be called a correction.

    The fixture makes this happen for real: the training split cannot hold a reference
    in all three batches, because a split keyed on the donor cannot split one donor's
    aliquots across splits.  B3 therefore lands outside the training split and comes
    back untouched, which is exactly the situation the report has to state.
    """

    _requires_cytonorm()
    project = _materialise(tmp_path / "proj")
    dataset, sampling = project["dataset"], project["sampling"]
    members = split_members(sampling.splits)
    all_batches = sorted({str(dataset.samples[s].batch) for s in dataset.samples})

    matrix, result = get_backend("cytonorm").run(
        dataset, "cytonorm_reference", dict(SMOKE_PARAMS), _context(project, tmp_path)
    )

    # The batches that were *fitted*, not the batches present in the training split:
    # biological samples of B3 are in the training split too, but only a measured
    # reference in the training split produces a quantile model for B3.
    fitted_batches = sorted({
        str(dataset.samples[s].batch) for s in members["train"]
        if dataset.samples[s].reference_role is not ReferenceRole.NONE
    })
    all_batches = sorted({str(dataset.samples[s].batch) for s in dataset.samples})
    unfitted = sorted(set(all_batches) - set(fitted_batches))
    assert unfitted, (
        f"this fixture is meant to leave a batch outside the fitted controls; all of "
        f"{all_batches} were fitted, so the partial-correction path is no longer exercised "
        "by this test"
    )
    assert result.resources["n_uncorrected_samples"] > 0
    assert result.resources["n_fitted_batches"] == float(len(fitted_batches))
    assert result.validation["all_samples_corrected"] is False

    joined = " ".join(result.warnings)
    assert "PARTIAL CORRECTION" in joined
    assert all(b in joined for b in unfitted), f"unfitted batch(es) {unfitted} not named"
    assert "bit-identical to the baseline" in joined

    # The uncorrected samples really were not corrected.  The comparison is against
    # the single-precision bound rather than exact equality, because the values make a
    # round trip through float32 FCS storage on the way in and on the way out; a change
    # smaller than what the format can represent is not a correction.
    obs = dataset.obs
    markers = dataset.marker_ids("correction")
    feature_pos = dataset.feature_index(markers)
    baseline = dataset.expression()
    for sid in sorted(dataset.samples):
        rows = np.flatnonzero(obs["sample_id"].to_numpy() == sid)
        rows = rows[np.argsort(
            obs["original_event_index"].to_numpy().astype(np.int64)[rows], kind="stable"
        )]
        block = baseline[rows][:, feature_pos]
        delta = float(np.max(np.abs(matrix[np.ix_(rows, feature_pos)] - block)))
        in_unfitted = str(dataset.samples[sid].batch) in unfitted
        if in_unfitted:
            assert delta <= round_trip_tolerance(block), (
                f"sample {sid!r} is in unfitted batch {unfitted} but changed by {delta:.4g}, "
                "so the record's claim that it is uncorrected would be false"
            )
        else:
            assert delta > round_trip_tolerance(block) * 10, (
                f"sample {sid!r} is in fitted batch(s) {fitted_batches} but barely changed "
                f"({delta:.4g}); the correction did not happen and the record does not say so"
            )
    # And a batch contrast spanning a corrected and an uncorrected batch is therefore
    # only half corrected, which is what the warning says.
    assert "batch contrast involving those batches is not corrected" in joined


@pytest.mark.requires_cytonorm
def test_the_cytonorm_script_refuses_a_transform_list_the_python_side_cannot_reach(
    tmp_path,
) -> None:
    """The R-side pin is reachable by any caller, so it is tested on its own terms."""

    _requires_cytonorm()
    with pytest.raises(RScriptFailed) as excinfo:
        run_r_script(
            "cytonorm_correct.R",
            {
                "train_files": ["a.fcs"], "train_labels": ["B1"],
                "apply_files": ["a.fcs"], "apply_labels": ["B1"],
                "channels": ["CD3"], "transform_list": "cytofTransform",
                "train_dir": str(tmp_path / "t"), "output_dir": str(tmp_path / "o"),
                "seed": 1, "nCells": 100, "xdim": 2, "ydim": 2, "nClus": 2,
                "nQ": 5, "scale": False,
            },
            tmp_path,
        )
    assert "refusing transform_list" in excinfo.value.log


@pytest.mark.requires_cytonorm
def test_the_cytonorm_script_refuses_a_fit_on_a_single_batch(tmp_path) -> None:
    """One batch has no per-batch quantile goal; that is reported, not worked around."""

    _requires_cytonorm()
    events = np.column_stack([np.arange(50, dtype=np.float64), np.ones((50, 2))])
    names = [IDENTITY_CHANNEL, "CD3", "CD4"]
    write_fcs(
        tmp_path / "a.fcs", events, names, names, scale="asinh_cofactor5",
        processing_note="x",
    )
    with pytest.raises(RScriptFailed) as excinfo:
        run_r_script(
            "cytonorm_correct.R",
            {
                "train_files": [str(tmp_path / "a.fcs")], "train_labels": ["B1"],
                "apply_files": [str(tmp_path / "a.fcs")], "apply_labels": ["B1"],
                "channels": ["CD3", "CD4"], "transform_list": TRANSFORM_LIST,
                "train_dir": str(tmp_path / "t"), "output_dir": str(tmp_path / "o"),
                "seed": 1, "nCells": 100, "xdim": 2, "ydim": 2, "nClus": 2,
                "nQ": 5, "scale": False,
            },
            tmp_path,
        )
    assert "needs at least two" in excinfo.value.log
    assert "not a defect in CytoNorm" in excinfo.value.log


def test_cytonorm_refuses_to_fit_on_a_partition_that_cannot_be_compared(tmp_path) -> None:
    """A fit on an arm-unbalanced partition produces a held-out number for nothing.

    The sampling manifest tolerates a design whose thinner arm holds fewer biospecimens
    than there are splits, because it also serves detectors and metrics.  A fit-and-apply
    backend may not: if the validation split cannot hold one arm, the contrast it would
    report as held out does not exist inside the held-out data.  This runs without
    CytoNorm, because the refusal happens before the package is called.
    """

    from batchdx.cyto.adapters.cytonorm import CytoNormBackend

    project = _materialise(tmp_path / "proj")
    dataset, sampling = project["dataset"], project["sampling"]

    # A partition that drops one arm from the validation split, as a design too thin to
    # stratify would produce.
    members = split_members(sampling.splits)
    unbalanced = {
        "split_unit": [split_unit_field(sampling.splits) or "donor_id"],
        "train": members["train"],
        "tune": members["tune"],
        "validate": members["validate"][:1],
    }
    thin = {s: dataset.samples[s] for s in unbalanced["validate"]}
    if len({dataset.samples[s].group for s in thin}) > 1:
        unbalanced["validate"] = [
            s for s in members["validate"] if dataset.samples[s].group == "case"
        ]
    if len(unbalanced["validate"]) == len(members["validate"]):
        pytest.skip("this fixture's validation split is single-arm already")

    assert missing_arms(unbalanced, dataset.samples), "the fixture must defeat the check"
    context = BackendContext(
        design=project["design"],
        sampling=sampling.model_copy(update={"splits": unbalanced}),
        work_dir=str(tmp_path / "w"), seed=473,
    )
    with pytest.raises(BackendError) as excinfo:
        CytoNormBackend()._split(dataset, context)
    message = str(excinfo.value)
    assert "not arm-balanced" in message
    assert "may not fit on a partition" in message
    # It is typed as not-available, not as a failure: the algorithm is fine, the data
    # is not.
    assert isinstance(excinfo.value, BackendNotAvailable)


def test_a_single_valued_r_payload_field_does_not_degrade_into_its_characters() -> None:
    """A bare string is not a one-element list, and iterating one iterates characters.

    ``jsonlite`` with ``auto_unbox = TRUE`` collapses a length-one character vector to a
    bare string.  Reading that as a list gives ``['3', 'B']`` for the batch name ``"B3"``,
    which is a wrong answer that looks right: the warning names batches that do not
    exist and omits the one that does.  The normaliser is the only thing standing
    between that and the report, so it is tested directly.
    """

    from batchdx.cyto.adapters.cytonorm import _as_list

    assert _as_list("B3") == ["B3"]
    assert _as_list(["B3"]) == ["B3"]
    assert _as_list(["B2", "B3"]) == ["B2", "B3"]
    assert _as_list(None) == []
    assert _as_list([]) == []
    # The exact shape jsonlite produces for a length-one vector.
    assert _as_list(json.loads('{"x": "B3"}')["x"]) == ["B3"]
