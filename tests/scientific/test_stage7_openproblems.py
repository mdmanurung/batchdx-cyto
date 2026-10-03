"""Stage 7, continued: the OpenProblems metric rows, reconciled against the pinned commit.

The nine metric names in the OpenProblems cytometry task at commit ``37151479`` are
reconciled one by one here.  Three of them are *extracted* -- the pinned Viash helper is
imported, run and compared -- and those tests are the only evidence behind the
``upstream_verified=True`` the rows carry.  ``scripts/audit_upstream.py`` fails the
build if such a row exists without the flag, so the flag and these tests have to agree.

Where the pinned checkout has not been fetched the parity tests skip with a visible
reason rather than passing vacuously, following ``test_stage0_audit.py``.  Everything
else -- the accounting, the record statuses, the case classifier, the frozen-partition
provenance -- runs either way, so the reconciliation cannot quietly become untested.
"""

from __future__ import annotations

import importlib.util
from dataclasses import replace
from pathlib import Path
from types import ModuleType

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import stats

from batchdx.cyto.metrics import openproblems
from batchdx.cyto.metrics.openproblems import (
    CLISI_SEED,
    EMD_MIN_CELLS_PER_SIDE,
    ILISI_SEED,
    LISI_NEIGHBOURS,
    PEAK_MIN_CELLS_PER_SIDE,
    R2_MIN_CELLS_PER_SIDE,
    _emd_records,
    _Panel,
    _peaks_unless_constant,
    bin_array,
    binned_emd,
    compute_openproblems_metrics,
    count_peaks,
    kde_density,
    standardise_marker_expression,
)
from batchdx.cyto.metrics.preservation import summarise
from batchdx.cyto.metrics.registry import (
    DEFERRED_METRICS,
    OPENPROBLEMS_METRIC_NAMES,
    REIMPLEMENTED_UPSTREAM_METRICS,
    UPSTREAM_METRICS_NOT_FOUND,
    registered_metrics,
)
from batchdx.cyto.schemas.contracts import RepresentationSettings, SamplingManifest
from batchdx.cyto.schemas.enums import LabelTier, MetricStatus

UPSTREAM_METRICS = "task_cyto_batch_integration/src/metrics"
PINNED_COMMIT = "37151479"

#: The metric names registered as extracted, and the pinned helper each is verified
#: against.  Written out rather than derived, so that a name quietly leaving the
#: namespace or losing its flag fails here instead of passing unnoticed.
EXTRACTED = {
    "average_batch_r2_ct": "average_batch_r2/helper.py",
    "cLisi": "lisi/helper.py",
    "emd_mean_ct_horiz": "emd/helper.py",
    "emd_mean_ct_vert": "emd/helper.py",
    "iLisi": "lisi/helper.py",
    "ratio_consistent_peaks": "ratio_consistent_peaks/helper.py",
}


def _pinned_helper(upstream_dir: Path, relative: str) -> ModuleType:
    """Import a pinned Viash ``helper.py`` by path, the way its component would."""

    path = upstream_dir / UPSTREAM_METRICS / relative
    if not path.exists():
        pytest.skip(
            f"the pinned OpenProblems checkout is absent ({path}); "
            "run scripts/fetch_upstream.sh"
        )
    spec = importlib.util.spec_from_file_location(
        f"pinned_{relative.replace('/', '_')}", path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sampling_manifest_stub() -> SamplingManifest:
    """A minimal frozen manifest, for the tests that exercise a classifier directly."""

    return SamplingManifest(
        dataset_id="fixture", seed=0, unit="sample", max_cells_per_sample=None,
        label_source=LabelTier.GATING, label_confidence="medium", label_column="label",
        representation=RepresentationSettings(),
    )


def _panel(values: np.ndarray, markers: list[str] | None = None) -> _Panel:
    """One donor, two batches, two samples, one population: the shape every unit test uses."""

    n = values.shape[0]
    half = n // 2
    return _Panel(
        expression=values,
        markers=markers or [f"M{i}" for i in range(values.shape[1])],
        population=np.array(["T"] * n, dtype=object),
        sample=np.array(["s1"] * half + ["s2"] * (n - half), dtype=object),
        donor=np.array(["D1"] * n, dtype=object),
        batch=np.array(["B1"] * half + ["B2"] * (n - half), dtype=object),
        rows=np.arange(n),
    )


def _with_expression(panel: _Panel, values: np.ndarray) -> _Panel:
    return _Panel(
        expression=values, markers=panel.markers, population=panel.population,
        sample=panel.sample, donor=panel.donor, batch=panel.batch, rows=panel.rows,
    )


def _bimodal(rng: np.random.Generator, n: int) -> np.ndarray:
    """Two well-separated modes, so the KDE has two peaks above the pinned thresholds."""

    half = n // 2
    return np.concatenate([rng.normal(loc=-3.0, size=half), rng.normal(loc=3.0, size=n - half)])


def _aggregate_record(records, metric_name: str):
    matches = [r for r in records if r.metric_name == metric_name and r.comparison_id is None]
    assert len(matches) == 1, f"expected exactly one aggregate row for {metric_name}"
    return matches[0]


def _aggregate_value(records, metric_name: str) -> float | None:
    return _aggregate_record(records, metric_name).value


# --------------------------------------------------------------------------------------
# The accounting: nine names, each resolved exactly once
# --------------------------------------------------------------------------------------


def test_every_openproblems_metric_name_is_resolved_exactly_once() -> None:
    """No name silently dropped, and none resolved twice.

    Each of the nine is either a registered ``compatibility`` spec under its own name, a
    name in ``REIMPLEMENTED_UPSTREAM_METRICS`` whose target is registered, or a row in
    ``DEFERRED_METRICS``.  Two of those at once is a contradiction -- a metric cannot be
    both deferred and registered -- and none of them is a row somebody forgot.
    """

    registry = registered_metrics()
    resolved: dict[str, list[str]] = {}
    for name in OPENPROBLEMS_METRIC_NAMES:
        where: list[str] = []
        bare = name.split(".", 1)[1]
        if ("compatibility", bare) in registry:
            where.append("registered")
        if name in REIMPLEMENTED_UPSTREAM_METRICS:
            where.append("reimplemented")
        if name in DEFERRED_METRICS:
            where.append("deferred")
        resolved[name] = where

    for name, where in resolved.items():
        assert len(where) == 1, f"{name} is resolved {where}; exactly one is required"


def test_a_reimplementation_is_registered_under_a_name_of_its_own() -> None:
    """The local spec has to exist, and it may not borrow the upstream metric's name.

    The mapping is empty at present, because the only reimplementation it held,
    ``average_batch_r2_ct``, became verifiable when scikit-learn was declared and was
    promoted.  The invariants still have to hold for the next one, so the loop is kept
    and the emptiness is asserted rather than assumed.
    """

    registry = registered_metrics()
    assert REIMPLEMENTED_UPSTREAM_METRICS == {}, (
        "a reimplementation is registered; check it against the rules below rather than "
        "letting this test rot"
    )
    for upstream_name, local in REIMPLEMENTED_UPSTREAM_METRICS.items():
        namespace, _, name = local.partition(".")
        assert (namespace, name) in registry, f"{local} is claimed but not registered"
        assert namespace == "common", f"{local} claims a namespace it was not verified in"
        assert name != upstream_name.split(".", 1)[1], (
            "reusing the upstream name would be a compatibility claim the reimplementation "
            "cannot support"
        )
        assert registry[(namespace, name)].upstream_verified is False


def test_the_lock_and_the_registry_agree_on_the_nine_names(repo_root: Path) -> None:
    lock = yaml.safe_load((repo_root / "config" / "upstream_lock.yaml").read_text("utf-8"))
    present = lock["resources"]["openproblems_cyto_batch_integration"][
        "metrics_present_at_commit"
    ]
    assert sorted(present) == sorted(n.split(".", 1)[1] for n in OPENPROBLEMS_METRIC_NAMES)


def test_cms_stays_recorded_as_absent_at_the_pinned_commit() -> None:
    """The CellMixS metric is not there, and that finding must not be quietly retired.

    It appears in the CHANGELOG but ``src/metrics/cms`` does not exist at ``37151479``,
    so the honest entry is "unavailable" -- not "deferred", and never a zero.
    """

    assert set(UPSTREAM_METRICS_NOT_FOUND) == {"openproblems.cms"}
    assert "absent at the pinned commit" in UPSTREAM_METRICS_NOT_FOUND["openproblems.cms"]
    assert "cms" not in {n.split(".", 1)[1] for n in OPENPROBLEMS_METRIC_NAMES}
    assert "openproblems.cms" not in DEFERRED_METRICS
    assert "openproblems.cms" not in REIMPLEMENTED_UPSTREAM_METRICS
    assert not any(name == "cms" for _, name in registered_metrics())


def test_each_extracted_row_names_its_pinned_helper_and_claims_verification() -> None:
    registry = registered_metrics()
    for name, helper in EXTRACTED.items():
        spec = registry[("compatibility", name)]
        assert spec.upstream_verified, name
        assert spec.namespace == "compatibility"
        assert PINNED_COMMIT in spec.version
        assert helper.split("/")[0] in spec.implementation
        assert f"openproblems.{name}" in OPENPROBLEMS_METRIC_NAMES


def test_each_spec_states_the_cell_minimum_the_code_actually_enforces() -> None:
    """A reader has to be able to tell from the spec which groups were excluded.

    The threshold in the spec text and the threshold in the computation are the same
    number, so the two cannot drift apart and quietly contradict each other.
    """

    minimums = {
        "average_batch_r2_ct": R2_MIN_CELLS_PER_SIDE,
        "cLisi": LISI_NEIGHBOURS + 1,
        "emd_mean_ct_horiz": EMD_MIN_CELLS_PER_SIDE,
        "emd_mean_ct_vert": EMD_MIN_CELLS_PER_SIDE,
        "iLisi": LISI_NEIGHBOURS + 1,
        "ratio_consistent_peaks": PEAK_MIN_CELLS_PER_SIDE,
    }
    registry = registered_metrics()
    for name, minimum in minimums.items():
        spec = registry[("compatibility", name)]
        assert str(minimum) in spec.eligible_population, (name, spec.eligible_population)


def test_the_promoted_row_carries_the_evidence_instead_of_an_apology() -> None:
    """``average_batch_r2_ct`` was a reimplementation; now it is a verified component.

    The spec has to say what was checked and what remains an adaptation, and it must no
    longer claim to be unreconciled -- an unreconciled row that carries a verification
    flag is exactly the kind of contradiction the audit cannot see.
    """

    spec = registered_metrics()[("compatibility", "average_batch_r2_ct")]
    assert spec.upstream_verified is True
    assert spec.namespace == "compatibility"
    assert PINNED_COMMIT in spec.version
    assert "reimplementation" not in spec.notes.lower()
    assert "not reconciled" not in spec.notes.lower()
    assert "ADAPTATION, not a copy" in spec.notes
    # and it must no longer be reachable under the old local name
    assert not any(name == "residual_batch_r2" for _, name in registered_metrics())


def test_every_remaining_deferred_row_names_the_obstacle_instead_of_a_tier() -> None:
    """A row that only says "Tier B, not implemented" has told nobody anything.

    The keyword per row is the mechanism that could not be resolved, not a paraphrase
    of the row itself.
    """

    expected = {
        "openproblems.flowsom_mean_mapping_similarity": "FlowSOM",
        "openproblems.functional_marker_preservation_wilcoxon": "retention ratio",
        "openproblems.functional_marker_preservation_cohens_d": "drift",
    }
    for name, keyword in expected.items():
        assert name in DEFERRED_METRICS, name
        assert keyword.lower() in DEFERRED_METRICS[name].lower(), (name, keyword)
        # the placeholder every one of these rows used to carry
        assert "Tier B" not in DEFERRED_METRICS[name], name


def test_the_functional_marker_rows_are_not_claimed_as_local_implementations() -> None:
    """The local statistics are the *inputs* to those metrics, not the metrics.

    Reading the pinned ``script.py`` is what established this: the Wilcoxon row is a
    fraction of significant pairs retained, and the Cohen's d row is effect-size drift
    against an uncorrected baseline.  Neither is the per-pair statistic Stage 8 computes.
    """

    for name in ("openproblems.functional_marker_preservation_wilcoxon",
                 "openproblems.functional_marker_preservation_cohens_d"):
        assert name not in REIMPLEMENTED_UPSTREAM_METRICS
        assert name in DEFERRED_METRICS
    registry = registered_metrics()
    assert ("common", "biological_contrast_rank_statistic") in registry
    assert ("common", "biological_contrast_cohens_d") in registry


# --------------------------------------------------------------------------------------
# Parity with the pinned emd/helper.py
# --------------------------------------------------------------------------------------


def test_the_pinned_bin_array_reproduces_the_extracted_one(upstream_dir: Path) -> None:
    pinned = _pinned_helper(upstream_dir, "emd/helper.py")
    rng = np.random.default_rng(11)
    cases = (
        rng.normal(size=500),
        rng.normal(loc=4.0, scale=0.2, size=500),
        np.concatenate([np.zeros(300), rng.normal(loc=1.0, size=200)]),
    )
    for values in cases:
        want_support, want_weights = pinned.bin_array(values)
        got_support, got_weights = bin_array(values)
        assert np.array_equal(got_support, want_support)
        assert np.array_equal(got_weights, want_weights)
        assert got_support.shape == (2000,), "the pinned grid is 2000 bins"
        assert np.isclose(got_weights.sum(), 1.0)
    # the 1e-7 edge shift the pinned helper adds, which keeps the left edge itself out
    # of the grid: a value at exactly -100.0 lands in no bin at all
    support, weights = bin_array(np.array([-100.0, -99.99]))
    assert weights[0] == 1.0
    assert np.isclose(weights.sum(), 1.0)
    assert support[0] == pytest.approx(-99.95, abs=1e-6)


def test_the_pinned_compute_emd_reproduces_the_extracted_distance(upstream_dir: Path) -> None:
    """The pinned per-marker table and the port, over the same expression values."""

    pinned = _pinned_helper(upstream_dir, "emd/helper.py")
    rng = np.random.default_rng(12)
    markers = ["CD3", "CD4", "CD8"]
    left = rng.normal(size=(400, 3))
    right = rng.normal(loc=1.5, size=(400, 3))
    view = ad.AnnData(
        np.concatenate([left, right]).astype("float32"),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(800)]),
        var=pd.DataFrame(index=pd.Index(markers, dtype=object)),
    )
    view.layers["integrated"] = view.X.copy()

    want = pinned.compute_emd(
        left_sample=view[:400], right_sample=view[400:], markers_to_assess=markers
    )
    for k, marker in enumerate(markers):
        assert binned_emd(left[:, k], right[:, k]) == pytest.approx(
            float(want[marker].iloc[0]), rel=1e-9, abs=1e-12
        )
    # one bin width of separation is one unit of this metric's scale
    assert binned_emd(left[:, 0], right[:, 0] + 0.1) - binned_emd(left[:, 0], right[:, 0]) \
        == pytest.approx(0.1, rel=1e-6, abs=1e-9)


def _shared_split_fixture() -> tuple[ad.AnnData, ad.AnnData, list[str], np.ndarray, np.ndarray]:
    """Two AnnData splits shaped the way the pinned emd component expects them.

    ``group`` is set equal to ``batch`` on purpose: that is exactly the substitution
    documented for the frozen partition, and it makes the pinned component's own
    grouping (every pair of samples within a group) and the local one (every pair of
    samples within a batch) the same set of comparisons, so the two runs can be
    compared number for number rather than approximately.
    """

    rng = np.random.default_rng(13)
    markers = ["CD3", "CD4", "CD8", "CD19"]
    cell_types = ["T_CD4", "Mono", "NK"]
    donors = ["D1", "D2", "D3"]
    per_group = 140
    blocks: list[tuple[str, str, str, str]] = []
    frames: list[np.ndarray] = []
    for batch in ("B1", "B2"):
        offset = 0.0 if batch == "B1" else 0.8
        for sample_index, donor in enumerate(donors):
            for replicate in range(2):
                for cell_type in cell_types:
                    frames.append(
                        rng.normal(
                            loc=offset + 0.4 * cell_types.index(cell_type), scale=1.0,
                            size=(per_group, len(markers)),
                        )
                    )
                    blocks.append(
                        (donor, batch, f"{batch}_s{sample_index}_{replicate}", cell_type)
                    )
    expression = np.concatenate(frames)
    # one obs row per cell, in the order the frames were concatenated
    per_cell = [block for block in blocks for _ in range(per_group)]
    obs = pd.DataFrame(
        {
            "donor": [b[0] for b in per_cell], "batch": [b[1] for b in per_cell],
            "sample": [b[2] for b in per_cell], "group": [b[1] for b in per_cell],
            "cell_type": [b[3] for b in per_cell],
        },
        index=[f"cell{i}" for i in range(len(per_cell))],
    )

    def _view(mask: np.ndarray) -> ad.AnnData:
        subset = ad.AnnData(
            expression[mask].astype("float32"), obs=obs.loc[mask].copy(),
            var=pd.DataFrame(index=pd.Index(markers, dtype=object)),
        )
        subset.layers["integrated"] = subset.X.copy()
        return subset

    first = obs["batch"].to_numpy() == "B1"
    return _view(first), _view(~first), markers, expression[first], expression[~first]


def _panel_from_expression(
    expression: np.ndarray, markers: list[str], obs: pd.DataFrame
) -> _Panel:
    return _Panel(
        expression=expression, markers=list(markers),
        population=obs["cell_type"].to_numpy().astype(object),
        sample=obs["sample"].to_numpy().astype(object),
        donor=obs["donor"].to_numpy().astype(object),
        batch=obs["batch"].to_numpy().astype(object),
        rows=np.arange(expression.shape[0]),
    )


def test_the_pinned_horizontal_emd_matches_the_extracted_component(upstream_dir: Path) -> None:
    """Both implementations run on the same fixture and their means must agree.

    The pinned component compares a donor's cells in split 1 with the same donor's cells
    in split 2; the local one compares a donor's batch A cells with its batch B cells.
    On this fixture those are the same comparison, and the recorded tolerance is 1e-9
    relative -- float summation order over the pooled markers, and nothing else.
    """

    pinned = _pinned_helper(upstream_dir, "emd/helper.py")
    split1, split2, markers, first, second = _shared_split_fixture()
    obs1, obs2 = split1.obs, split2.obs
    panel = _Panel(
        expression=np.concatenate([first, second]), markers=markers,
        population=np.concatenate([
            obs1["cell_type"].to_numpy().astype(object),
            obs2["cell_type"].to_numpy().astype(object),
        ]),
        sample=np.concatenate([
            obs1["sample"].to_numpy().astype(object),
            obs2["sample"].to_numpy().astype(object),
        ]),
        donor=np.concatenate([
            obs1["donor"].to_numpy().astype(object),
            obs2["donor"].to_numpy().astype(object),
        ]),
        batch=np.concatenate([
            obs1["batch"].to_numpy().astype(object),
            obs2["batch"].to_numpy().astype(object),
        ]),
        rows=np.arange(first.shape[0] + second.shape[0]),
    )

    want = pinned.calculate_horizontal_emd(
        i_split1_adata=split1, i_split2_adata=split2, markers_to_assess=markers,
        donor_list=list(split1.obs["donor"].unique()),
    )[pinned.KEY_MEAN_EMD_CT]
    got = _aggregate_value(
        openproblems._horizontal_emd(panel, _sampling_manifest_stub(), "run", "cand"),
        "emd_mean_ct_horiz",
    )

    assert got == pytest.approx(float(want), rel=1e-9, abs=1e-12)
    assert got > 0.0, "the fixture carries a real batch offset, so this cannot be zero"


def test_the_pinned_vertical_emd_matches_the_extracted_component(upstream_dir: Path) -> None:
    """Same fixture, same comparison set, same mean."""

    pinned = _pinned_helper(upstream_dir, "emd/helper.py")
    split1, split2, markers, first, second = _shared_split_fixture()
    obs1, obs2 = split1.obs, split2.obs
    panel = _Panel(
        expression=np.concatenate([first, second]), markers=markers,
        population=np.concatenate([
            obs1["cell_type"].to_numpy().astype(object),
            obs2["cell_type"].to_numpy().astype(object),
        ]),
        sample=np.concatenate([
            obs1["sample"].to_numpy().astype(object),
            obs2["sample"].to_numpy().astype(object),
        ]),
        donor=np.concatenate([
            obs1["donor"].to_numpy().astype(object),
            obs2["donor"].to_numpy().astype(object),
        ]),
        batch=np.concatenate([
            obs1["batch"].to_numpy().astype(object),
            obs2["batch"].to_numpy().astype(object),
        ]),
        rows=np.arange(first.shape[0] + second.shape[0]),
    )

    want = pinned.calculate_vertical_emd(
        i_split1_adata=split1, i_split2_adata=split2, markers_to_assess=markers
    )[pinned.KEY_MEAN_EMD_CT]
    got = _aggregate_value(
        openproblems._vertical_emd(panel, _sampling_manifest_stub(), "run", "cand"),
        "emd_mean_ct_vert",
    )

    assert got == pytest.approx(float(want), rel=1e-9, abs=1e-12)


# --------------------------------------------------------------------------------------
# Parity with the pinned ratio_consistent_peaks/helper.py, and the case transcription
# --------------------------------------------------------------------------------------


def test_the_pinned_peak_helpers_reproduce_the_extracted_ones(upstream_dir: Path) -> None:
    """Both helpers run on the same arrays, so the peak counts are identical."""

    pinned = _pinned_helper(upstream_dir, "ratio_consistent_peaks/helper.py")
    rng = np.random.default_rng(14)
    cases = (
        rng.normal(size=400),
        _bimodal(rng, 400),
        np.concatenate([rng.normal(size=300), rng.normal(loc=6.0, size=100)]),
    )
    for values in cases:
        a, b = values, values * 2.0 + 1.0
        for got, expected in zip(
            standardise_marker_expression(a, b),
            pinned.standardise_marker_expression(a, b),
            strict=True,
        ):
            assert np.array_equal(got, expected)
        assert np.allclose(kde_density(values), pinned.get_kde_density(values),
                           rtol=0, atol=1e-12)
        assert count_peaks(kde_density(values)) == pinned.call_peaks(
            pinned.get_kde_density(values)
        )
    # the pinned helper cannot peak-call a constant side at all -- its covariance is
    # singular -- which is why the component guards for that before calling it, and why
    # the local guard reports the triple instead of raising
    with pytest.raises(np.linalg.LinAlgError):
        pinned.get_kde_density(np.zeros(400))
    assert _peaks_unless_constant(np.zeros(400)) is None


def test_the_peak_case_classifier_calls_a_flattened_marker_broken(upstream_dir: Path) -> None:
    """The transcription of the case logic, exercised on a triple that changes.

    The pinned component's own script cannot be run here -- it reads ``.h5ad`` resources
    the checkout does not carry -- so the classifier is tested directly: a correction
    that collapses a bimodal marker to one mode is case 2, and one that leaves the shape
    alone is case 1.
    """

    pinned = _pinned_helper(upstream_dir, "ratio_consistent_peaks/helper.py")
    manifest = _sampling_manifest_stub()
    rng = np.random.default_rng(15)
    baseline = np.concatenate([_bimodal(rng, 200), _bimodal(rng, 200)])
    panel = _panel(baseline[:, None], ["CD3"])

    # the counts the classifier will use are the pinned helper's own counts, taken on
    # the same pooled standardisation the component applies
    scaled, _ = pinned.standardise_marker_expression(baseline[:200], baseline[200:])
    assert pinned.call_peaks(pinned.get_kde_density(scaled)) == 2
    assert openproblems._standardised_peaks(
        {"B1": baseline[:200], "B2": baseline[200:]}
    ) == {"B1": 2, "B2": 2}

    flattened = _panel(rng.normal(size=400)[:, None], ["CD3"])
    broken = openproblems._ratio_consistent_peaks(panel, flattened, manifest, "r", "c")
    assert _aggregate_value(broken, "ratio_consistent_peaks") == 0.0
    assert any("case 2 (broken)" in (record.reason or "") for record in broken)

    shifted = _with_expression(panel, (baseline + 3.0)[:, None])
    preserved = openproblems._ratio_consistent_peaks(panel, shifted, manifest, "r", "c")
    assert _aggregate_value(preserved, "ratio_consistent_peaks") == 1.0
    assert any("case 1 (preserved)" in (record.reason or "") for record in preserved)


def test_a_triple_that_already_disagreed_is_excluded_from_the_ratio() -> None:
    """Upstream's case NGT: the uncorrected data was already inconsistent.

    It must count neither as a success nor as a failure, or the ratio would reward a
    correction for a problem it did not cause -- and with nothing left to divide, the
    metric is undefined rather than 1.0.
    """

    rng = np.random.default_rng(17)
    baseline = np.concatenate([_bimodal(rng, 200), rng.normal(size=200)])
    assert openproblems._standardised_peaks(
        {"B1": baseline[:200], "B2": baseline[200:]}
    ) == {"B1": 2, "B2": 1}
    panel = _panel(baseline[:, None], ["CD3"])
    candidate = _with_expression(panel, (baseline + 1.0)[:, None])
    records = openproblems._ratio_consistent_peaks(
        panel, candidate, _sampling_manifest_stub(), "r", "c"
    )
    aggregate = _aggregate_record(records, "ratio_consistent_peaks")
    assert aggregate.value is None
    assert aggregate.status is MetricStatus.UNDEFINED
    assert any("case NGT" in (record.reason or "") for record in records)


def test_a_triple_the_correction_flattened_is_broken_not_excluded() -> None:
    """Upstream's zero-variance guard: a flat corrected side beside a structured one."""

    rng = np.random.default_rng(18)
    baseline = np.concatenate([_bimodal(rng, 200), _bimodal(rng, 200)])
    panel = _panel(baseline[:, None], ["CD3"])
    flattened = _with_expression(panel, np.full(400, 5.0)[:, None])
    records = openproblems._ratio_consistent_peaks(
        panel, flattened, _sampling_manifest_stub(), "r", "c"
    )
    assert _aggregate_value(records, "ratio_consistent_peaks") == 0.0
    assert any("case 2 (flattened)" in (record.reason or "") for record in records)
    assert _peaks_unless_constant(np.full(200, 3.0)) is None
    assert _peaks_unless_constant(_bimodal(rng, 200)) == 2


# --------------------------------------------------------------------------------------
# The functional-marker helpers the local statistics share with the pinned component
# --------------------------------------------------------------------------------------


def _per_sample_means(rng: np.random.Generator) -> pd.DataFrame:
    """The frame ``functional_marker_preservation/helper.py`` consumes."""

    rows = []
    for marker in ("CD25", "CD69"):
        for cell_type in ("T_CD4", "Mono"):
            for sample_index in range(4):
                group = "A" if sample_index < 2 else "B"
                rows.append({
                    "cell_type": cell_type, "sample": f"s{sample_index}", "group": group,
                    "batch": "B1", "marker": marker,
                    "mean_expr": float(rng.normal(loc=1.0 if group == "A" else 0.0)),
                })
    return pd.DataFrame(rows)


def test_the_pinned_wilcoxon_helper_matches_the_local_rank_statistic(
    upstream_dir: Path,
) -> None:
    """The statistic Stage 8 computes is the pinned component's, on the pinned input.

    The pinned component's *metric* is a retention ratio and is still deferred; the
    Mann-Whitney U that feeds it is not, and this is the evidence that
    ``common.biological_contrast_rank_statistic`` computes the same number.
    """

    pinned = _pinned_helper(upstream_dir, "functional_marker_preservation/helper.py")
    frame = _per_sample_means(np.random.default_rng(19))
    groups = ad.AnnData(
        np.zeros((len(frame), 1), dtype="float32"),
        obs=frame[["group"]].astype("category"),
        var=pd.DataFrame(index=pd.Index(["CD3"], dtype=object)),
    )
    group_a, group_b = pinned.infer_groups(groups)
    assert (group_a, group_b) == ("A", "B"), "the pinned helper sorts, so the sign is stable"
    want = pinned.run_wilcoxon_group_tests(frame, group_a, group_b)

    for marker in ("CD25", "CD69"):
        for cell_type in ("T_CD4", "Mono"):
            subset = frame[(frame["marker"] == marker) & (frame["cell_type"] == cell_type)]
            a = subset[subset["group"] == group_a]["mean_expr"].to_numpy()
            b = subset[subset["group"] == group_b]["mean_expr"].to_numpy()
            row = want[(want["marker"] == marker) & (want["cell_type"] == cell_type)]
            assert len(row) == 1

            local = summarise(np.column_stack([a, b]), model="unpaired")
            assert local.rank_statistic == pytest.approx(
                float(row["u_statistic"].iloc[0]), rel=1e-12
            )
            assert local.rank_pvalue == pytest.approx(
                float(row["pval"].iloc[0]), rel=1e-12
            )
            assert local.rank_statistic == pytest.approx(
                float(stats.mannwhitneyu(a, b, alternative="two-sided").statistic),
                rel=1e-12,
            )


def test_the_pinned_cohens_d_helper_matches_the_local_cohens_d(upstream_dir: Path) -> None:
    """The same pooled-SD form, verified against the pinned helper.

    The metric the helper feeds is effect-size drift against an uncorrected baseline,
    which is still deferred; what is checked here is the per-pair statistic, because
    that is the one the local metric claims to compute.
    """

    pinned = _pinned_helper(upstream_dir, "functional_marker_preservation/helper.py")
    rng = np.random.default_rng(20)
    a = rng.normal(loc=1.0, scale=2.0, size=6)
    b = rng.normal(loc=0.0, scale=1.0, size=5)
    want = pinned._cohens_d(a, b)

    pooled = np.sqrt(
        ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
        / (len(a) + len(b) - 2)
    )
    assert want == pytest.approx((a.mean() - b.mean()) / pooled, rel=1e-12)
    # the donor-level matrix pads the shorter arm with NaN, which is what the two
    # levels of unequal size look like there
    matrix = np.full((max(len(a), len(b)), 2), np.nan)
    matrix[: len(a), 0], matrix[: len(b), 1] = a, b
    assert summarise(matrix, model="unpaired").cohens_d == pytest.approx(float(want), rel=1e-9)
    assert "openproblems.functional_marker_preservation_cohens_d" in DEFERRED_METRICS


# --------------------------------------------------------------------------------------
# The frozen partition: statuses and provenance
# --------------------------------------------------------------------------------------


def test_every_record_carries_the_frozen_partition_provenance(frozen) -> None:
    _, dataset, design, sampling = frozen
    result = compute_openproblems_metrics(
        dataset, dataset.expression(), sampling, design,
        run_id="run-openproblems", candidate_id="uncorrected",
    )
    assert result.records
    for record in result.records:
        assert record.label_tier is sampling.label_source, record.metric_name
        assert record.label_confidence == sampling.label_confidence
        assert record.representation == sampling.representation
        assert record.run_id == "run-openproblems"
        assert record.candidate_id == "uncorrected"
        assert record.status in set(MetricStatus)


def test_a_thin_comparison_reports_a_reason_and_never_a_zero(frozen) -> None:
    """A group below the cell minimum is ``insufficient_data``, not a score of zero.

    Zero would say the two batches agree perfectly, which was never measured.
    """

    _, dataset, design, sampling = frozen
    result = compute_openproblems_metrics(
        dataset, dataset.expression(), sampling, design, run_id="r", candidate_id="c",
    )
    minimums = {
        "average_batch_r2_ct": R2_MIN_CELLS_PER_SIDE,
        "cLisi": LISI_NEIGHBOURS + 1,
        "emd_mean_ct_horiz": EMD_MIN_CELLS_PER_SIDE,
        "emd_mean_ct_vert": EMD_MIN_CELLS_PER_SIDE,
        "iLisi": LISI_NEIGHBOURS + 1,
        "ratio_consistent_peaks": PEAK_MIN_CELLS_PER_SIDE,
    }
    thin = [r for r in result.records if r.status is MetricStatus.INSUFFICIENT_DATA]
    assert thin, "the fixture has populations well below the cell minimum on one side"
    for record in thin:
        assert record.value is None
        assert record.reason
        assert str(minimums[record.metric_name]) in record.reason, record.reason

    # a metric that was refused must not be reported as a score of zero either
    refused = [r for r in result.records if r.status is MetricStatus.UNDEFINED]
    for record in refused:
        assert record.value is None
        assert record.reason
    assert refused, (
        "this fixture confounds at least one population with a single batch, so iLisi is "
        "withheld; if that ever stops being true the test is no longer checking the "
        "refusal path and should say so"
    )


def test_a_design_with_no_within_donor_pair_says_so_instead_of_scoring_zero(
    frozen,
) -> None:
    """Donors nested inside batches give the horizontal component nothing to compare."""

    _, dataset, design, sampling = frozen
    panel = openproblems._panel(dataset, sampling)
    # Make the refusal case explicit even when the shared fixture acquires
    # cross-batch donor references for other scientific tests.
    panel = replace(panel, donor=np.array([
        f"{donor}:{batch}" for donor, batch in zip(panel.donor, panel.batch, strict=True)
    ], dtype=object))
    assert all(len(b) == 1 for b in openproblems._donors_batches(panel).values())
    records = openproblems._horizontal_emd(panel, sampling, "r", "c")
    aggregate = _aggregate_record(records, "emd_mean_ct_horiz")
    assert aggregate.value is None
    assert aggregate.status is MetricStatus.INSUFFICIENT_DATA
    assert "spans two batches" in aggregate.reason


def _two_batches_two_samples_each(
    rng: np.random.Generator, between_batch: float, within_batch: float
) -> _Panel:
    """One donor, two batches, two samples each: the shape both EMD components need.

    ``between_batch`` shifts every cell of B2 relative to B1; ``within_batch`` shifts the
    second sample of each batch relative to the first, without touching the batches'
    relative position.
    """

    n = 140
    per_marker = np.concatenate([
        rng.normal(size=n),                                   # B1 / s1
        rng.normal(loc=within_batch, size=n),                # B1 / s2
        rng.normal(loc=between_batch, size=n),                # B2 / t1
        rng.normal(loc=between_batch + within_batch, size=n),  # B2 / t2
    ])
    return _Panel(
        expression=per_marker[:, None], markers=["CD3"],
        population=np.array(["T"] * (4 * n), dtype=object),
        sample=np.array(["s1"] * n + ["s2"] * n + ["t1"] * n + ["t2"] * n, dtype=object),
        donor=np.array(["D1"] * (4 * n), dtype=object),
        batch=np.array(["B1"] * (2 * n) + ["B2"] * (2 * n), dtype=object),
        rows=np.arange(4 * n),
    )


def test_the_horizontal_component_responds_to_a_removed_batch_offset() -> None:
    """Otherwise the component is a constant and its presence is a decoration."""

    manifest = _sampling_manifest_stub()
    panel = _two_batches_two_samples_each(np.random.default_rng(22), 3.0, 0.0)
    aligned = _with_expression(panel, panel.expression.copy())
    aligned.expression[2 * panel.n_cells() // 4:, 0] -= 3.0

    before = _aggregate_value(
        openproblems._horizontal_emd(panel, manifest, "r", "c"), "emd_mean_ct_horiz"
    )
    after = _aggregate_value(
        _emd_records(
            registered_metrics()[("compatibility", "emd_mean_ct_horiz")], manifest, "r", "c",
            openproblems._donor_batch_comparisons(aligned), aligned,
            EMD_MIN_CELLS_PER_SIDE, "donor",
        ),
        "emd_mean_ct_horiz",
    )
    assert before is not None and after is not None
    assert before > 2.0, "three bins of offset is at least three units of this scale"
    assert after < 0.5, "removing the offset should leave almost nothing to measure"


def test_the_vertical_component_sees_what_the_horizontal_one_cannot() -> None:
    """A correction that only aligns across batches leaves within-batch pairs alone.

    This is the reason the two components are kept apart, and the claim is in both
    specs' notes, so it is tested rather than asserted in prose.
    """

    manifest = _sampling_manifest_stub()
    panel = _two_batches_two_samples_each(np.random.default_rng(23), 3.0, 4.0)
    aligned = _with_expression(panel, panel.expression.copy())
    aligned.expression[2 * panel.n_cells() // 4:, 0] -= 3.0

    horizontal_before = _aggregate_value(
        openproblems._horizontal_emd(panel, manifest, "r", "c"), "emd_mean_ct_horiz"
    )
    horizontal_after = _aggregate_value(
        _emd_records(
            registered_metrics()[("compatibility", "emd_mean_ct_horiz")], manifest, "r", "c",
            openproblems._donor_batch_comparisons(aligned), aligned,
            EMD_MIN_CELLS_PER_SIDE, "donor",
        ),
        "emd_mean_ct_horiz",
    )
    vertical_before = _aggregate_value(
        openproblems._vertical_emd(panel, manifest, "r", "c"), "emd_mean_ct_vert"
    )
    vertical_after = _aggregate_value(
        _emd_records(
            registered_metrics()[("compatibility", "emd_mean_ct_vert")], manifest, "r", "c",
            openproblems._sample_pair_comparisons(aligned), aligned,
            EMD_MIN_CELLS_PER_SIDE, "batch",
        ),
        "emd_mean_ct_vert",
    )

    assert horizontal_after < horizontal_before
    assert vertical_before == pytest.approx(vertical_after, rel=0.5)
    assert vertical_after > 2.0, (
        "the two samples inside a batch still differ by four bins, and only the "
        "horizontal component was blind to that"
    )


def test_every_extracted_component_is_produced_for_the_frozen_partition(frozen) -> None:
    """All four rows come out of one call, each with an aggregate of its own."""

    _, dataset, design, sampling = frozen
    result = compute_openproblems_metrics(
        dataset, dataset.expression(), sampling, design, run_id="r", candidate_id="c",
    )
    aggregates = {
        r.metric_name for r in result.records if r.comparison_id is None
    }
    assert aggregates == {
        "average_batch_r2_ct", "cLisi", "emd_mean_ct_horiz", "emd_mean_ct_vert",
        "iLisi", "ratio_consistent_peaks",
    }
    assert result.status == "ok"
    assert not any(r.status is MetricStatus.FAILED for r in result.records)


def test_a_removed_batch_shift_moves_the_horizontal_component(frozen) -> None:
    """A real shift in one batch of one donor has to move the metric.

    Otherwise the component is a constant and its presence is a decoration.
    """

    _, dataset, design, sampling = frozen
    panel = openproblems._panel(dataset, sampling)
    spans = openproblems._donors_batches(panel)
    donors = [donor for donor, batches in spans.items() if len(batches) > 1]
    if not donors:
        pytest.skip("no donor in this design spans two batches")
    donor = donors[0]
    reference = spans[donor][0]

    before = _aggregate_value(
        openproblems._horizontal_emd(panel, sampling, "r", "c"), "emd_mean_ct_horiz"
    )
    assert before is not None
    shifted = dataset.expression().copy()
    rows = panel.rows[(panel.donor == donor) & (panel.batch == reference)]
    shifted[rows, 0] += 6.0
    moved = openproblems._candidate_panel(panel, dataset, shifted)
    after = _aggregate_value(
        _emd_records(
            registered_metrics()[("compatibility", "emd_mean_ct_horiz")], sampling, "r", "c",
            openproblems._donor_batch_comparisons(moved), moved,
            EMD_MIN_CELLS_PER_SIDE, "donor",
        ),
        "emd_mean_ct_horiz",
    )
    assert after > before


def test_a_constant_marker_is_undefined_in_the_batch_r2_not_zero() -> None:
    """The batch covariate cannot explain any of a zero total variance."""

    values = np.full(400, 2.0)
    panel = _panel(values[:, None], ["CD3"])
    records = openproblems._residual_batch_r2(
        panel, panel, _sampling_manifest_stub(), "r", "c"
    )
    assert _aggregate_value(records, "average_batch_r2_ct") is None
    assert all(record.value is None for record in records)
    assert any("constant" in (record.reason or "") for record in records)


def test_the_batch_r2_falls_when_the_batch_offset_is_removed() -> None:
    rng = np.random.default_rng(21)
    values = np.concatenate([rng.normal(size=200), rng.normal(loc=3.0, size=200)])
    panel = _panel(values[:, None], ["CD3"])
    aligned = values.copy()
    aligned[200:] -= 3.0
    manifest = _sampling_manifest_stub()

    assert _aggregate_value(
        openproblems._residual_batch_r2(panel, panel, manifest, "r", "c"),
        "average_batch_r2_ct",
    ) > 0.5
    assert _aggregate_value(
        openproblems._residual_batch_r2(
            panel, _with_expression(panel, aligned[:, None]), manifest, "r", "c"
        ),
        "average_batch_r2_ct",
    ) < 0.05


# --------------------------------------------------------------------------------------
# Parity with the pinned average_batch_r2/helper.py
# --------------------------------------------------------------------------------------


def _batch_r2_fixture(
    rng: np.random.Generator, batches: dict[str, np.ndarray], markers: list[str]
) -> tuple[ad.AnnData, ad.AnnData, list[str]]:
    """Two splits of one donor, the way the pinned component receives them.

    ``batches`` maps a batch label to that batch's cells.  Split 1 is the first batch and
    split 2 is everything else, which is how the component's script pairs them.
    """

    def _view(parts: dict[str, np.ndarray], tag: str) -> ad.AnnData:
        expression = np.concatenate(list(parts.values()))
        obs = pd.DataFrame(
            {
                "batch": [b for b, values in parts.items() for _ in range(len(values))],
                "donor": "D1",
                "cell_type": "T",
                "group": "case",
            },
            index=[f"{tag}{i}" for i in range(expression.shape[0])],
        )
        view = ad.AnnData(
            expression.astype("float64"),
            obs=obs.astype({"batch": "category", "cell_type": "category", "group": "category"}),
            var=pd.DataFrame(index=pd.Index(markers, dtype=object)),
        )
        view.layers["integrated"] = view.X.copy()
        return view

    first = list(batches)[0]
    return (
        _view({first: batches[first]}, "s1"),
        _view({b: v for b, v in batches.items() if b != first}, "s2"),
        markers,
    )


def _batch_r2_panel(parts: dict[str, np.ndarray], markers: list[str]) -> _Panel:
    expression = np.concatenate(list(parts.values()))
    n = expression.shape[0]
    sizes = [len(v) for v in parts.values()]
    return _Panel(
        expression=expression, markers=list(markers),
        population=np.array(["T"] * n, dtype=object),
        sample=np.array(["s1"] * n, dtype=object),
        donor=np.array(["D1"] * n, dtype=object),
        batch=np.concatenate([
            np.array([b] * size, dtype=object) for b, size in zip(parts, sizes, strict=True)
        ]),
        rows=np.arange(n),
    )


def test_the_pinned_batch_r2_helper_matches_the_extracted_component(
    upstream_dir: Path,
) -> None:
    """Both fits run on the same values and must agree, per marker."""

    pinned = _pinned_helper(upstream_dir, "average_batch_r2/helper.py")
    markers = ["CD3", "CD4", "CD8"]
    rng = np.random.default_rng(24)
    parts = {"B1": rng.normal(size=(120, 3)), "B2": rng.normal(loc=1.4, size=(150, 3))}
    split1, split2, _ = _batch_r2_fixture(rng, parts, markers)
    panel = _batch_r2_panel(parts, markers)

    want, want_names = pinned.batch_r2(split1, split2)
    records = openproblems._residual_batch_r2(
        panel, panel, _sampling_manifest_stub(), "r", "c"
    )
    got = {r.marker: r.value for r in records
           if r.metric_name == "average_batch_r2_ct" and r.marker is not None}
    assert want_names == markers
    for marker, expected in zip(markers, want, strict=True):
        assert got[marker] == pytest.approx(float(expected), rel=1e-9, abs=1e-12)
    assert _aggregate_value(records, "average_batch_r2_ct") == pytest.approx(
        float(np.mean(want)), rel=1e-9, abs=1e-12
    )


def test_the_batch_r2_fits_one_dummy_per_batch_and_not_a_single_contrast(
    upstream_dir: Path,
) -> None:
    """Three batches, so a one-column indicator would be a different model.

    ``pd.get_dummies(batch, drop_first=True)`` on three batches gives two columns and two
    degrees of freedom.  Collapsing them to one contrast throws away half the batch
    effect and lowers the R^2; the test pins the two apart so the fit cannot quietly
    regress to the wrong one.
    """

    pinned = _pinned_helper(upstream_dir, "average_batch_r2/helper.py")
    markers = ["CD3"]
    rng = np.random.default_rng(25)
    parts = {
        "B1": rng.normal(size=(120, 1)),
        "B2": rng.normal(loc=1.0, size=(110, 1)),
        "B3": rng.normal(loc=-1.0, size=(100, 1)),
    }
    split1, split2, _ = _batch_r2_fixture(rng, parts, markers)
    panel = _batch_r2_panel(parts, markers)
    want, _ = pinned.batch_r2(split1, split2)

    records = openproblems._residual_batch_r2(
        panel, panel, _sampling_manifest_stub(), "r", "c"
    )
    got = _aggregate_value(records, "average_batch_r2_ct")

    # the wrong model, on purpose: one contrast against B1
    values = panel.expression[:, 0]
    dummies = np.array([0.0 if b == "B1" else 1.0 for b in panel.batch])
    wrong = openproblems._empirical_r2(values, dummies)
    assert got == pytest.approx(float(want[0]), rel=1e-9, abs=1e-12)
    assert abs(got - wrong) > 0.05, (
        "if the single-contrast fit agreed with the pinned three-batch fit, this test "
        "would not be distinguishing the two models"
    )
    assert _batch_dummies_for(panel) == pytest.approx(np.column_stack([
        (panel.batch == "B2").astype(float), (panel.batch == "B3").astype(float)
    ]))


def _batch_dummies_for(panel: _Panel) -> np.ndarray:
    return openproblems._batch_dummies(panel.batch, sorted(set(panel.batch)))


# --------------------------------------------------------------------------------------
# Parity with the pinned lisi/helper.py
# --------------------------------------------------------------------------------------


def _lisi_fixture(
    rng: np.random.Generator, n: int = 400, n_markers: int = 6
) -> tuple[ad.AnnData, np.ndarray, np.ndarray, np.ndarray]:
    """One donor-free AnnData with the three obs columns the pinned helper reads.

    Every group is placed in both batches, so the pinned confounding check passes and the
    comparison is of the score rather than of the refusal.
    """

    batches = np.array(["B1", "B2"])[rng.integers(0, 2, size=n)]
    groups = np.array(["case", "control"])[rng.integers(0, 2, size=n)]
    populations = np.array(["T_CD4", "Mono", "NK"])[rng.integers(0, 3, size=n)]
    batch_effect = np.where(batches == "B2", 0.8, 0.0)
    expression = rng.normal(
        loc=np.where(populations[:, None] == "T_CD4", 2.5, 0.0) + batch_effect[:, None],
        scale=1.0, size=(n, n_markers),
    )
    view = ad.AnnData(
        expression.astype("float64"),
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(batches),
                "group": pd.Categorical(groups),
                "cell_type": pd.Categorical(populations),
            },
            index=[f"c{i}" for i in range(n)],
        ),
        var=pd.DataFrame(index=pd.Index(
            [f"M{i}" for i in range(n_markers)], dtype=object
        )),
    )
    view.layers["integrated"] = view.X.copy()
    return view, expression, batches, populations


def test_the_pinned_lisi_helpers_match_the_extracted_ones(upstream_dir: Path) -> None:
    """Per-cell LISI and both normalisations, on the same graph and the same seed."""

    pinned = _pinned_helper(upstream_dir, "lisi/helper.py")
    view, expression, batches, populations = _lisi_fixture(np.random.default_rng(26))
    assert not pinned._check_batch_group_confounding(view)

    want_ilisi, want_per_cell = pinned.compute_ilisi(view, int(len(set(batches))))
    got_per_cell = openproblems.local_inverse_simpson_index(
        expression, batches, seed=ILISI_SEED
    )
    assert np.allclose(got_per_cell, want_per_cell, rtol=0, atol=0)
    assert (float(np.nanmedian(want_per_cell)) - 1.0) / (
        len(set(batches)) - 1
    ) == pytest.approx(float(want_ilisi), rel=1e-12, abs=1e-15)

    want_clisi, want_clisi_per_cell = pinned.compute_clisi(view)
    got_clisi_per_cell = openproblems.local_inverse_simpson_index(
        expression, populations, seed=CLISI_SEED
    )
    assert np.allclose(got_clisi_per_cell, want_clisi_per_cell, rtol=0, atol=0)
    assert (len(set(populations)) - float(np.nanmedian(want_clisi_per_cell))) / (
        len(set(populations)) - 1
    ) == pytest.approx(float(want_clisi), rel=1e-12, abs=1e-15)


def test_the_confounding_check_fires_on_the_same_groups_as_the_pinned_one(
    upstream_dir: Path,
) -> None:
    """The pinned helper returns NaN for the whole metric; this refuses with a reason.

    The refusal is not a zero: a mixing score over a population that was only ever
    measured in one batch would describe the design, and reporting it as a number would
    let a candidate look well mixed because of how it was split.
    """

    pinned = _pinned_helper(upstream_dir, "lisi/helper.py")
    rng = np.random.default_rng(27)
    n = 400
    batches = np.array(["B1", "B2"])[rng.integers(0, 2, size=n)]
    groups = np.array(["case", "control"])[rng.integers(0, 2, size=n)]
    populations = np.array(["T_CD4", "Mono"])[rng.integers(0, 2, size=n)]
    # 20 cells of one group, and of one population, confined to a single batch.  The
    # pinned helper guards on group versus batch and the local one on population versus
    # batch, so the same construction confounds both and the two can be compared.
    groups[:20] = "confounded"
    batches[:20] = "B1"
    # every Mono cell is in B1, so the population is confounded with the batch; a few
    # cells would not do it, the other Mono cells would still reach B2
    batches[populations == "Mono"] = "B1"
    expression = rng.normal(size=(n, 6))
    view = ad.AnnData(
        expression.astype("float64"),
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(batches), "group": pd.Categorical(groups),
                "cell_type": pd.Categorical(populations),
            },
            index=[f"c{i}" for i in range(n)],
        ),
        var=pd.DataFrame(index=pd.Index([f"M{i}" for i in range(6)], dtype=object)),
    )
    view.layers["integrated"] = view.X.copy()

    assert pinned._check_batch_group_confounding(view) is True
    want_ilisi, want_per_cell = pinned.compute_ilisi(view, 2)
    assert np.isnan(want_ilisi) and want_per_cell.size == 0

    assert openproblems.confounded_groups(populations, batches) == ["Mono"]
    panel = _Panel(
        expression=expression, markers=[f"M{i}" for i in range(6)],
        population=populations.astype(object), sample=np.array(["s1"] * n, dtype=object),
        donor=np.array(["D1"] * n, dtype=object), batch=batches.astype(object),
        rows=np.arange(n),
    )
    records = openproblems._ilisi(panel, panel, _sampling_manifest_stub(), "r", "c")
    record = records[0]
    assert record.status is MetricStatus.UNDEFINED
    assert record.value is None
    assert "Mono" in record.reason and "confounded" in record.reason


def test_i_lisi_rises_when_the_batches_are_mixed_and_c_lisi_falls() -> None:
    """The two components must move in opposite directions, or one of them is broken.

    A correction that aligns batches also tends to blur populations, so a candidate that
    scores well on both has done something, and a candidate that scores well on one while
    the other collapses has destroyed structure.  That trade-off is only visible if both
    are computed on the same cells, which is why they are one function call.
    """

    manifest = _sampling_manifest_stub()
    rng = np.random.default_rng(28)
    n = 600
    batches = np.array(["B1", "B2"])[rng.integers(0, 2, size=n)]
    populations = np.array(["T_CD4", "Mono", "NK"])[rng.integers(0, 3, size=n)]
    # batch carries no signal at all, so the neighbourhood structure is free to mix
    well_mixed = rng.normal(loc=1.5, size=(n, 4))
    scrambled = rng.normal(size=(n, 4))
    panel = _Panel(
        expression=scrambled, markers=[f"M{i}" for i in range(4)],
        population=populations.astype(object),
        sample=np.array([f"s{i % 3}" for i in range(n)], dtype=object),
        donor=np.array(["D1"] * n, dtype=object), batch=batches.astype(object),
        rows=np.arange(n),
    )

    def _scores(values: np.ndarray) -> tuple[float, float]:
        shaped = _with_expression(panel, values)
        ilisi = openproblems._ilisi(panel, shaped, manifest, "r", "c")[0]
        clisi = openproblems._clisi(panel, shaped, manifest, "r", "c")[0]
        return ilisi.value, clisi.value

    mixed_ilisi, mixed_clisi = _scores(well_mixed)
    blurred_ilisi, blurred_clisi = _scores(panel.expression)

    assert blurred_ilisi is not None and mixed_ilisi is not None
    assert mixed_ilisi > blurred_ilisi, "mixing the batches must raise iLisi"
    assert mixed_clisi is not None and blurred_clisi is not None
    assert mixed_clisi > blurred_clisi, "destroying the populations must lower cLisi"
    assert 0.0 <= mixed_ilisi <= 1.0 and 0.0 <= mixed_clisi <= 1.0


def test_a_single_batch_or_a_single_population_is_undefined_rather_than_zero() -> None:
    """Both normalisations divide by (n_levels - 1), which has no value at one level."""

    manifest = _sampling_manifest_stub()
    n = 200
    values = np.random.default_rng(29).normal(size=(n, 3))
    populations = np.array(["T"] * n, dtype=object)
    one_batch = _Panel(
        expression=values, markers=["A", "B", "C"], population=populations,
        sample=np.array([f"s{i % 2}" for i in range(n)], dtype=object),
        donor=np.array(["D1"] * n, dtype=object), batch=np.array(["B1"] * n, dtype=object),
        rows=np.arange(n),
    )
    one_population = _Panel(
        expression=values, markers=["A", "B", "C"], population=populations,
        sample=np.array([f"s{i % 2}" for i in range(n)], dtype=object),
        donor=np.array(["D1"] * n, dtype=object),
        batch=np.array(["B1", "B2"] * (n // 2), dtype=object), rows=np.arange(n),
    )
    for panel, metric in ((one_batch, "iLisi"), (one_population, "cLisi")):
        compute = openproblems._ilisi if metric == "iLisi" else openproblems._clisi
        record = compute(panel, panel, manifest, "r", "c")[0]
        assert record.metric_name == metric
        assert record.status is MetricStatus.UNDEFINED
        assert record.value is None
        assert record.reason


def test_too_few_cells_for_the_neighbourhood_is_insufficient_data() -> None:
    """The pinned components never say what happens below 101 cells; this does."""

    n = LISI_NEIGHBOURS
    panel = _Panel(
        expression=np.random.default_rng(30).normal(size=(n, 3)),
        markers=["A", "B", "C"], population=np.array(["T"] * n, dtype=object),
        sample=np.array(["s1"] * n, dtype=object), donor=np.array(["D1"] * n, dtype=object),
        batch=np.array(["B1", "B2"] * (n // 2), dtype=object), rows=np.arange(n),
    )
    record = openproblems._ilisi(panel, panel, _sampling_manifest_stub(), "r", "c")[0]
    assert record.status is MetricStatus.INSUFFICIENT_DATA
    assert record.value is None
    assert str(LISI_NEIGHBOURS) in record.reason
