"""Multi-seed checks for the claims the package makes about determinism and stability.

Before this file the whole scientific suite ran at exactly one seed, 473, and that seed
was repeated as a literal in eight places.  A single seed cannot distinguish two
different kinds of claim:

* **Determinism** — the same seed must reproduce the same run.  This is checkable at one
  seed, and a bug that broke it would be caught there.
* **Stability** — a claim that holds for the fixture the maintainer happened to look at.
  This is *only* checkable across seeds.  One seed cannot tell a property of the method
  from a property of the draw.

So the sweep is over stability, and over exactly three things:

1. the seed is threaded — the configuration a run reads declares the seed its data was
   generated with (:func:`test_the_emitted_configuration_records_the_seed_the_data_was_generated_with`),
   and the seed is part of the analytical fingerprint;
2. the draw is a function of the seed, and reproducible from it;
3. the two claims that would be worthless if they held only at 473 — a correct offset
   removal reduces alignment error on the markers that carry the injected offset, and a
   confounded design refuses to recommend anything.

**What this deliberately does not do.** It does not re-run
``test_stage7_metrics.py``/``test_stage8_preservation.py``/
``test_stage9_decisions.py`` at each seed.  That is not a cost saving — ``tests/scientific``
takes 1m47s, so five passes would be under ten minutes, and the whole suite is affordable.
The blocker is structural instead: those files take their data from the ``frozen`` fixture
in ``tests/conftest.py``, which calls ``materialise_example`` without a seed, so their
data is 473 regardless of what a parametrisation says.  Sweeping them needs a
seed-aware conftest, which is a change to a shared file and therefore a central change,
not a change this file can make for itself.  Until then the honest scope of the sweep is
the three claims above, and ``docs/validation_results.md`` says so.

**Cost.** One materialised 15-sample, 3,000-cell example per (example, seed) pair, built
once per session and reused by every test that needs it.  Measured on the environment of
record: 25s for the whole file, against 1m47s for the single-seed scientific suite it
sits beside.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import yaml

from batchdx.cyto.config import DEFAULT_SEED, parse_config
from batchdx.cyto.decisions import assess_candidates
from batchdx.cyto.design import audit_design
from batchdx.cyto.diagnostics.base import DetectorContext, eligibility_for
from batchdx.cyto.examples import materialise_example
from batchdx.cyto.io import load_dataset
from batchdx.cyto.metrics import compute_common_metrics, compute_reductions
from batchdx.cyto.sampling import build_sampling_manifest, resolve_labels
from batchdx.cyto.schemas.enums import MetricStatus

#: Five seeds, the first being the historical default so the shipped examples stay
#: covered.  The rest are deliberately far apart, because nearby integers share
#: generator state in ways that make a near-miss look like a pass.
SEEDS = (473, 7, 101, 20250930, 65537)

#: Cells per sample.  Small, because nothing here is a power calculation: the claims are
#: about agreement across seeds, and a smaller fixture buys more seeds for the same time.
CELLS_PER_SAMPLE = 200

#: ``paired_reference_spec`` puts a technical offset on exactly these two markers
#: (B2: CD4 +0.8, HLA-DR -0.5; B3: CD4 -0.6, HLA-DR +0.7).  A correct correction must
#: reduce alignment error on them.  The 0.5 ratio is the tolerance already declared in
#: ``docs/validation_results.md`` ("offset removal on affected markers"), restated here
#: rather than re-derived, so a change to one is visible in the other.
AFFECTED_MARKERS = ("CD4", "HLADR")
MAX_REMAINING_FRACTION = 0.5

#: The evaluation draw only consumes the seed when it actually subsamples.  The shipped
#: examples cap at 1200 cells per sample against 200 generated here, so without this
#: override the manifest would be the *whole* dataset at every seed and "different seeds
#: draw differently" would be vacuously true of a draw that never happened.
DRAW_CAP = 120


@dataclass(frozen=True)
class World:
    """One materialised example, loaded and frozen, at one seed."""

    seed: int
    root: Path
    config: object
    dataset: object
    design: object
    sampling: object

    def context(self) -> DetectorContext:
        return DetectorContext(
            run_id="t", design=self.design, sampling=self.sampling,
            out_dir=str(self.root), seed=self.config.project.seed,
        )


@pytest.fixture(scope="session")
def worlds(tmp_path_factory: pytest.TempPathFactory):
    """Build each (example, seed) world once per session and memoise it."""

    cache: dict[tuple[str, int, int | None], World] = {}

    def get(name: str, seed: int, cap: int | None = None) -> World:
        key = (name, seed, cap)
        if key not in cache:
            root = tmp_path_factory.mktemp(f"{name}-{seed}-{cap}")
            materialise_example(name, root, cells_per_sample=CELLS_PER_SAMPLE, seed=seed)
            raw = yaml.safe_load((root / "config" / "config.yaml").read_text(encoding="utf-8"))
            if cap is not None:
                raw["evaluation"]["sampling"]["max_cells_per_sample"] = cap
            config = parse_config(raw, source=str(root / "config" / "config.yaml"))
            dataset = load_dataset(config, root)
            design = audit_design(
                list(dataset.samples.values()),
                dataset_id=dataset.manifest.dataset_id,
                correction_factor=config.design.correction_batch,
                protected=config.design.protect,
                donor_factor=config.design.donor,
                minimum_batches=config.design.minimum_batches,
                require_shared_reference=config.design.require_shared_reference,
            )
            labels, tier, confidence = resolve_labels(dataset, config, root)
            sampling = build_sampling_manifest(dataset, config, design, labels, tier, confidence)
            cache[key] = World(seed, root, config, dataset, design, sampling)
        return cache[key]

    return get


def _median_emd(result, marker: str) -> float:
    values = [
        float(r.value) for r in result.records
        if r.metric_name == "emd_within_population" and r.marker == marker
        and r.status is MetricStatus.OK and r.value is not None
    ]
    assert values, f"no ok within-population EMD for {marker}"
    return float(np.median(values))


# --------------------------------------------------------------------------------------
# 1. The seed is threaded
# --------------------------------------------------------------------------------------


def test_the_seed_has_exactly_one_default() -> None:
    """One literal, one home.

    This is the guard on the defect the sweep was written to prevent: 473 was repeated as
    a literal in eight places, and nothing would have failed if one of them were changed
    on its own.  ``ProjectConfig.seed``, ``SimulationSpec.seed``, the three named
    scenarios, ``materialise_example`` and ``DetectorContext.seed`` must all resolve to
    the same value, and that value must be the one ``config.DEFAULT_SEED`` names.
    """
    from batchdx.cyto.config import ProjectConfig
    from batchdx.cyto.diagnostics.base import DetectorContext as DC
    from batchdx.cyto.examples import EXAMPLES, _config_dict
    from batchdx.cyto.simulate import (
        SimulationSpec,
        confounded_spec,
        overlapping_no_reference_spec,
        paired_reference_spec,
    )

    observed = {
        "ProjectConfig.seed": ProjectConfig(id="x").seed,
        "SimulationSpec.seed": SimulationSpec(samples=[]).seed,
        "paired_reference_spec": paired_reference_spec().seed,
        "confounded_spec": confounded_spec().seed,
        "overlapping_no_reference_spec": overlapping_no_reference_spec().seed,
        "materialise_example": _config_dict("paired_reference", EXAMPLES["paired_reference"])["project"]["seed"],
        "DetectorContext.seed": DC(run_id="r", design=None, sampling=None, out_dir="o").seed,
    }
    assert set(observed.values()) == {DEFAULT_SEED}, observed


def test_the_emitted_configuration_records_the_seed_the_data_was_generated_with(
    worlds,
) -> None:
    """The regression this whole item exists for.

    ``_config_dict`` used to write ``"seed": 473`` as a literal while
    ``materialise_example(seed=7)`` generated the data with seed 7.  Every stage reads
    ``config.project.seed``, so the shipped run configuration described a run that was not
    the run: right data, wrong recorded provenance, and a fingerprint that did not
    correspond to the fixture.  Both numbers must be the one seed.
    """
    for seed in SEEDS:
        world = worlds("paired_reference", seed)
        emitted = json.loads((world.root / "config" / "example_manifest.json").read_text("utf-8"))
        assert emitted["seed"] == seed
        assert world.config.project.seed == seed, (
            f"seed {seed}: the configuration declares {world.config.project.seed} but the "
            f"data was generated with {seed}"
        )
        assert world.sampling.max_cells_per_sample == world.config.evaluation.sampling.max_cells_per_sample


def test_the_seed_is_part_of_the_analytical_fingerprint(worlds) -> None:
    """Two runs differing only in the seed are different analyses.

    ``Config.analytical_fingerprint`` excludes ``report`` and nothing else, so the seed is
    already covered.  That is a load-bearing property and this asserts it rather than
    assuming it: if someone ever widens the exclusion, two runs that resampled
    differently would begin sharing a fingerprint and reusing each other's cached
    correction.
    """
    seen = {}
    for seed in SEEDS:
        world = worlds("paired_reference", seed)
        seen[seed] = world.config.analytical_fingerprint()
    assert len(set(seen.values())) == len(SEEDS), (
        f"seeds collapsed onto {len(set(seen.values()))} fingerprints: {seen}"
    )


# --------------------------------------------------------------------------------------
# 2. The draw is a function of the seed, and reproducible from it
# --------------------------------------------------------------------------------------


def test_the_same_seed_reproduces_the_data_and_the_frozen_draw(worlds) -> None:
    """Determinism, which one seed could already have checked.

    Asserted as exact equality, not a tolerance: both the generator and the manifest
    builder take an explicit ``default_rng``, so anything short of bit-identity is a
    defect rather than sampling noise.
    """
    for seed in SEEDS:
        first = worlds("paired_reference", seed, cap=DRAW_CAP)
        second = worlds("paired_reference", seed, cap=DRAW_CAP)
        assert first is second, "the session cache failed to memoise; cost, not correctness"
        # A second, independent build at the same seed must agree with the cached one.
        rebuilt_root = first.root.parent / f"rebuild-{seed}"
        materialise_example("paired_reference", rebuilt_root,
                            cells_per_sample=CELLS_PER_SAMPLE, seed=seed)
        rebuilt_raw = yaml.safe_load(
            (rebuilt_root / "config" / "config.yaml").read_text(encoding="utf-8")
        )
        rebuilt_raw["evaluation"]["sampling"]["max_cells_per_sample"] = DRAW_CAP
        rebuilt = parse_config(rebuilt_raw)
        rebuilt_dataset = load_dataset(rebuilt, rebuilt_root)
        rebuilt_design = audit_design(
            list(rebuilt_dataset.samples.values()),
            dataset_id=rebuilt_dataset.manifest.dataset_id,
            correction_factor=rebuilt.design.correction_batch,
            protected=rebuilt.design.protect,
            donor_factor=rebuilt.design.donor,
            minimum_batches=rebuilt.design.minimum_batches,
            require_shared_reference=rebuilt.design.require_shared_reference,
        )
        labels, tier, confidence = resolve_labels(rebuilt_dataset, rebuilt, rebuilt_root)
        rebuilt_sampling = build_sampling_manifest(
            rebuilt_dataset, rebuilt, rebuilt_design, labels, tier, confidence
        )
        assert np.array_equal(
            first.dataset.expression(), rebuilt_dataset.expression()
        ), f"seed {seed}: the generator is not reproducible"
        assert first.sampling.labels_by_cell == rebuilt_sampling.labels_by_cell, (
            f"seed {seed}: the evaluation draw is not reproducible from the seed"
        )
        assert first.sampling.per_sample_counts == rebuilt_sampling.per_sample_counts


def test_different_seeds_give_different_data_and_a_different_draw(worlds) -> None:
    """The sweep is not vacuous.

    If every seed produced the same data and the same draw, the other four seeds would
    prove nothing while looking like four passing tests.  This is the assertion that makes
    the sweep mean anything, and it is why the draw is forced through a cap.
    """
    draws = {seed: worlds("paired_reference", seed, cap=DRAW_CAP) for seed in SEEDS}
    expressions = {s: w.dataset.expression() for s, w in draws.items()}
    reference = SEEDS[0]
    for other in SEEDS[1:]:
        assert not np.array_equal(expressions[reference], expressions[other]), (
            f"seeds {reference} and {other} produced identical data; the generator ignores "
            f"the seed, so the sweep is testing nothing"
        )
        symmetric_difference = (
            set(draws[reference].sampling.labels_by_cell)
            ^ set(draws[other].sampling.labels_by_cell)
        )
        assert symmetric_difference, (
            f"seeds {reference} and {other} produced an identical evaluation draw despite "
            f"the cap; the draw is not a function of the seed"
        )
        # Each sample must still be capped, or the draw above was a no-op.
        assert set(draws[other].sampling.per_sample_counts.values()) == {DRAW_CAP}


# --------------------------------------------------------------------------------------
# 3. Stability of the claims that would be worthless at one seed
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("seed", SEEDS)
def test_the_correct_correction_reduces_error_on_the_markers_that_carry_the_offset(
    worlds, tmp_path: Path, seed: int
) -> None:
    """The headline scientific claim, re-checked per seed.

    ``fixture_offset_removal`` is the *correct* answer for this fixture: it subtracts each
    batch's offset as measured on its reference aliquot.  It must therefore leave less
    alignment error than the uncorrected baseline on exactly the markers
    ``paired_reference_spec`` shifted.  The tolerance is the one already declared in
    ``docs/validation_results.md``; the worst of the five seeds measured here leaves
    0.31 of the original error, so the bound is not close to being the thing being tested.
    """
    from batchdx.cyto.adapters import get_backend

    world = worlds("paired_reference", seed)
    context = world.context()
    baseline_matrix, _ = get_backend("none").run(
        world.dataset, "uncorrected", {}, context
    )
    corrected_matrix, _ = get_backend("fixture_offset_removal").run(
        world.dataset, "offset", {}, context
    )
    baseline = compute_common_metrics(
        world.dataset, baseline_matrix, world.sampling, world.design,
        run_id="t", candidate_id="uncorrected",
    )
    corrected = compute_common_metrics(
        world.dataset, corrected_matrix, world.sampling, world.design,
        run_id="t", candidate_id="offset",
    )
    for marker in AFFECTED_MARKERS:
        before = _median_emd(baseline, marker)
        after = _median_emd(corrected, marker)
        assert after < MAX_REMAINING_FRACTION * before, (
            f"seed {seed}: {marker} EMD {before:.4f} -> {after:.4f}, which is "
            f"{after / before:.3f} of the uncorrected baseline"
        )


@pytest.mark.parametrize("seed", SEEDS)
def test_a_metric_on_the_same_data_is_identical_across_repeated_runs(
    worlds, seed: int
) -> None:
    """A metric is a function of its inputs, not of when it was called.

    Exact equality on every record's value and status.  Anything less would mean a metric
    is reading uninitialised memory, a shared RNG, or iteration over a set — the class of
    defect that a single run hides and a second run exposes.
    """
    from batchdx.cyto.adapters import get_backend

    world = worlds("paired_reference", seed)
    context = world.context()
    matrix, _ = get_backend("none").run(world.dataset, "uncorrected", {}, context)
    first = compute_common_metrics(
        world.dataset, matrix, world.sampling, world.design,
        run_id="t", candidate_id="uncorrected",
    )
    second = compute_common_metrics(
        world.dataset, matrix, world.sampling, world.design,
        run_id="t", candidate_id="uncorrected",
    )
    assert len(first.records) == len(second.records) > 0
    for a, b in zip(first.records, second.records, strict=True):
        assert a.metric_name == b.metric_name
        assert a.status is b.status
        assert a.value == b.value, f"seed {seed}: {a.metric_name} changed between runs"


@pytest.mark.parametrize("seed", SEEDS)
def test_a_confounded_design_refuses_to_recommend_a_correction(worlds, seed: int) -> None:
    """Abstention is the claim most likely to be a property of one fixture.

    ``confounded_spec`` puts every case sample in B1 and every control in B2, so batch and
    group are the same variable.  There is no analysis that can separate them, and the
    correct behaviour is to say so rather than to report a correction.  Measured across
    these five seeds the refusal is identical every time, which is the point: a
    non-identifiable design must not become identifiable on an unlucky draw.
    """
    from batchdx.cyto.adapters import get_backend

    world = worlds("confounded", seed)
    assert world.design.confounded_factors == ["group"]
    assert world.design.below_design_floor
    supports = {(c.batch_a, c.batch_b): c.support for c in world.design.comparisons}
    assert supports, "the confounded design reported no comparisons at all"
    for (batch_a, batch_b), support in supports.items():
        assert support.value == "non_identifiable", (
            f"seed {seed}: {batch_a} vs {batch_b} came out {support.value!r} on a design "
            f"where batch and group are the same variable"
        )
    cycombine = eligibility_for(world.design, "cycombine")
    assert not cycombine.eligible
    assert cycombine.reason_code == "batch_not_identifiable"
    assert not cycombine.override_allowed, (
        "a non-identifiable contrast must not be waveable through by an exploratory "
        "override, whatever the seed"
    )

    context = world.context()
    matrix, correction = get_backend("none").run(
        world.dataset, "uncorrected", {}, context
    )
    metrics = compute_common_metrics(
        world.dataset, matrix, world.sampling, world.design,
        run_id="t", candidate_id="uncorrected",
    )
    reductions = compute_reductions(
        {"uncorrected": metrics}, world.sampling, run_id="t",
        identity_candidate_id="uncorrected",
    )
    decision = assess_candidates(
        world.design, world.sampling, {"uncorrected": correction},
        {"uncorrected": metrics}, reductions,
        run_id="t", policy_name="conservative", identity_candidate_id="uncorrected",
    )
    assert decision.recommended_candidate_id is None, (
        f"seed {seed}: a confounded design recommended {decision.recommended_candidate_id!r}"
    )


@pytest.mark.parametrize("seed", SEEDS)
def test_a_calibrated_policy_can_still_recommend_at_every_seed(worlds, seed: int) -> None:
    """A gate that never opens is not a gate.

    The companion to the abstention test: the refusal above must come from the design
    being unsupportable, not from the policy being permanently closed.  A policy whose
    thresholds declare a real calibration is exercised at every seed and must reach
    ``recommended_under_policy`` on the reference-supported design — the same asymmetry
    check, in the other direction.
    """
    from batchdx.cyto.adapters import get_backend
    from batchdx.cyto.decisions import POLICIES, Policy

    world = worlds("paired_reference", seed)
    context = world.context()
    baseline_matrix, baseline_correction = get_backend("none").run(
        world.dataset, "uncorrected", {}, context
    )
    corrected_matrix, corrected = get_backend("fixture_offset_removal").run(
        world.dataset, "offset", {}, context
    )
    baseline = compute_common_metrics(
        world.dataset, baseline_matrix, world.sampling, world.design,
        run_id="t", candidate_id="uncorrected",
    )
    corrected_metrics = compute_common_metrics(
        world.dataset, corrected_matrix, world.sampling, world.design,
        run_id="t", candidate_id="offset",
    )
    # The reductions are computed over *both* candidates: the reduction is measured
    # against the identity candidate, so passing only the corrected one would compare a
    # candidate against nothing.
    reductions = compute_reductions(
        {"uncorrected": baseline, "offset": corrected_metrics},
        world.sampling, run_id="t", identity_candidate_id="uncorrected",
    )
    # Registered rather than passed inline: ``assess_candidates`` takes a policy *name*,
    # and the uncalibrated-origin wording is what distinguishes this from the default.
    name = f"calibrated_multiseed_{seed}"
    POLICIES[name] = Policy(
        name=name,
        threshold_origins={
            k: "calibrated against the synthetic paired-reference fixture"
            for k in Policy(name="x").threshold_origins
        },
    )
    try:
        decision = assess_candidates(
            world.design, world.sampling,
            {"uncorrected": baseline_correction, "offset": corrected},
            {"uncorrected": baseline, "offset": corrected_metrics},
            reductions, run_id="t", policy_name=name,
            identity_candidate_id="uncorrected",
        )
    finally:
        POLICIES.pop(name, None)
    states = {a.candidate_id: a.state.value for a in decision.assessments}
    assert states["offset"] == "recommended_under_policy", (
        f"seed {seed}: a calibrated policy over a correct correction on a "
        f"reference-supported design gave {states['offset']!r}; if the policy is the "
        f"problem, this test and the abstention test above are not measuring different "
        f"things"
    )
    assert decision.recommended_candidate_id == "offset"
