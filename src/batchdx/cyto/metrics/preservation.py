"""Prespecified biological contrasts, sample-level inference, and uncertainty (Stage 8).

The scientific rule this module exists to enforce
---------------------------------------------------
**Cells are not replicates.**  A tube of ten thousand events is one observation.  Every
number emitted here is computed from *per-sample means* of the frozen evaluation cells,
aggregated to donors, and tested with a model whose replicate unit is the donor.  A test
that ran on cells would report an effect size orders of magnitude too large and a
p-value orders of magnitude too small, and it would do so silently, because nothing in a
p-value tells you which unit fed it.

Uncertainty
-----------
Every record that reports a value reports an interval for it, and says how the interval
was obtained:

* ``donor_block_bootstrap_percentile_95`` — donors (or samples, when no donor factor is
  declared) resampled with replacement; the interval is the percentile interval of the
  statistic itself.
* ``restricted_permutation_null_95`` — the group label is permuted **inside its own
  block**, so the null keeps the design's structure; the interval is the central 95 % of
  the null distribution of the statistic.
* ``permutation_p_monte_carlo_95`` — the Clopper-Pearson interval on the fraction of
  resamples at least as extreme, which is the only uncertainty a Monte Carlo p has.

Refusal, not silence
--------------------
A contrast the design cannot support is refused with a reason.  It is never skipped,
never scored zero, and never replaced by a different algorithm.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator
from scipy import stats

from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.metrics.records import metric_record
from batchdx.cyto.metrics.registry import (
    BIOLOGICAL_CONTRAST_ADMISSIBILITY,
    BIOLOGICAL_CONTRAST_COHENS_D,
    BIOLOGICAL_CONTRAST_EFFECT_SIZE,
    BIOLOGICAL_CONTRAST_PERMUTATION_P,
    BIOLOGICAL_CONTRAST_RANK_STATISTIC,
    BIOLOGICAL_CONTRAST_TTEST_STATISTIC,
    MetricSpec,
)
from batchdx.cyto.schemas.contracts import (
    DesignAudit,
    MetricRecord,
    MetricResult,
    SampleRecord,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import MetricStatus

#: A (population, sample) pair below this many frozen evaluation cells contributes no
#: per-sample mean.  The floor the OpenProblems functional-marker component uses.
MIN_CELLS_PER_SAMPLE_MEAN = 5
#: Minimum replicates per level for a rank sum, and the floor below which a record
#: reports ``insufficient_data`` rather than a number.
MIN_SAMPLES_PER_LEVEL = 2
#: A donor has to contribute both levels before the contrast is treated as paired.
MIN_PAIRED_DONORS = 3
#: Below this many reachable label assignments the restricted permutation has almost no
#: null distribution, so it is reported as insufficient data rather than as a p-value.
MIN_PERMUTATION_ASSIGNMENTS = 20
#: Resamples per interval.  Reported in every record's reason, together with the seed.
DEFAULT_N_RESAMPLES = 2000
#: The smallest denominator a ratio is formed with, in units of the declared expression
#: scale.  Below this a standard deviation is not a measurement but the residue of
#: floating-point summation, and a ratio built on it is a number with no meaning.  This is
#: the same absolute floor the dispersion metric and the relative reduction use, on the
#: same scale, and it is applied rather than skipped: a denominator below it yields
#: ``undefined`` instead of a plausible-looking quotient.
ROUND_OFF_FLOOR = 1e-9

#: The three interval kinds, in a closed vocabulary.  The percentile level is appended.
BOOTSTRAP_KIND = "donor_block_bootstrap_percentile"
PERMUTATION_NULL_KIND = "restricted_permutation_null"
PERMUTATION_P_KIND = "permutation_p_monte_carlo"

#: The records emitted per (contrast, population, marker, model), in a fixed order.
CONTRASTAT_RECORDS: tuple[MetricSpec, ...] = (
    BIOLOGICAL_CONTRAST_EFFECT_SIZE,
    BIOLOGICAL_CONTRAST_COHENS_D,
    BIOLOGICAL_CONTRAST_TTEST_STATISTIC,
    BIOLOGICAL_CONTRAST_RANK_STATISTIC,
    BIOLOGICAL_CONTRAST_PERMUTATION_P,
)


class ContrastSpec(BaseModel):
    """One prespecified biological contrast.

    A contrast is written down *before* the data are looked at: which factor, which two
    levels, in which order, optionally which way the effect is expected to point, and
    which covariates the claim is conditional on.  Whatever the design cannot support is
    refused by :func:`resolve_contrast` rather than quietly dropped.

    ``level_a`` is always the alphabetically first of the two levels.  The normalisation
    is deliberate and matches the OpenProblems functional-marker component: with a
    sorted order a positive effect size means the same biological thing in every run and
    a sign cannot flip because two samples happened to be listed the other way round.
    ``expected_direction`` is stated in that normalised order.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    contrast_id: str
    factor: str
    level_a: str
    level_b: str
    expected_direction: Literal["none", "increase_in_a", "increase_in_b"] = "none"
    #: Covariates the claim is conditional on.  One nested inside a contrast level
    #: confounds the contrast with the covariate and refuses the contrast.
    protected_covariates: tuple[str, ...] = ()
    #: Restrict to these markers, usually the panel's functional markers.  ``None`` means
    #: every evaluation marker on the frozen partition.
    markers: tuple[str, ...] | None = None
    #: Restrict to these frozen populations.  ``None`` means every frozen population.
    populations: tuple[str, ...] | None = None
    #: ``auto`` matches the samples inside a donor when the donor factor allows it.
    model: Literal["auto", "unpaired", "paired"] = "auto"
    min_samples_per_level: int = MIN_SAMPLES_PER_LEVEL
    n_resamples: int = DEFAULT_N_RESAMPLES
    alpha: float = 0.05
    seed: int = 0

    @model_validator(mode="before")
    @classmethod
    def _normalise_level_order(cls, data):
        if isinstance(data, dict) and data.get("level_a") and data.get("level_b"):
            first, second = sorted((str(data["level_a"]), str(data["level_b"])))
            data = {**data, "level_a": first, "level_b": second}
        return data

    @model_validator(mode="after")
    def _check_declared(self):
        if self.level_a == self.level_b:
            raise ValueError(
                f"contrast {self.contrast_id!r}: a contrast needs two different levels, "
                f"got {self.level_a!r} twice"
            )
        if self.n_resamples < 1:
            raise ValueError("n_resamples must be at least 1")
        if not 0.0 < self.alpha < 1.0:
            raise ValueError("alpha must lie strictly between 0 and 1")
        return self

    @property
    def interval_level(self) -> int:
        """The percentile level of every interval this contrast reports, as a percentage."""

        return int(round(100.0 * (1.0 - self.alpha)))


@dataclass(frozen=True)
class ContrastPlan:
    """A contrast that survived validation, with its replicate unit already decided."""

    spec: ContrastSpec
    level_a: str
    level_b: str
    model: str  # "paired" or "unpaired"
    block_unit: str  # "donor_id" or "sample_id"
    rows: tuple[str, ...]  # replicate-unit keys, in the order they enter every matrix
    blocks: tuple[tuple[str, tuple[int, ...]], ...]  # permutation blocks -> row indices
    reachable_permutations: int
    contributing: tuple[str, ...]
    block_of: dict[str, str]
    level_of: dict[str, str]
    markers: tuple[str, ...]
    populations: tuple[str, ...]

    @property
    def n_samples(self) -> int:
        return len(self.contributing)

    @property
    def n_donors(self) -> int:
        return len(self.rows)

    @property
    def comparison_id(self) -> str:
        return (
            f"contrast:{self.spec.contrast_id}:{self.spec.factor}:"
            f"{self.level_a}__vs__{self.level_b}:{self.model}"
        )


@dataclass(frozen=True)
class ContrastRefusal:
    """A contrast that may not be run, with the reason it may not."""

    spec: ContrastSpec
    reason_code: str
    reason: str

    @property
    def comparison_id(self) -> str:
        return (
            f"contrast:{self.spec.contrast_id}:{self.spec.factor}:"
            f"{self.spec.level_a}__vs__{self.spec.level_b}:refused"
        )


@dataclass(frozen=True)
class ContrastStatistics:
    """Every scalar one (contrast, population, marker) needs, each in its own units."""

    effect: float | None
    cohens_d: float | None
    t_statistic: float | None
    t_pvalue: float | None
    rank_statistic: float | None
    rank_pvalue: float | None
    n_level_a: int
    n_level_b: int
    note: str | None = None


# --------------------------------------------------------------------------------------
# Sample-table access
# --------------------------------------------------------------------------------------


def _factor_value(sample: SampleRecord, factor: str) -> str | None:
    """One factor's value for a sample, from the field or the technical-factor table."""

    if hasattr(sample, factor):
        value = getattr(sample, factor)
        return None if value is None else str(value)
    return sample.technical_factors.get(factor)


def _factor_values(dataset: Dataset, factor: str) -> dict[str, str | None]:
    return {sid: _factor_value(rec, factor) for sid, rec in dataset.samples.items()}


# --------------------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------------------


def resolve_contrast(
    dataset: Dataset,
    design: DesignAudit,
    sampling: SamplingManifest,
    spec: ContrastSpec,
) -> ContrastPlan | ContrastRefusal:
    """Decide whether a prespecified contrast may be run, and with which replicate unit.

    Every branch that refuses names the condition that failed.  A contrast is never
    dropped without a record, and never reported as a score of zero.
    """

    values = _factor_values(dataset, spec.factor)
    levels_present = sorted({v for v in values.values() if v is not None})
    if not levels_present:
        return ContrastRefusal(
            spec, "unknown_factor",
            f"no sample declares {spec.factor!r}, so the contrast has no data to run on; "
            f"the declared levels were {spec.level_a!r} and {spec.level_b!r}",
        )

    missing = [lvl for lvl in (spec.level_a, spec.level_b) if lvl not in levels_present]
    if missing:
        return ContrastRefusal(
            spec, "levels_absent",
            f"level(s) {missing} of contrast {spec.contrast_id!r} are not present in "
            f"{spec.factor!r}; the sample table carries {levels_present}. A contrast is "
            "prespecified, not discovered after the fact.",
        )

    # -- identifiability ---------------------------------------------------------------
    if spec.factor not in design.protected_factors:
        return ContrastRefusal(
            spec, "factor_not_protected",
            f"{spec.factor!r} is not among the design's declared protected factors "
            f"{list(design.protected_factors)}, so a difference in it is not separately "
            "estimable from the batch effect; declare it in design.protect first",
        )
    if spec.factor in design.confounded_factors:
        return ContrastRefusal(
            spec, "non_identifiable",
            f"the design audit reports {spec.factor!r} perfectly confounded with "
            f"{design.correction_factor!r} (confounded_factors="
            f"{list(design.confounded_factors)}); the batch term and the {spec.factor!r} "
            "term span the same subspace, so no covariate adjustment separates them and no "
            "preservation claim can rest on this contrast",
        )
    batch_values = _factor_values(dataset, design.correction_factor)
    level_by_batch: dict[str, set[str]] = {}
    for sid, level in values.items():
        if level in {spec.level_a, spec.level_b} and batch_values.get(sid):
            level_by_batch.setdefault(str(batch_values[sid]), set()).add(str(level))
    if len(level_by_batch) > 1 and all(len(v) == 1 for v in level_by_batch.values()):
        return ContrastRefusal(
            spec, "non_identifiable",
            f"the two levels of {spec.factor!r} occupy disjoint batches "
            f"({ {b: sorted(v)[0] for b, v in sorted(level_by_batch.items())} }); the "
            "biological contrast and the batch contrast are the same contrast",
        )

    # -- declared covariates -----------------------------------------------------------
    for covariate in spec.protected_covariates:
        covariate_values = _factor_values(dataset, covariate)
        present = {v for v in covariate_values.values() if v is not None}
        if not present:
            return ContrastRefusal(
                spec, "covariate_absent",
                f"protected covariate {covariate!r} is declared on this contrast but no "
                "sample carries it, so the contrast would not be conditional on anything",
            )
        if present <= {spec.level_a, spec.level_b}:
            return ContrastRefusal(
                spec, "covariate_nested_in_contrast",
                f"protected covariate {covariate!r} is nested inside {spec.factor!r}: its "
                f"levels {sorted(present)} are a subset of "
                f"{{{spec.level_a!r}, {spec.level_b!r}}}, so the contrast and the "
                "covariate cannot be separated",
            )

    # -- scope -------------------------------------------------------------------------
    evaluation_markers = tuple(sampling.evaluation_markers)
    markers = evaluation_markers if spec.markers is None else tuple(spec.markers)
    if not markers:
        return ContrastRefusal(
            spec, "no_evaluation_markers",
            "the frozen partition names no marker flagged use_for_evaluation, so there is "
            "nothing to test",
        )
    absent_markers = [m for m in markers if m not in evaluation_markers]
    if absent_markers:
        return ContrastRefusal(
            spec, "marker_absent",
            f"marker(s) {absent_markers} are not flagged use_for_evaluation on the frozen "
            f"partition, which carries {list(evaluation_markers)}",
        )
    frozen_populations = sorted(set(sampling.labels_by_cell.values()))
    populations = (
        tuple(frozen_populations) if spec.populations is None else tuple(spec.populations)
    )
    absent_populations = [p for p in populations if p not in frozen_populations]
    if not populations or absent_populations:
        return ContrastRefusal(
            spec, "population_absent",
            f"population(s) {absent_populations or list(populations)} are not on the frozen "
            f"partition, which carries {frozen_populations}",
        )

    # -- replicate unit ----------------------------------------------------------------
    contributing = tuple(sorted(
        sid for sid, level in values.items() if level in {spec.level_a, spec.level_b}
    ))
    donors = {sid: rec.donor_id for sid, rec in dataset.samples.items()}
    has_donor = all(donors.get(sid) for sid in contributing)
    # The replicate unit is the donor when every contributing sample declares one, and
    # the sample itself otherwise.  Mixing the two would make a resampled block
    # sometimes a donor and sometimes a sample.
    block_of = {sid: (str(donors[sid]) if has_donor else sid) for sid in contributing}
    block_unit = "donor_id" if has_donor else "sample_id"
    level_of = {sid: str(values[sid]) for sid in contributing}

    levels_in_block: dict[str, set[str]] = {}
    for sid in contributing:
        levels_in_block.setdefault(block_of[sid], set()).add(level_of[sid])
    complete = sorted(
        block for block, levels in levels_in_block.items()
        if levels == {spec.level_a, spec.level_b}
    )
    paired = has_donor and len(complete) >= MIN_PAIRED_DONORS
    if spec.model == "paired" and not paired:
        return ContrastRefusal(
            spec, "pairing_unsupported",
            f"the contrast declares a paired model but only {len(complete)} donor(s) "
            f"contribute both {spec.level_a!r} and {spec.level_b!r}; at least "
            f"{MIN_PAIRED_DONORS} complete donor pairs are required, and a donor factor is "
            + ("declared for every contributing sample" if has_donor
               else "not declared for every contributing sample"),
        )
    model = "paired" if paired else "unpaired"
    rows = tuple(complete) if paired else tuple(sorted(levels_in_block))
    index = {row: i for i, row in enumerate(rows)}

    if paired:
        # Each block is one donor holding both levels, so relabelling inside it is a swap
        # of that donor's two values, and every donor has exactly two such assignments.
        blocks = tuple((f"donor:{row}", (index[row],)) for row in rows)
        reachable = 2 ** len(rows)
    else:
        # Each block is a batch, and relabelling inside it permutes the labels among that
        # batch's own samples, so the number of each level in the batch cannot change.
        grouped: dict[str, list[int]] = {}
        for sid in contributing:
            grouped.setdefault(str(batch_values.get(sid)), []).append(index[block_of[sid]])
        blocks = tuple((f"batch:{b}", tuple(sorted(v))) for b, v in sorted(grouped.items()))
        holds_a = {
            row: any(level_of[s] == spec.level_a for s in contributing if block_of[s] == row)
            for row in rows
        }
        reachable = 1
        for _, indices in blocks:
            reachable *= math.comb(len(indices), sum(1 for i in indices if holds_a[rows[i]]))
            reachable = min(reachable, 10**9)

    return ContrastPlan(
        spec=spec,
        level_a=spec.level_a,
        level_b=spec.level_b,
        model=model,
        block_unit=block_unit,
        rows=rows,
        blocks=blocks,
        reachable_permutations=min(reachable, 10**9),
        contributing=contributing,
        block_of=block_of,
        level_of=level_of,
        markers=markers,
        populations=populations,
    )


# --------------------------------------------------------------------------------------
# Donor-level summaries
# --------------------------------------------------------------------------------------


def _frozen_rows_by_sample(
    dataset: Dataset, sampling: SamplingManifest
) -> dict[str, dict[str, list[int]]]:
    """``population -> sample_id -> dataset row indices`` over the frozen evaluation cells."""

    obs = dataset.obs
    obs = obs.loc[obs[CELL_ID].isin(sampling.labels_by_cell.keys())].copy()
    obs["_label"] = obs[CELL_ID].map(sampling.labels_by_cell).astype(str)
    positions = {cid: i for i, cid in enumerate(dataset.obs[CELL_ID])}
    obs["_row"] = obs[CELL_ID].map(positions).to_numpy()

    out: dict[str, dict[str, list[int]]] = {}
    for (population, sample_id), block in obs.groupby(["_label", SAMPLE_ID], sort=True):
        out.setdefault(str(population), {})[str(sample_id)] = sorted(
            int(r) for r in block["_row"].to_numpy()
        )
    return out


def donor_level_matrix(
    plan: ContrastPlan, samples: list[str], values: np.ndarray
) -> np.ndarray:
    """``(n_replicates, 2)`` matrix of donor-level summaries; column 0 is level_a.

    ``values`` holds one per-sample mean per entry of ``samples``.  A donor with several
    samples in one level has them averaged here, before anything else happens, so a donor
    with three timepoints cannot outvote a donor with one.  A donor that contributes no
    sample in a level is ``nan`` there, never a fabricated zero.
    """

    matrix = np.full((len(plan.rows), 2), np.nan, dtype=float)
    for column, level in enumerate((plan.level_a, plan.level_b)):
        for i, row in enumerate(plan.rows):
            selected = [
                k for k, sid in enumerate(samples)
                if plan.block_of.get(sid) == row and plan.level_of.get(sid) == level
            ]
            if selected:
                matrix[i, column] = float(values[selected].mean())
    return matrix


# --------------------------------------------------------------------------------------
# Statistics.  Each returns None for a quantity that is undefined, never zero.
# --------------------------------------------------------------------------------------


def _safe(func, *args, **kwargs):
    """Run a scipy test, returning its p-value, or ``None`` when it is undefined.

    Degenerate inputs are probed on purpose — a resample can collapse a variance to zero
    — and translated into ``None`` here, so scipy's warning about the same thing is
    silenced rather than printed once per resample.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = func(*args, **kwargs)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            return None
    pvalue = getattr(result, "pvalue", None)
    if pvalue is None or not np.isfinite(pvalue):
        return None
    return float(pvalue)


def welch_p(a: np.ndarray, b: np.ndarray) -> float | None:
    """The two-sided p-value of Welch's unequal-variance t-test on replicate-level values."""

    return _safe(stats.ttest_ind, a, b, equal_var=False)


def paired_p(a: np.ndarray, b: np.ndarray) -> float | None:
    """The two-sided p-value of the paired t-test on matched replicate-level values."""

    return _safe(stats.ttest_rel, a, b)


def rank_sum_p(a: np.ndarray, b: np.ndarray) -> float | None:
    """The two-sided p-value of the Wilcoxon rank-sum test (Mann-Whitney U)."""

    return _safe(stats.mannwhitneyu, a, b, alternative="two-sided")


def signed_rank_p(a: np.ndarray, b: np.ndarray) -> float | None:
    """The two-sided p-value of the Wilcoxon signed-rank test on matched values."""

    return _safe(stats.wilcoxon, a - b, alternative="two-sided")


def _entry_ranks(matrix: np.ndarray, model: str) -> np.ndarray:
    """Average ranks of every finite entry, computed once and reused by every resample.

    Relabelling moves an entry from one column to the other; it never changes the
    multiset of values.  The ranks can therefore be computed once, which is what makes
    the rank statistics affordable over thousands of resamples.
    """

    if model == "paired":
        differences = matrix[:, 0] - matrix[:, 1]
        ranks = np.full(matrix.shape[0], np.nan, dtype=float)
        finite = ~np.isnan(differences)
        ranks[finite] = stats.rankdata(np.abs(differences[finite]))
        return ranks
    finite = ~np.isnan(matrix)
    ranks = np.full(matrix.shape, np.nan, dtype=float)
    ranks[finite] = stats.rankdata(matrix[finite])
    return ranks


def _column(draws: np.ndarray, index: int) -> np.ndarray:
    """One column of a ``(n_draws, n_replicates, 2)`` stack, masked to its finite entries."""

    values = draws[:, :, index]
    return np.where(np.isnan(values), 0.0, values)


def _present(draws: np.ndarray, index: int) -> np.ndarray:
    return ~np.isnan(draws[:, :, index])


def _mean_and_sd(draws: np.ndarray, index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    present = _present(draws, index)
    n = present.sum(axis=1)
    values = _column(draws, index)
    mean = values.sum(axis=1) / np.maximum(n, 1)
    deviation = np.where(present, draws[:, :, index] - mean[:, None], 0.0)
    sd = np.sqrt((deviation ** 2).sum(axis=1) / np.maximum(n - 1, 1))
    return mean, sd, n


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """``numerator / denominator``, or ``nan`` where the denominator is round-off."""

    usable = np.isfinite(denominator) & (denominator > ROUND_OFF_FLOOR)
    return np.where(usable, numerator / np.where(usable, denominator, 1.0), np.nan)


def _draw_statistics(
    draws: np.ndarray, model: str, ranks: np.ndarray
) -> dict[str, np.ndarray]:
    """Every statistic, for a whole stack of resampled matrices at once.

    The observed matrix is passed through this same function as a single draw, so a
    record's value and the interval around it are computed by identical arithmetic and
    cannot drift apart.
    """

    if model == "paired":
        complete = _present(draws, 0) & _present(draws, 1)
        n = complete.sum(axis=1)
        differences = np.where(complete, draws[:, :, 0] - draws[:, :, 1], np.nan)
        effect = np.nansum(np.where(complete, differences, 0.0), axis=1) / np.maximum(n, 1)
        deviation = np.where(complete, differences - effect[:, None], 0.0)
        sd = np.sqrt((deviation ** 2).sum(axis=1) / np.maximum(n - 1, 1))
        usable = n >= MIN_SAMPLES_PER_LEVEL
        # Cohen's d for a matched design standardises the within-replicate differences by
        # their own SD.  W+ is the sum of the ranks of the *positive* differences; the
        # p-value reported alongside it is scipy's two-sided statistic, which is the
        # smaller tail, so the direction travels with the effect size, not with the p.
        positive = complete & (differences > 0)
        w_plus = (np.where(positive, np.nan_to_num(ranks)[None, :], 0.0)).sum(axis=1)
        return {
            "effect": np.where(usable, effect, np.nan),
            "cohens_d": _safe_ratio(effect, sd),
            "t_statistic": _safe_ratio(effect, sd / np.sqrt(np.maximum(n, 1))),
            "rank_statistic": np.where(usable, w_plus, np.nan),
            "n_level_a": n.astype(float), "n_level_b": n.astype(float),
        }

    mean_a, sd_a, n_a = _mean_and_sd(draws, 0)
    mean_b, sd_b, n_b = _mean_and_sd(draws, 1)
    usable = (n_a >= MIN_SAMPLES_PER_LEVEL) & (n_b >= MIN_SAMPLES_PER_LEVEL)
    effect = mean_a - mean_b
    # Cohen's d across two arms: the SD weighted by degrees of freedom, which is the
    # denominator the OpenProblems functional-marker component uses.
    degrees = np.maximum(n_a + n_b - 2, 1)
    pooled_sd = np.sqrt(((n_a - 1) * sd_a ** 2 + (n_b - 1) * sd_b ** 2) / degrees)
    standard_error = np.sqrt(sd_a ** 2 / np.maximum(n_a, 1) + sd_b ** 2 / np.maximum(n_b, 1))
    # U for level_a: the sum of the pooled ranks it collects, less what it collects by
    # construction.  A large U means level_a sits above level_b.
    present_a = _present(draws, 0)
    rank_sum = (np.where(present_a, np.nan_to_num(ranks)[:, 0][None, :], 0.0)).sum(axis=1)
    return {
        "effect": np.where(usable, effect, np.nan),
        "cohens_d": _safe_ratio(effect, pooled_sd),
        "t_statistic": _safe_ratio(effect, standard_error),
        "rank_statistic": np.where(usable, rank_sum - n_a * (n_a + 1) / 2.0, np.nan),
        "n_level_a": n_a.astype(float), "n_level_b": n_b.astype(float),
    }


def summarise(matrix: np.ndarray, model: str) -> ContrastStatistics:
    """Point estimates from a ``(n_replicates, 2)`` donor-level matrix."""

    values = _draw_statistics(matrix[None, :, :], model, _entry_ranks(matrix, model))
    n_a, n_b = int(values["n_level_a"][0]), int(values["n_level_b"][0])
    if model == "paired":
        a, b = matrix[:n_a, 0], matrix[:n_a, 1]
        rank_p, t_p = signed_rank_p(a, b), paired_p(a, b)
    else:
        a, b = matrix[~np.isnan(matrix[:, 0]), 0], matrix[~np.isnan(matrix[:, 1]), 1]
        rank_p, t_p = rank_sum_p(a, b), welch_p(a, b)
    return ContrastStatistics(
        effect=_first(values["effect"]), cohens_d=_first(values["cohens_d"]),
        t_statistic=_first(values["t_statistic"]), t_pvalue=t_p,
        rank_statistic=_first(values["rank_statistic"]), rank_pvalue=rank_p,
        n_level_a=n_a, n_level_b=n_b,
        note=(
            None if values["effect"][0] == values["effect"][0] else
            (f"only {n_a} complete donor pair(s)" if model == "paired"
             else f"the two levels carry {n_a} and {n_b} replicate(s)")
            + f", {MIN_SAMPLES_PER_LEVEL} each required"
        ),
    )


def _first(values: np.ndarray) -> float | None:
    value = float(values[0])
    return None if not np.isfinite(value) else value


# --------------------------------------------------------------------------------------
# Resampling
# --------------------------------------------------------------------------------------


def donor_block_bootstrap(
    matrix: np.ndarray, model: str, n_resamples: int, rng: np.random.Generator
) -> dict[str, np.ndarray]:
    """Bootstrap every statistic by resampling the replicate unit with replacement.

    The replicate unit is the donor, or the sample when no donor factor is declared.  The
    resampled matrix is fed through the same :func:`_draw_statistics` as the observed one,
    so an interval cannot drift away from the value it is meant to bracket.
    """

    n_replicates = matrix.shape[0]
    draws = matrix[rng.integers(0, n_replicates, size=(n_resamples, n_replicates))]
    return _draw_statistics(draws, model, _entry_ranks(matrix, model))


def restricted_permutation_null(
    matrix: np.ndarray, plan: ContrastPlan, n_resamples: int, rng: np.random.Generator
) -> dict[str, np.ndarray]:
    """Null distribution of every statistic under permutation of the label in-block.

    A paired block is one donor holding both levels, so relabelling inside it swaps that
    donor's two values.  An unpaired block is a batch, and relabelling inside it permutes
    the labels among that batch's own one-armed samples, so the number of each level in the
    batch cannot change and a batch that happens to hold mostly one group cannot
    manufacture a difference.  Permuting across the whole dataset would assume an
    exchangeability the design never promised; permuting cells would treat them as
    replicates.
    """

    complete = ~np.isnan(matrix[:, 0]) & ~np.isnan(matrix[:, 1])
    holds_a = ~np.isnan(matrix[:, 0])
    swaps = np.zeros((n_resamples, matrix.shape[0]), dtype=bool)
    if plan.model == "paired":
        for _, (row,) in plan.blocks:
            swaps[:, row] = rng.random(n_resamples) < 0.5
    else:
        for _, indices in plan.blocks:
            # A replicate that already carries both levels is not relabelled: it has no
            # single label to permute.
            block = [i for i in indices if not complete[i]]
            if len(block) > 1:
                swaps[:, block] = rng.permuted(
                    np.broadcast_to(holds_a[block], (n_resamples, len(block))), axis=1
                )
    draws = np.where(swaps[:, :, None], matrix[None, :, ::-1], matrix[None, :, :])
    return _draw_statistics(draws, plan.model, _entry_ranks(matrix, plan.model))


def permutation_p(observed: float, null: np.ndarray, n_resamples: int) -> float:
    """The +1 empirical p-value: a resampled p is never printed as zero."""

    finite = null[np.isfinite(null)]
    if len(finite) == 0:
        return float("nan")
    as_extreme = int(np.sum(np.abs(finite) >= abs(observed) - 1e-12))
    return float((1 + as_extreme) / (1 + n_resamples))


def clopper_pearson(successes: int, trials: int, alpha: float) -> tuple[float, float]:
    """Exact binomial interval on a Monte Carlo proportion; never an epsilon."""

    if trials < 1:
        return float("nan"), float("nan")
    low = 0.0 if successes == 0 else float(
        stats.beta.ppf(alpha / 2.0, successes, trials - successes + 1)
    )
    high = 1.0 if successes == trials else float(
        stats.beta.ppf(1.0 - alpha / 2.0, successes + 1, trials - successes)
    )
    return low, high


def percentile_interval(values: np.ndarray, alpha: float) -> tuple[float, float] | None:
    """The central ``1 - alpha`` percentile interval, or ``None`` if nothing is finite."""

    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return None
    low, high = np.percentile(finite, [100 * alpha / 2.0, 100 * (1 - alpha / 2.0)])
    return float(low), float(high)


# --------------------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------------------


def _context(plan: ContrastPlan, population: str, marker: str, run_id: str,
             candidate_id: str) -> dict:
    return dict(
        run_id=run_id, candidate_id=candidate_id, marker=marker, population=population,
        comparison_id=plan.comparison_id, n_samples=plan.n_samples, n_donors=plan.n_donors,
        aggregation_weight=1.0,
    )


def _provenance(plan: ContrastPlan, statistics: ContrastStatistics, n_cells: int,
                coverage: float, scope: str) -> str:
    return "; ".join([
        f"model={plan.model}, block_unit={plan.block_unit}, level_a={plan.level_a!r}, "
        f"level_b={plan.level_b!r} (sorted order, so a positive value always means "
        f"{plan.level_a!r} is the higher arm)",
        f"{scope}; the two levels carry {statistics.n_level_a} and "
        f"{statistics.n_level_b} replicate(s) over {n_cells} frozen evaluation cell(s), "
        f"coverage {coverage:.3f}",
        f"cells are reduced to a per-sample mean first and then to a {plan.block_unit} "
        f"mean, so the replicate unit is the {plan.block_unit} and never the cell",
        f"n_resamples={plan.spec.n_resamples}, seed={plan.spec.seed}",
    ] + ([statistics.note] if statistics.note else []))


def _expectation_note(spec: ContrastSpec, interval: tuple[float, float] | None,
                      effect: float) -> str:
    """A machine-greppable verdict on the prespecified expected direction."""

    if spec.expected_direction == "none":
        return "expected_direction=none: no direction was prespecified"
    if spec.expected_direction == "increase_in_a":
        held = interval is None or interval[0] > 0.0
    else:
        held = interval is None or interval[1] < 0.0
    shown = "unavailable" if interval is None else f"[{interval[0]:+.4g}, {interval[1]:+.4g}]"
    return (f"expected_direction={spec.expected_direction}: "
            f"{'supported' if held else 'not_supported'} "
            f"(effect {effect:+.4g}, interval {shown})")


def _record(
    registered: MetricSpec, plan: ContrastPlan, sampling: SamplingManifest,
    population: str, marker: str, run_id: str, candidate_id: str, *,
    value: float | None, interval: tuple[float, float] | None, kind: str | None,
    reason: str, n_cells: int, coverage: float,
) -> MetricRecord:
    low, high = interval if interval is not None else (None, None)
    return metric_record(
        registered, sampling,
        **_context(plan, population, marker, run_id, candidate_id),
        value=None if value is None else float(value),
        status=MetricStatus.OK if value is not None else MetricStatus.UNDEFINED,
        reason=reason, n_cells=n_cells, coverage=coverage,
        uncertainty_low=low, uncertainty_high=high,
        uncertainty_kind=None if kind is None else f"{kind}_{plan.spec.interval_level}",
    )


def _no_value_records(
    plan: ContrastPlan, sampling: SamplingManifest, population: str, marker: str,
    run_id: str, candidate_id: str, provenance: str, n_cells: int, coverage: float,
) -> list[MetricRecord]:
    """Nothing is invented when a group is too thin: no value, and the reason on each row.

    The two effect sizes are ``undefined`` in the mathematical sense (the quantity does
    not exist); the three test rows are ``insufficient_data``, because there was never
    enough data to run the test.  Neither is a failure, and neither is a zero.
    """

    undefined = (BIOLOGICAL_CONTRAST_EFFECT_SIZE, BIOLOGICAL_CONTRAST_COHENS_D)
    records = []
    for registered in CONTRASTAT_RECORDS:
        if registered in undefined:
            records.append(_record(
                registered, plan, sampling, population, marker, run_id, candidate_id,
                value=None, interval=None, kind=None,
                reason=f"{provenance}; the effect size itself is undefined for this group, "
                       "which is reported as undefined and never as a value of zero",
                n_cells=n_cells, coverage=coverage,
            ))
        else:
            records.append(metric_record(
                registered, sampling,
                **_context(plan, population, marker, run_id, candidate_id),
                status=MetricStatus.INSUFFICIENT_DATA,
                reason=f"{provenance}; the test was not run because the effect size does "
                       "not exist for this group, and no substitute was reported",
                n_cells=n_cells, coverage=coverage,
            ))
    return records


def _one_group(
    matrix: np.ndarray, plan: ContrastPlan, sampling: SamplingManifest,
    population: str, marker: str, run_id: str, candidate_id: str,
    n_cells: int, coverage: float, scope: str,
) -> list[MetricRecord]:
    """The records emitted for one (contrast, population, marker, model)."""

    spec = plan.spec
    statistics = summarise(matrix, plan.model)
    provenance = _provenance(plan, statistics, n_cells, coverage, scope)
    if statistics.effect is None:
        return _no_value_records(plan, sampling, population, marker, run_id, candidate_id,
                                 provenance, n_cells, coverage)

    bootstrap = donor_block_bootstrap(
        matrix, plan.model, spec.n_resamples, np.random.default_rng(spec.seed)
    )
    effect_interval = percentile_interval(bootstrap["effect"], spec.alpha)
    d_interval = percentile_interval(bootstrap["cohens_d"], spec.alpha)
    provenance += (
        f"; {int(np.isfinite(bootstrap['effect']).sum())} of {spec.n_resamples} bootstrap "
        "resamples produced a finite effect in this group"
    )
    reason = f"{provenance}; {_expectation_note(spec, effect_interval, statistics.effect)}"

    records = [
        _record(BIOLOGICAL_CONTRAST_EFFECT_SIZE, plan, sampling, population, marker,
                run_id, candidate_id, value=statistics.effect, interval=effect_interval,
                kind=BOOTSTRAP_KIND, reason=f"{reason}; effect = level_a minus level_b",
                n_cells=n_cells, coverage=coverage),
        _record(BIOLOGICAL_CONTRAST_COHENS_D, plan, sampling, population, marker,
                run_id, candidate_id, value=statistics.cohens_d, interval=d_interval,
                kind=BOOTSTRAP_KIND if d_interval is not None else None,
                reason=reason + "; "
                + ("pooled-SD standardisation of the same effect across the two arms"
                   if plan.model == "unpaired"
                   else "standardised by the SD of the within-donor differences"),
                n_cells=n_cells, coverage=coverage),
    ]

    if plan.reachable_permutations < MIN_PERMUTATION_ASSIGNMENTS:
        powerless = (
            f"the restricted permutation reaches only {plan.reachable_permutations} label "
            f"assignment(s) across blocks {[name for name, _ in plan.blocks]}, below the "
            f"{MIN_PERMUTATION_ASSIGNMENTS} a null distribution needs; no p-value and no "
            "null interval are reported rather than a degenerate one"
        )
        for registered in CONTRASTAT_RECORDS[2:]:
            records.append(metric_record(
                registered, sampling,
                **_context(plan, population, marker, run_id, candidate_id),
                status=MetricStatus.INSUFFICIENT_DATA, reason=f"{provenance}; {powerless}",
                n_cells=n_cells, coverage=coverage,
            ))
        return records

    nulls = restricted_permutation_null(
        matrix, plan, spec.n_resamples, np.random.default_rng(spec.seed + 1)
    )
    p_value = permutation_p(statistics.effect, nulls["effect"], spec.n_resamples)
    null_interval = percentile_interval(nulls["effect"], spec.alpha)
    finite_null = nulls["effect"][np.isfinite(nulls["effect"])]
    as_extreme = int(np.sum(np.abs(finite_null) >= abs(statistics.effect) - 1e-12))
    p_interval = clopper_pearson(as_extreme, len(finite_null), spec.alpha)
    permutation_note = "; ".join([
        f"the {plan.level_a}/{plan.level_b} label was permuted inside "
        f"{[name for name, _ in plan.blocks]} and never across the dataset",
        f"permutation p={p_value:.4g} on the effect over {len(finite_null)} finite "
        f"resample(s) of {spec.n_resamples}, {as_extreme} of them at least as extreme",
        "null interval of the effect "
        + (f"[{null_interval[0]:+.4g}, {null_interval[1]:+.4g}]" if null_interval
           else "unavailable"),
    ])
    paired = plan.model == "paired"
    analytic_t = (
        f"analytic {'paired t' if paired else 'Welch t'} "
        f"p={statistics.t_pvalue:.4g}" if statistics.t_pvalue is not None
        else f"the analytic {'paired t' if paired else 'Welch t'} statistic is undefined "
             "(zero variance in one arm)"
    )
    analytic_rank = (
        f"analytic {'signed-rank' if paired else 'rank-sum'} p={statistics.rank_pvalue:.4g}"
        if statistics.rank_pvalue is not None
        else f"the analytic {'signed-rank' if paired else 'rank-sum'} statistic is "
             "undefined (no replicate pair is discordant)"
    )

    records.append(_record(
        BIOLOGICAL_CONTRAST_TTEST_STATISTIC, plan, sampling, population, marker,
        run_id, candidate_id, value=statistics.t_statistic,
        interval=percentile_interval(nulls["t_statistic"], spec.alpha),
        kind=PERMUTATION_NULL_KIND, reason=f"{provenance}; {analytic_t}; {permutation_note}",
        n_cells=n_cells, coverage=coverage,
    ))
    records.append(_record(
        BIOLOGICAL_CONTRAST_RANK_STATISTIC, plan, sampling, population, marker,
        run_id, candidate_id, value=statistics.rank_statistic,
        interval=percentile_interval(nulls["rank_statistic"], spec.alpha),
        kind=PERMUTATION_NULL_KIND,
        reason=f"{provenance}; {analytic_rank}; the reported statistic is "
               + ("the sum of the ranks of the positive within-donor differences"
                  if paired else "U for level_a")
               + f", so a larger value means {plan.level_a!r} sits higher; "
               + permutation_note,
        n_cells=n_cells, coverage=coverage,
    ))
    records.append(_record(
        BIOLOGICAL_CONTRAST_PERMUTATION_P, plan, sampling, population, marker,
        run_id, candidate_id, value=p_value, interval=p_interval, kind=PERMUTATION_P_KIND,
        # This record is the one that most needs its provenance: the p-value is
        # meaningless without the resample count, the seed, and the fact that the
        # replicate unit is the block and not the cell.  Carrying the full
        # provenance here also keeps every record in the family self-describing,
        # so a row read in isolation states how it was produced.
        reason=f"{provenance}; {permutation_note}", n_cells=n_cells, coverage=coverage,
    ))
    return records


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------


def _refusal_record(
    refusal: ContrastRefusal, sampling: SamplingManifest, run_id: str, candidate_id: str
) -> MetricRecord:
    return metric_record(
        BIOLOGICAL_CONTRAST_ADMISSIBILITY, sampling,
        run_id=run_id, candidate_id=candidate_id, comparison_id=refusal.comparison_id,
        status=MetricStatus.INELIGIBLE,
        reason=f"{refusal.reason_code}: {refusal.reason}",
        n_samples=0, n_donors=0, n_cells=0, coverage=0.0, aggregation_weight=1.0,
    )


def _admissible_record(
    plan: ContrastPlan, sampling: SamplingManifest, run_id: str, candidate_id: str
) -> MetricRecord:
    return metric_record(
        BIOLOGICAL_CONTRAST_ADMISSIBILITY, sampling,
        run_id=run_id, candidate_id=candidate_id, comparison_id=plan.comparison_id,
        status=MetricStatus.OK, value=1.0, coverage=1.0, aggregation_weight=1.0,
        n_samples=plan.n_samples, n_donors=plan.n_donors,
        reason=(
            f"admissible: {plan.spec.factor!r} is a declared protected factor, the design "
            f"audit does not call it non_identifiable, and the design supports a "
            f"{plan.model} model with the {plan.block_unit} as the replicate unit "
            f"({plan.n_donors} replicate(s), {plan.reachable_permutations} reachable "
            "restricted label assignments)"
        ),
    )


def compute_preservation_metrics(
    dataset: Dataset,
    candidate_matrix: np.ndarray,
    sampling: SamplingManifest,
    design: DesignAudit,
    contrasts: list[ContrastSpec],
    *,
    run_id: str,
    candidate_id: str,
) -> MetricResult:
    """Biological preservation for one candidate, on the frozen partition.

    A refused contrast becomes one ``biological_contrast_admissibility`` record carrying
    its reason, so a contrast that could not be run is visible in the record set instead
    of absent from it.  Every record that reports a value also reports the interval and
    the method that produced it.
    """

    records: list[MetricRecord] = []
    frozen = _frozen_rows_by_sample(dataset, sampling)
    feature_idx = {str(name): i for i, name in enumerate(dataset.adata.var_names)}

    for spec in contrasts:
        outcome = resolve_contrast(dataset, design, sampling, spec)
        if isinstance(outcome, ContrastRefusal):
            records.append(_refusal_record(outcome, sampling, run_id, candidate_id))
            continue
        plan = outcome
        records.append(_admissible_record(plan, sampling, run_id, candidate_id))

        for population in plan.populations:
            by_sample = frozen.get(population, {})
            present = [sid for sid in plan.contributing if sid in by_sample]
            usable = [
                sid for sid in present if len(by_sample[sid]) >= MIN_CELLS_PER_SAMPLE_MEAN
            ]
            scope = (
                f"population {population!r} contributes {len(usable)} of "
                f"{plan.n_samples} sample(s) with at least "
                f"{MIN_CELLS_PER_SAMPLE_MEAN} frozen evaluation cells"
            )
            if len(present) < len(plan.contributing):
                scope += (
                    f"; {len(plan.contributing) - len(present)} carry no cell of this "
                    "population on the frozen partition"
                )
            if len(present) > len(usable):
                scope += f"; {len(present) - len(usable)} fell below that cell floor"
            n_cells = sum(len(by_sample[sid]) for sid in usable)
            coverage = len(usable) / plan.n_samples if plan.n_samples else 0.0
            for marker in plan.markers:
                column = feature_idx[marker]
                values = np.array(
                    [candidate_matrix[np.asarray(by_sample[sid], dtype=int), column].mean()
                     for sid in usable], dtype=float
                )
                matrix = donor_level_matrix(plan, usable, values)
                records.extend(_one_group(
                    matrix, plan, sampling, population, marker, run_id, candidate_id,
                    n_cells, coverage, scope,
                ))

    status = "ok" if all(
        r.status in (MetricStatus.OK, MetricStatus.INELIGIBLE,
                     MetricStatus.INSUFFICIENT_DATA, MetricStatus.UNDEFINED)
        for r in records
    ) else "partial"
    return MetricResult(
        run_id=run_id, candidate_id=candidate_id, metric_group="biological_preservation",
        records=records, status=status,
        reason=None if records else "no biological contrast was prespecified for this run",
    )
