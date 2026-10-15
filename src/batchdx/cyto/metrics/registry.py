"""Metric registry.

Namespaces are kept apart on purpose:

``native``
    produced by running the upstream implementation itself.

``compatibility``
    a local implementation that has been checked against a pinned upstream
    implementation within a documented tolerance.

``common``
    written here, computed on the frozen common partition.  A metric in this namespace
    is never described as an upstream component, whatever it resembles.
"""

from __future__ import annotations

from dataclasses import dataclass

from batchdx.cyto.schemas.enums import MetricDirection


@dataclass(frozen=True)
class MetricSpec:
    """Everything a reader needs to interpret one metric."""

    name: str
    namespace: str
    version: str
    stage: str
    grouping: str
    formula: str
    implementation: str
    direction: MetricDirection
    scale: str
    eligible_population: str
    aggregation: str
    requirements: tuple[str, ...] = ()
    notes: str = ""
    upstream_verified: bool = False


_REGISTRY: dict[tuple[str, str], MetricSpec] = {}


def register_metric(spec: MetricSpec) -> MetricSpec:
    _REGISTRY[(spec.namespace, spec.name)] = spec
    return spec


def get_metric(namespace: str, name: str) -> MetricSpec:
    key = (namespace, name)
    if key not in _REGISTRY:
        raise KeyError(f"unknown metric {namespace}.{name}; registered: {sorted(_REGISTRY)}")
    return _REGISTRY[key]


def registered_metrics() -> dict[tuple[str, str], MetricSpec]:
    return dict(_REGISTRY)


# --------------------------------------------------------------------------------------
# Common-namespace metrics implemented in this build
# --------------------------------------------------------------------------------------

EMD_WITHIN_POPULATION = register_metric(MetricSpec(
    name="emd_within_population",
    namespace="common",
    version="1.0.0",
    stage="distribution_alignment",
    grouping="(population, marker, batch pair) on the frozen common partition",
    formula="first Wasserstein distance between the two batches' empirical distributions",
    implementation="scipy.stats.wasserstein_distance on the declared transformed scale",
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="units of the declared expression scale",
    eligible_population="populations with at least the declared minimum cells in both batches",
    aggregation="median across comparisons, equally weighted per (population, marker, pair)",
    requirements=("frozen evaluation labels", "at least one eligible batch pair"),
    notes=(
        "an empirical-distance implementation in data units. It is deliberately named "
        "differently from cyCombine's binned EMD, which works in bin units on histogram "
        "densities and is a different number."
    ),
))

MAD_SIGNED_CHANGE = register_metric(MetricSpec(
    name="mad_signed_change",
    namespace="common",
    version="1.0.0",
    stage="dispersion",
    grouping="(population, marker, batch)",
    formula="MAD(candidate) - MAD(baseline), with the sign retained",
    implementation="scipy.stats.median_abs_deviation with scale='normal' (constant 1.4826)",
    direction=MetricDirection.NEUTRAL,
    scale="units of the declared expression scale",
    eligible_population="populations with at least the declared minimum cells in the batch",
    aggregation="median across groups; negative values mean the candidate compressed spread",
    requirements=("frozen evaluation labels",),
    notes=(
        "signed on purpose: cyCombine's MAD score takes an absolute difference, which "
        "cannot distinguish a correction that restores spread from one that destroys it."
    ),
))

MAD_COMPRESSION = register_metric(MetricSpec(
    name="mad_compression_fraction",
    namespace="common",
    version="1.0.0",
    stage="dispersion",
    grouping="(population, marker, batch)",
    formula="fraction of groups whose MAD fell by more than the declared tolerance",
    implementation="derived from mad_signed_change",
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="fraction in [0, 1]",
    eligible_population="as mad_signed_change",
    aggregation="single value per candidate",
    requirements=("frozen evaluation labels",),
))

POPULATION_FREQUENCY_SHIFT = register_metric(MetricSpec(
    name="population_frequency_shift",
    namespace="common",
    version="1.0.0",
    stage="composition",
    grouping="(population, batch pair)",
    formula="absolute difference of mean sample-level frequency between two batches",
    implementation="sample-level frequencies averaged within batch, then differenced",
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="fraction of evaluation cells",
    eligible_population="all frozen populations",
    aggregation="median across pairs, equally weighted per sample",
    requirements=("frozen evaluation labels",),
    notes="reported separately from intensity alignment; a frequency shift is not an "
          "intensity artefact",
))

JOINT_STRUCTURE = register_metric(MetricSpec(
    name="joint_structure_distance",
    namespace="common",
    version="1.0.0",
    stage="joint_structure",
    grouping="(population, batch)",
    formula=(
        "mean absolute difference between the candidate's and the baseline's "
        "within-population marker correlation matrices (upper triangle)"
    ),
    implementation="numpy Pearson correlation over evaluation markers within population",
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="correlation units in [0, 2]",
    eligible_population="populations with at least the declared minimum cells",
    aggregation="median across (population, batch)",
    requirements=("frozen evaluation labels", "at least two evaluation markers"),
    notes=(
        "required so that a corruption which preserves every marginal distribution while "
        "destroying within-population co-expression is detectable at all."
    ),
))

EMD_REDUCTION = register_metric(MetricSpec(
    name="emd_relative_reduction",
    namespace="common",
    version="1.0.0",
    stage="distribution_alignment",
    grouping="candidate versus the identity candidate, over identical comparisons",
    formula="(sum EMD_identity - sum EMD_candidate) / sum EMD_identity",
    implementation="computed over the intersection of comparisons that succeeded for both",
    direction=MetricDirection.HIGHER_IS_BETTER,
    scale="fraction",
    eligible_population="as emd_within_population",
    aggregation="single value per candidate",
    requirements=("the identity candidate must have been evaluated",),
    notes=(
        "numerator and denominator use the same eligible comparisons. A near-zero "
        "denominator returns status 'undefined'; it is never stabilised with an epsilon."
    ),
))

BIOLOGICAL_CONTRAST = register_metric(MetricSpec(
    name="biological_contrast_change",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(prespecified contrast, population) on the frozen common partition",
    formula=(
        "effect_after - effect_before, where the effect is the group-to-group difference "
        "of donor-level (or sample-level) marker means; the contrast is preserved when "
        "|change| is within the predeclared preservation tolerance"
    ),
    implementation=(
        "donor-paired difference or donor-level pooled difference (numpy), as the design "
        "declares; the model used is stated in the record reason"
    ),
    direction=MetricDirection.NEUTRAL,
    scale="units of the declared expression scale",
    eligible_population="populations with enough cells and enough independent units in "
                        "both groups",
    aggregation="one record per prespecified contrast and population; never aggregated "
                "across contrasts",
    requirements=("declared biological contrasts in the configuration", "frozen "
                  "evaluation labels"),
    notes=(
        "a contrast confounded with the correction factor is reported unassessable, "
        "never preserved. The verdict (preserved / not preserved / unassessable) lives "
        "in the record reason together with the tolerance and its declared origin."
    ),
))


# --------------------------------------------------------------------------------------
# Biological preservation (Stage 8)
#
# Every metric below is computed on *per-donor* summaries of *per-sample* means.  A cell
# is not a replicate: ten thousand cells from one tube are one observation, and a test
# that treated them as ten thousand would report an effect size the design cannot
# support.  Each record carries a resampling interval, so a preservation claim is always
# a claim with a stated width.
# --------------------------------------------------------------------------------------

BIOLOGICAL_CONTRAST_ADMISSIBILITY = register_metric(MetricSpec(
    name="biological_contrast_admissibility",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="one record per prespecified contrast, per candidate",
    formula="(no value; the record exists to say whether the contrast may be run at all)",
    implementation="validation of the contrast against the sample table and the design audit",
    direction=MetricDirection.NEUTRAL,
    scale="admissible or not",
    eligible_population="all",
    aggregation="one record per contrast",
    requirements=("frozen evaluation labels", "design audit"),
    notes=(
        "a contrast whose levels are absent, whose factor is not a declared protected "
        "factor, whose design the audit calls non_identifiable, or whose declared "
        "covariate is nested inside one of its own levels is REFUSED with a reason here. "
        "It is never silently dropped, and it never reports a score of zero."
    ),
))

BIOLOGICAL_CONTRAST_EFFECT_SIZE = register_metric(MetricSpec(
    name="biological_contrast_effect_size",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(contrast, population, marker, model)",
    formula="mean(level_a) - mean(level_b) on the donor-level summaries of per-sample means",
    implementation=(
        "numpy over per-sample means of the frozen evaluation cells, averaged within donor"
    ),
    direction=MetricDirection.NEUTRAL,
    scale="units of the declared expression scale",
    eligible_population="populations with at least the declared minimum cells in a sample",
    aggregation="one value per (contrast, population, marker, model); never pooled over markers",
    requirements=("frozen evaluation labels", "a prespecified admissible contrast"),
    notes=(
        "the level order is the sorted order, so a positive value means the same thing in "
        "every run. uncertainty_low/high are the percentile interval from resampling the "
        "block unit (donor where a donor factor exists, otherwise sample) with replacement."
    ),
))

BIOLOGICAL_CONTRAST_COHENS_D = register_metric(MetricSpec(
    name="biological_contrast_cohens_d",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(contrast, population, marker, model)",
    formula=(
        "unpaired: (mean_a - mean_b) / pooled SD, with the SD weighted by degrees of "
        "freedom; paired: mean(a - b) / SD(a - b) over complete donor pairs"
    ),
    implementation=(
        "the same standardised-mean-difference form the OpenProblems functional-marker "
        "component computes, evaluated on this package's donor-level summaries"
    ),
    direction=MetricDirection.NEUTRAL,
    scale="standardised units",
    eligible_population="as biological_contrast_effect_size",
    aggregation="one value per (contrast, population, marker, model)",
    requirements=("frozen evaluation labels", "a prespecified admissible contrast"),
    notes=(
        "direction-consistent because the two levels are sorted, not taken in the order "
        "they happen to appear. A zero pooled SD is UNDEFINED, never zero: there is no "
        "effect size to report when both arms are constant."
    ),
))

BIOLOGICAL_CONTRAST_TTEST_STATISTIC = register_metric(MetricSpec(
    name="biological_contrast_ttest_statistic",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(contrast, population, marker, model)",
    formula="Welch's t when the samples are unpaired, the paired t when they are matched",
    implementation="scipy.stats.ttest_ind(equal_var=False) or scipy.stats.ttest_rel on "
                   "donor-level summaries",
    direction=MetricDirection.NEUTRAL,
    scale="t units",
    eligible_population="as biological_contrast_effect_size",
    aggregation="one value per (contrast, population, marker, model)",
    requirements=("frozen evaluation labels", "a prespecified admissible contrast"),
    notes=(
        "the analytic p-value is reported in the record's reason; the interval on the "
        "statistic is the restricted-permutation null, because a parametric interval on "
        "n=6 donors is not an honest uncertainty statement."
    ),
))

BIOLOGICAL_CONTRAST_RANK_STATISTIC = register_metric(MetricSpec(
    name="biological_contrast_rank_statistic",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(contrast, population, marker, model)",
    formula=(
        "unpaired: Mann-Whitney U; paired: the Wilcoxon signed-rank W, reduced to the "
        "larger tail so that a positive value always means 'level_a shifted upwards'"
    ),
    implementation="scipy.stats.mannwhitneyu or scipy.stats.wilcoxon on donor-level summaries",
    direction=MetricDirection.NEUTRAL,
    scale="rank units",
    eligible_population="as biological_contrast_effect_size",
    aggregation="one value per (contrast, population, marker, model)",
    requirements=("frozen evaluation labels", "a prespecified admissible contrast"),
    notes=(
        "the distribution-free counterpart of the t statistic, included because a "
        "Wilcoxon rank-sum over two or three samples per arm cannot assume normality. "
        "With all paired differences exactly zero the statistic is UNDEFINED."
    ),
))

BIOLOGICAL_CONTRAST_PERMUTATION_P = register_metric(MetricSpec(
    name="biological_contrast_permutation_p",
    namespace="common",
    version="1.0.0",
    stage="biological_preservation",
    grouping="(contrast, population, marker, model)",
    formula="(1 + #{resamples whose |effect| >= the observed |effect|}) / (1 + n_resamples)",
    implementation=(
        "the group label is permuted inside its own block (donor when matched, batch "
        "otherwise) with numpy.random.default_rng(seed); never across cells"
    ),
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="probability in [0, 1]",
    eligible_population="as biological_contrast_effect_size",
    aggregation="one value per (contrast, population, marker, model)",
    requirements=(
        "frozen evaluation labels", "a prespecified admissible contrast",
        "enough block sizes for the null to be reachable",
    ),
    notes=(
        "the +1 form means the smallest reportable value is 1/(1+n_resamples) and a p of "
        "zero is never printed. The interval is the Clopper-Pearson interval on the Monte "
        "Carlo estimate, which is the only uncertainty a resampled p-value has."
    ),
))

NATIVE_CYCOMBINE_EMD = register_metric(MetricSpec(
    name="cycombine_evaluate_emd",
    namespace="native",
    version="upstream",
    stage="distribution_alignment",
    grouping="cyCombine's own partition (label column as supplied to evaluate_emd)",
    formula="cyCombine::evaluate_emd; reduction = sum(Reduction) / sum(Uncorrected) after "
            "filtering pairs below filter_limit in both states",
    implementation="the upstream R function, executed as-is",
    direction=MetricDirection.HIGHER_IS_BETTER,
    scale="bin units; reduction is a fraction",
    eligible_population="cyCombine's own labels",
    aggregation="upstream",
    requirements=("R", "cyCombine", "emdist"),
    notes=(
        "native compatibility mode only. cyCombine's benchmarking procedure lets each "
        "correction re-cluster, which changes the uncorrected baseline; that behaviour is "
        "preserved here and must not be mixed with the common-partition numbers."
    ),
))

NATIVE_CYCOMBINE_MAD = register_metric(MetricSpec(
    name="cycombine_evaluate_mad",
    namespace="native",
    version="upstream",
    stage="dispersion",
    grouping="cyCombine's own partition",
    formula="median of |MAD(uncorrected) - MAD(corrected)| over (population, batch, marker)",
    implementation="the upstream R function, executed as-is",
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="units of the declared expression scale",
    eligible_population="cyCombine's own labels",
    aggregation="upstream",
    requirements=("R", "cyCombine"),
    notes="unsigned by construction; pair it with common.mad_signed_change to tell "
          "restoration from compression",
))

# --------------------------------------------------------------------------------------
# OpenProblems components
#
# The metrics below are the ones this build computes itself, on the frozen common
# partition, from the Viash components pinned at commit 37151479.  The distinction
# between the three of them is not cosmetic and is not decided here:
#
# *   `compatibility` means the pinned helper was *run* on a shared fixture and this
#     code produced the same number, within the tolerance recorded in
#     docs/validation_results.md.  scripts/audit_upstream.py fails the build if a row
#     claims that without upstream_verified.
# *   `common` means the *form* was reproduced here without ever running the pinned
#     code, and the spec says so in `notes`.  It is never called a component.
# --------------------------------------------------------------------------------------

EMD_MEAN_CT_HORIZ = register_metric(MetricSpec(
    name="emd_mean_ct_horiz",
    namespace="compatibility",
    version="37151479",
    stage="distribution_alignment",
    grouping="(donor, batch pair, population), markers pooled within a comparison",
    formula=(
        "for each cell type and marker, the Wasserstein distance between two "
        "distributions obtained by histogramming the expression into 2000 bins over "
        "[-100, 100.1); the metric is the mean of those distances"
    ),
    implementation=(
        "extracted from src/metrics/emd/helper.py:bin_array and :compute_emd at the "
        "pinned commit, and compared against the pinned helper on a shared fixture"
    ),
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="bin units; 1.0 is one bin width",
    eligible_population="populations with at least 50 frozen evaluation cells on each side",
    aggregation="mean over every (donor, batch pair, population, marker) distance",
    requirements=("frozen evaluation labels", "a donor measured in two batches"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy. Upstream compares the two halves of one donor's split; "
        "the frozen partition has no split, so a donor's cells are grouped by batch and "
        "every pair of batches that donor appears in is compared. A donor measured in "
        "only one batch yields insufficient_data with that reason, never a zero. The "
        "50-cell-per-side threshold, the bin grid, the weighted Wasserstein call and the "
        "mean-over-markers aggregation are the pinned ones, and the tests run the pinned "
        "helper itself and compare. Upstream reports a point statistic with no interval, "
        "so none is attached here."
    ),
))

EMD_MEAN_CT_VERT = register_metric(MetricSpec(
    name="emd_mean_ct_vert",
    namespace="compatibility",
    version="37151479",
    stage="distribution_alignment",
    grouping="(batch, sample pair, population), markers pooled within a comparison",
    formula=(
        "the same binned Wasserstein distance as emd_mean_ct_horiz, over every pair of "
        "samples from the same group"
    ),
    implementation=(
        "extracted from src/metrics/emd/helper.py:bin_array, :compute_emd and the "
        "combination loop of :get_vert_emd_for_integrated_adata at the pinned commit, and "
        "compared against the pinned helper on a shared fixture"
    ),
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="bin units; 1.0 is one bin width",
    eligible_population="populations with at least 50 frozen evaluation cells on each side",
    aggregation="mean over every (batch, sample pair, population, marker) distance",
    requirements=("frozen evaluation labels", "two samples in one batch"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy. Upstream's group is a biological condition (WT versus "
        "KO) and it compares every pair of samples within one; the frozen partition has no "
        "group factor, so the batch stands in for it. That is a harder test, not an "
        "easier one: a correction which only aligns across batches leaves a within-batch "
        "sample pair untouched, so this component penalises a failure mode the horizontal "
        "one cannot see. The distance, the bin grid and the 50-cell threshold are the "
        "pinned ones. Upstream raises on an all-thin input; here that is "
        "insufficient_data with a reason."
    ),
))

RATIO_CONSISTENT_PEAKS = register_metric(MetricSpec(
    name="ratio_consistent_peaks",
    namespace="compatibility",
    version="37151479",
    stage="distributional_shape",
    grouping="(donor, population, marker), one case per triple",
    formula=(
        "case 1 / (case 1 + case 2), where case 1 is a triple whose expression peaks "
        "number the same in every batch of the donor after correction as they did in "
        "every batch before it, and case 2 is a triple that agreed before and does not "
        "after -- including one the correction flattened to a constant; triples that "
        "already disagreed are excluded from the ratio"
    ),
    implementation=(
        "extracted from src/metrics/ratio_consistent_peaks/helper.py:standardise_marker_"
        "expression, :get_kde_density and :call_peaks at the pinned commit, with the case "
        "classification transcribed from that component's script.py, and compared against "
        "the pinned helper on a shared fixture"
    ),
    direction=MetricDirection.HIGHER_IS_BETTER,
    scale="fraction in [0, 1]",
    eligible_population=(
        "populations with at least 100 frozen evaluation cells in every batch of the "
        "donor, which is the threshold in the pinned script"
    ),
    aggregation="a single value per candidate, with one case row per triple",
    requirements=("frozen evaluation labels", "a donor measured in two batches"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy. Upstream compares two batches and tests equality between "
        "them; a donor here may span more than two, so consistency means all of them "
        "agreeing. Two honest limits. First, only the helper is verified against an "
        "upstream run: the component's script.py reads .h5ad resources the pinned "
        "checkout does not carry, so the case classification here is a transcription that "
        "the tests exercise directly. Second, peak calling at height=0.1 is an absolute "
        "threshold on a KDE density, so on data with no clearly bimodal marker every "
        "triple reads as zero peaks and the ratio is 1.0 for a reason that has nothing to "
        "do with the correction. The per-triple rows carry the case, so a trivial 1.0 is "
        "visible rather than hidden in the aggregate."
    ),
))

AVERAGE_BATCH_R2_CT = register_metric(MetricSpec(
    name="average_batch_r2_ct",
    namespace="compatibility",
    version="37151479",
    stage="distribution_alignment",
    grouping="(donor, population, marker), markers pooled within a donor and population",
    formula=(
        "for each donor, population and marker, 1 - SS_res/SS_tot of a least-squares fit "
        "of the expression on the batch dummies; the metric is the mean of those fits"
    ),
    implementation=(
        "extracted from src/metrics/average_batch_r2/helper.py:concat_paired_samples and "
        ":fit_r2 (scikit-learn's LinearRegression on the batch dummies, read through "
        "model.score) at the pinned commit, and compared against the pinned helper on a "
        "shared fixture"
    ),
    direction=MetricDirection.LOWER_IS_BETTER,
    scale="fraction of variance; numerically it may fall below zero",
    eligible_population="populations with at least 20 frozen evaluation cells in each batch",
    aggregation="mean over every (donor, population, marker) fit",
    requirements=("frozen evaluation labels", "a donor measured in two batches", "scikit-learn"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy. Upstream walks donors, subsets one (donor, cell type) "
        "view and requires at least 20 cells in each *split*; the frozen partition has no "
        "split, so a donor's cells are pooled across the batches it appears in and the "
        "threshold is per batch. The fit itself is the pinned one, dummies included: a "
        "donor spanning three batches is fitted with two dummy columns, as "
        "pd.get_dummies(drop_first=True) produces, not with a single contrast. The two "
        "assertions the pinned helper makes -- equal marker counts, identical donor sets "
        "in the two splits -- have no analogue here and are not reproduced. A marker that "
        "is constant within the donor is undefined, never 0.0: the batch covariate cannot "
        "explain any of a zero total variance."
    ),
))

I_LISI = register_metric(MetricSpec(
    name="iLisi",
    namespace="compatibility",
    version="37151479",
    stage="local_structure",
    grouping="one record per candidate, over every frozen evaluation cell",
    formula=(
        "(median per-cell LISI - 1) / (n_batches - 1), where a cell's LISI is the "
        "inverse Simpson index of the batch labels among its 100 nearest neighbours"
    ),
    implementation=(
        "extracted from src/metrics/lisi/helper.py:compute_ilisi and "
        ":_check_batch_group_confounding at the pinned commit, over "
        "scib_metrics.nearest_neighbors.pynndescent and scib_metrics.lisi_knn, and "
        "compared against the pinned helper on a shared fixture"
    ),
    direction=MetricDirection.HIGHER_IS_BETTER,
    scale="0 to 1; 0 = every neighbourhood is one batch, 1 = uniformly mixed",
    eligible_population=(
        "all frozen evaluation cells, and at least 101 of them so a 100-neighbour "
        "graph exists"
    ),
    aggregation="one value per candidate: the median over cells, normalised",
    requirements=("frozen evaluation labels", "scib-metrics", "pynndescent"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy, and this one is stricter than upstream. The pinned helper "
        "tests whether batch is confounded by obs['group'], a biological condition the "
        "frozen partition does not carry; the frozen population stands in for it, so the "
        "check here is whether the batch is confounded by the population, and the "
        "confounded population is named on the record. A gating population is finer than a "
        "study group, so the substitution can refuse where the pinned component would have "
        "answered: on the paired-reference fixture a single 'ungated' cell in one batch is "
        "enough to withhold the metric, which is a property of the labelling, not of the "
        "correction. Reporting a number there would be reporting the design. The pinned helper returns NaN for the whole metric when any group is "
        "confounded; here that is 'undefined' with the reason, never a zero, and a "
        "confounded run cannot be reported as a candidate that mixes well. The pinned "
        "seed (42) is kept, which also makes the approximate neighbour search "
        "reproducible. The pinned script averages the two splits; the frozen partition "
        "holds one candidate, so this is the single-split computation. Upstream reports a "
        "point statistic and no interval, so none is attached."
    ),
))

C_LISI = register_metric(MetricSpec(
    name="cLisi",
    namespace="compatibility",
    version="37151479",
    stage="local_structure",
    grouping="one record per candidate, over every frozen evaluation cell",
    formula=(
        "(n_populations - median per-cell LISI) / (n_populations - 1), with the "
        "population labels in place of the batch labels"
    ),
    implementation=(
        "extracted from src/metrics/lisi/helper.py:compute_clisi at the pinned commit, "
        "over scib_metrics.nearest_neighbors.pynndescent and scib_metrics.lisi_knn, and "
        "compared against the pinned helper on a shared fixture"
    ),
    direction=MetricDirection.HIGHER_IS_BETTER,
    scale="0 to 1; 0 = populations fully mixed, 1 = neighbourhoods pure",
    eligible_population=(
        "all frozen evaluation cells, at least 101 of them, and at least two populations"
    ),
    aggregation="one value per candidate: the median over cells, normalised",
    requirements=("frozen evaluation labels", "scib-metrics", "pynndescent"),
    upstream_verified=True,
    notes=(
        "ADAPTATION, not a copy, and the same two as iLisi: the frozen population is the "
        "component's cell_type, and there is one candidate rather than two splits. The "
        "pinned seed for this component is 0, not 42, and keeping it is what makes the "
        "approximate neighbour search reproducible. A partition with a single population "
        "has no denominator and is undefined, not 0.0. iLisi and cLisi answer different "
        "questions -- whether batches mixed and whether populations stayed apart -- and "
        "one of them improving by destroying the other is a distortion, not a correction."
    ),
))

#: Metric names that exist in the OpenProblems task at the pinned commit 37151479.  The
#: nine are reconciled one by one: each is either registered here, reimplemented under a
#: name of its own, or carries a row in :data:`DEFERRED_METRICS` saying why not.  A name
#: in none of those is a row someone forgot.
OPENPROBLEMS_METRIC_NAMES: tuple[str, ...] = (
    "openproblems.average_batch_r2_ct",
    "openproblems.cLisi",
    "openproblems.emd_mean_ct_horiz",
    "openproblems.emd_mean_ct_vert",
    "openproblems.flowsom_mean_mapping_similarity",
    "openproblems.functional_marker_preservation_cohens_d",
    "openproblems.functional_marker_preservation_wilcoxon",
    "openproblems.iLisi",
    "openproblems.ratio_consistent_peaks",
)

#: Upstream metric names this build computes itself, mapped to the local spec that does
#: the work.  A reimplementation listed here is a local metric, never a component.
#: Empty at present: the only candidate, ``average_batch_r2_ct``, stopped being a
#: reimplementation once scikit-learn became available and was promoted to a verified
#: ``compatibility`` row.  The mapping stays because the next reimplementation should be
#: recorded here rather than given the upstream metric's name.
REIMPLEMENTED_UPSTREAM_METRICS: dict[str, str] = {}

#: Metrics named in the specification that this build does not implement.  They are
#: listed so that the report can say "unavailable" instead of silently omitting them.
DEFERRED_METRICS: dict[str, str] = {
    "openproblems.functional_marker_preservation_wilcoxon":
        "not implemented, and the earlier description of this row was wrong. The pinned "
        "script.py does not report a Wilcoxon statistic: it reports a RETENTION RATIO, "
        "the fraction of (functional marker, cell type) pairs still significant at "
        "p <= 0.1 in both integrated splits, out of those significant in both "
        "uncorrected batches. common.biological_contrast_rank_statistic is the "
        "Mann-Whitney U that feeds that test, not its result, so calling it an "
        "implementation of this row was a category error. Two things are missing and "
        "neither is code: a significance baseline (the ratio's denominator needs an "
        "uncorrected run significant in BOTH batches, and this package evaluates one "
        "candidate against a frozen partition, not two), and the third factor the "
        "component filters on (var['marker_type'] == 'functional', which the frozen "
        "partition's marker set does not carry). The shared statistic is verified "
        "against the pinned helper by "
        "tests/scientific/test_stage7_openproblems.py::"
        "test_the_pinned_wilcoxon_helper_matches_the_local_rank_statistic.",
    "openproblems.functional_marker_preservation_cohens_d":
        "not implemented, and the earlier description of this row was wrong in the same "
        "way. The pinned script.py reports mean(|d_split1 - d_unintegrated|, "
        "|d_split2 - d_unintegrated|) over the pairs significant in both uncorrected "
        "batches -- effect-size DRIFT against an uncorrected baseline -- not a Cohen's d. "
        "common.biological_contrast_cohens_d computes the per-pair d, which the pinned "
        "helper._cohens_d also computes and which the tests now verify against it, but a "
        "per-pair statistic is not the drift. The missing inputs are the same as for the "
        "Wilcoxon row: two integrated splits and an uncorrected baseline significance "
        "set. Reporting the statistic as the metric would have hidden the difference.",
    "openproblems.flowsom_mean_mapping_similarity":
        "BLOCKED IN renv/, AND NOT COMPUTABLE FROM THE PINNED SOURCE ALONE. Two separate "
        "obstacles, both still true. (1) FlowSOM is not installable in renv/: it is not "
        "on CRAN (available.packages() has no FlowSOM row), it is Bioconductor-only, "
        "renv/ has no BiocManager, and its flowCore dependency is likewise "
        "Bioconductor-only; adding Bioconductor repositories to the isolated renv "
        "library is a central environment change, not a metrics change. (2) Even with "
        "FlowSOM installed, the metric needs three benchmark parameters read from the "
        "unintegrated file (uns['parameter_num_clusters'], uns['parameter_som_xdim'], "
        "uns['parameter_som_ydim']) and two FlowSOM trees per donor with cells mapped "
        "between them; this package's frozen partition has one candidate and no benchmark "
        "parameter, so there is no honest way to produce the number without inventing "
        "those three values. The dissimilarity arithmetic in helper.R (weighted sum of "
        "absolute cell-type proportion differences between two clusterings) is portable, "
        "but it consumes a FlowSOM clustering this package cannot produce, so porting it "
        "would be code with no caller.",
    "cytobatchflagr.flag_resolution":
        "Tier B; requires two CytoBatchFlagR runs, which need a working upstream checkout",
}

#: Metrics this plan named that do not exist upstream at the pinned commit.
UPSTREAM_METRICS_NOT_FOUND: dict[str, str] = {
    "openproblems.cms": (
        "the CellMixS-based 'cms' metric appears in the task CHANGELOG but "
        "src/metrics/cms is absent at the pinned commit; treat it as unavailable rather "
        "than as a metric to reconcile against"
    ),
}
