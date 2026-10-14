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

#: Metrics named in the specification that this build does not implement.  They are
#: listed so that the report can say "unavailable" instead of silently omitting them.
DEFERRED_METRICS: dict[str, str] = {
    "openproblems.iLisi": "Tier B; local-structure metric, not implemented",
    "openproblems.cLisi": "Tier B; local-structure metric, not implemented",
    "openproblems.ratio_consistent_peaks": "Tier B; peak-consistency metric, not implemented",
    "openproblems.average_batch_r2_ct": "Tier B; residual batch association, not implemented",
    "openproblems.functional_marker_preservation_wilcoxon":
        "Tier B; biological-preservation compatibility component, not implemented",
    "openproblems.functional_marker_preservation_cohens_d":
        "Tier B; biological-preservation compatibility component, not implemented",
    "openproblems.emd_mean_ct_horiz": "Tier B; upstream EMD component, not implemented",
    "openproblems.emd_mean_ct_vert": "Tier B; upstream EMD component, not implemented",
    "openproblems.flowsom_mean_mapping_similarity":
        "Tier B; local-structure metric, not implemented",
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
