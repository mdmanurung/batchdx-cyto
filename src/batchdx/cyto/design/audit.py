"""Sample-level design audit, identifiability, and the method eligibility manifest.

The audit runs on samples, never on cells: cells from one donor are not independent
replicates, and every statement here is about which *contrasts* the design can
support.  Perfect confounding disables recommendations for the affected contrast
while leaving descriptive reporting intact.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from batchdx.cyto.schemas.contracts import (
    ComparisonAssessment,
    DesignAudit,
    MethodEligibility,
    SampleRecord,
)
from batchdx.cyto.schemas.enums import ComparisonSupport, ReferenceRole

#: Reference roles that count as *measured technical* evidence.  A biological
#: pseudo-control never counts, whatever the configuration says.
MEASURED_REFERENCE_ROLES = {ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR}


@dataclass(frozen=True)
class EligibilityRule:
    """A method constraint with a stated origin, so no threshold is a hidden constant."""

    name: str
    kind: str
    origin: str
    min_batches: int | None = None
    max_batches: int | None = None
    requires_reference_in_every_batch: bool = False
    requires_identifiable_batch: bool = False


#: Constraints established by reading the pinned upstream sources in Stage 0.
ELIGIBILITY_RULES: dict[str, EligibilityRule] = {
    "cycombine_cluster": EligibilityRule(
        name="cycombine_cluster",
        kind="detector",
        min_batches=3,
        max_batches=30,
        origin=(
            "cyCombine::detect_batch_effect computes EMD and cluster-percentage outlier "
            "analysis only when 3 <= n_batches <= 30 (R/detect_batch_effect.R); outside "
            "that range it emits UMAPs alone."
        ),
    ),
    "cycombine_express": EligibilityRule(
        name="cycombine_express",
        kind="detector",
        min_batches=2,
        origin=(
            "cyCombine::detect_batch_effect_express stops with an error below 2 batches "
            "(R/detect_batch_effect.R)."
        ),
    ),
    "cytobatchflagr": EligibilityRule(
        name="cytobatchflagr",
        kind="detector",
        min_batches=4,
        requires_reference_in_every_batch=True,
        origin=(
            "CytoBatchFlagR flags a batch when a per-control per-marker value falls beyond "
            "1.5 IQR of the quartiles taken across batches (detect_outlier in "
            "R/CytoBatchFlagR_functions.R). The rule requires a control sample in every "
            "batch, and for n <= 3 values the Tukey fence can never be exceeded, so the "
            "assessment is degenerate below 4 batches."
        ),
    ),
    "acquisition_drift": EligibilityRule(
        name="acquisition_drift",
        kind="detector",
        min_batches=1,
        origin="within-file diagnostic; needs no batch structure at all",
    ),
    "composition": EligibilityRule(
        name="composition",
        kind="detector",
        min_batches=2,
        origin="frequency comparison needs at least two batches",
    ),
    "none": EligibilityRule(
        name="none",
        kind="backend",
        min_batches=1,
        origin="the identity candidate is always eligible and is always benchmarked",
    ),
    "cycombine": EligibilityRule(
        name="cycombine",
        kind="backend",
        min_batches=2,
        requires_identifiable_batch=True,
        origin=(
            "cyCombine::correct_data runs ComBat per SOM node; the batch term must be "
            "estimable given the covariates supplied (cyCombine:::check_confound)."
        ),
    ),
    "cytonorm": EligibilityRule(
        name="cytonorm",
        kind="backend",
        min_batches=2,
        requires_reference_in_every_batch=True,
        origin=(
            "CytoNorm.train learns per-batch quantile splines from control files; the "
            "package DESCRIPTION states a control sample is taken along in each batch. "
            "This is the strongest evidence tier batchdx offers for CytoNorm: the goal "
            "quantiles come from measured references present in every batch."
        ),
    ),
    "cytonorm_control_free": EligibilityRule(
        name="cytonorm_control_free",
        kind="backend",
        min_batches=2,
        origin=(
            "CytoNorm 2.0 supports control-free normalisation: QuantileNorm.train "
            "accepts goal='mean' (the default), which learns the goal quantile "
            "distribution as the mean over all batches, or a numeric goal distribution "
            "supplied by the analyst (R/QuantileNorm.R, goal parameter; pinned 2.0.12). "
            "No measured reference is required for that route, so the evidence is "
            "weaker than the reference-based 'cytonorm' tier and the two are reported "
            "separately."
        ),
    ),
}


def _levels(samples: list[SampleRecord], factor: str) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for s in samples:
        out[s.sample_id] = _factor_value(s, factor)
    return out


def _factor_value(sample: SampleRecord, factor: str) -> str | None:
    if hasattr(sample, factor):
        value = getattr(sample, factor)
        if isinstance(value, ReferenceRole):
            return value.value
        return None if value is None else str(value)
    return sample.technical_factors.get(factor)


def build_design_matrix(
    samples: list[SampleRecord], factors: list[str]
) -> tuple[np.ndarray, list[str]]:
    """Treatment-coded design matrix with an intercept, over sample-level factors."""

    columns: list[np.ndarray] = [np.ones(len(samples))]
    names = ["(Intercept)"]
    for factor in factors:
        values = [_factor_value(s, factor) for s in samples]
        levels = sorted({v for v in values if v is not None})
        for level in levels[1:]:  # first level is the reference
            columns.append(np.array([1.0 if v == level else 0.0 for v in values]))
            names.append(f"{factor}::{level}")
    return np.column_stack(columns), names


def _perfectly_confounded(
    samples: list[SampleRecord], batch_factor: str, factor: str
) -> bool:
    """True when every batch contains exactly one level of ``factor``.

    Under that layout the batch term and the factor term span the same subspace, so no
    amount of covariate adjustment separates them.
    """

    by_batch: dict[str, set[str]] = {}
    for s in samples:
        b = _factor_value(s, batch_factor)
        v = _factor_value(s, factor)
        if b is None or v is None:
            continue
        by_batch.setdefault(b, set()).add(v)
    if len(by_batch) < 2:
        return False
    if any(len(v) != 1 for v in by_batch.values()):
        return False
    return len({next(iter(v)) for v in by_batch.values()}) > 1


def reference_components(
    samples: list[SampleRecord], batch_factor: str
) -> tuple[dict[str, list[str]], list[list[str]]]:
    """Measured-reference coverage per batch and the connected components of batches.

    Connectivity is reported because it is real, but no method is assumed able to
    propagate a correction across an indirect connection; each adapter states the
    layout it supports.
    """

    coverage: dict[str, list[str]] = {}
    ref_to_batches: dict[str, set[str]] = {}
    for s in samples:
        batch = _factor_value(s, batch_factor)
        if batch is None:
            continue
        coverage.setdefault(batch, [])
        if s.reference_role in MEASURED_REFERENCE_ROLES and s.reference_id:
            coverage[batch].append(s.reference_id)
            ref_to_batches.setdefault(s.reference_id, set()).add(batch)

    adjacency: dict[str, set[str]] = {b: set() for b in coverage}
    for batches in ref_to_batches.values():
        for a, b in combinations(sorted(batches), 2):
            adjacency[a].add(b)
            adjacency[b].add(a)

    seen: set[str] = set()
    components: list[list[str]] = []
    for node in sorted(adjacency):
        if node in seen:
            continue
        stack, comp = [node], []
        seen.add(node)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in sorted(adjacency[cur]):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        components.append(sorted(comp))
    return {b: sorted(set(v)) for b, v in coverage.items()}, components


def _classify_pair(
    batch_a: str,
    batch_b: str,
    samples: list[SampleRecord],
    batch_factor: str,
    protected: list[str],
    coverage: dict[str, list[str]],
    population: str | None = None,
    population_counts: dict[str, dict[str, int]] | None = None,
) -> ComparisonAssessment:
    in_a = [s for s in samples if _factor_value(s, batch_factor) == batch_a]
    in_b = [s for s in samples if _factor_value(s, batch_factor) == batch_b]

    if population is not None and population_counts is not None:
        in_a = [s for s in in_a if population_counts.get(s.sample_id, {}).get(population, 0) > 0]
        in_b = [s for s in in_b if population_counts.get(s.sample_id, {}).get(population, 0) > 0]

    shared_refs = sorted(set(coverage.get(batch_a, [])) & set(coverage.get(batch_b, [])))
    comparison_id = f"{batch_a}__vs__{batch_b}" + (f"@{population}" if population else "")
    base = dict(
        comparison_id=comparison_id, batch_a=batch_a, batch_b=batch_b, population=population,
        n_samples_a=len(in_a), n_samples_b=len(in_b),
        n_shared_reference_ids=len(shared_refs), shared_reference_ids=shared_refs,
    )

    if not in_a or not in_b:
        return ComparisonAssessment(
            **base, support=ComparisonSupport.DESCRIPTIVE_ONLY,
            reason=(
                f"population {population!r} has no cells in at least one batch of the pair; "
                "the contrast is not supported within this population"
                if population else "one batch of the pair contains no samples"
            ),
        )

    overlapping: list[str] = []
    disjoint: list[str] = []
    for factor in protected:
        levels_a = {_factor_value(s, factor) for s in in_a} - {None}
        levels_b = {_factor_value(s, factor) for s in in_b} - {None}
        if levels_a & levels_b:
            overlapping.extend(f"{factor}={v}" for v in sorted(levels_a & levels_b))
        elif levels_a and levels_b:
            disjoint.append(factor)

    # A shared reference makes the *technical* contrast measurable, but never rescues a
    # protected factor that takes disjoint levels across the two batches: the batch
    # contrast and that biological contrast remain the same contrast.  Identifiability
    # is checked first so a confounded design can never present as reference-supported.
    if disjoint:
        if shared_refs:
            return ComparisonAssessment(
                **base, support=ComparisonSupport.NON_IDENTIFIABLE,
                reason=(
                    f"measured reference(s) {shared_refs} are present in both batches, but "
                    f"protected factor(s) {disjoint} take disjoint levels across them; the "
                    "batch contrast and the biological contrast are not separately "
                    "estimable, and the shared reference does not change that. Adding the "
                    "same factor as a covariate does not resolve this."
                ),
                overlapping_biology=overlapping,
            )
        return ComparisonAssessment(
            **base, support=ComparisonSupport.NON_IDENTIFIABLE,
            reason=(
                f"protected factor(s) {disjoint} take disjoint levels across the two batches; "
                "the batch contrast and the biological contrast are not separately estimable. "
                "Adding the same factor as a covariate does not resolve this."
            ),
            overlapping_biology=overlapping,
        )
    if shared_refs:
        return ComparisonAssessment(
            **base, support=ComparisonSupport.REFERENCE_SUPPORTED,
            reason=f"measured reference identity/identities {shared_refs} present in both batches",
            overlapping_biology=overlapping,
        )
    if overlapping:
        return ComparisonAssessment(
            **base, support=ComparisonSupport.BIOLOGICALLY_STRATIFIED_CONDITIONAL,
            reason=(
                "no measured reference is shared; evidence is conditional on comparable "
                f"biological samples ({', '.join(overlapping)}) behaving identically"
            ),
            overlapping_biology=overlapping,
        )
    return ComparisonAssessment(
        **base, support=ComparisonSupport.DESCRIPTIVE_ONLY,
        reason="no measured reference and no declared protected factor to stratify on",
    )


def _evaluate_rule(
    rule: EligibilityRule,
    n_batches: int,
    batches_with_reference: int,
    batch_identifiable: bool,
) -> MethodEligibility:
    requirements: dict[str, bool] = {}
    failures: list[tuple[str, str]] = []

    if rule.min_batches is not None:
        ok = n_batches >= rule.min_batches
        requirements[f"n_batches >= {rule.min_batches}"] = ok
        if not ok:
            failures.append((
                "batch_count_below_minimum",
                f"{n_batches} batch(es) present, {rule.min_batches} required",
            ))
    if rule.max_batches is not None:
        ok = n_batches <= rule.max_batches
        requirements[f"n_batches <= {rule.max_batches}"] = ok
        if not ok:
            failures.append((
                "batch_count_above_maximum",
                f"{n_batches} batch(es) present, at most {rule.max_batches} supported",
            ))
    if rule.requires_reference_in_every_batch:
        ok = n_batches > 0 and batches_with_reference == n_batches
        requirements["measured reference in every batch"] = ok
        if not ok:
            failures.append((
                "missing_reference_in_some_batch",
                f"{batches_with_reference} of {n_batches} batch(es) contain a measured "
                "reference sample",
            ))
    if rule.requires_identifiable_batch:
        requirements["batch effect estimable given protected factors"] = batch_identifiable
        if not batch_identifiable:
            failures.append((
                "batch_not_identifiable",
                "at least one protected factor is perfectly confounded with the batch factor",
            ))

    if failures:
        # Every failed requirement is reported: a method can fail several at once
        # (e.g. too few batches AND a missing reference), and a lossy reason would
        # hide what else must be fixed.
        codes = [code for code, _ in failures]
        details = "; ".join(detail for _, detail in failures)
        return MethodEligibility(
            name=rule.name, kind=rule.kind, eligible=False,
            reason_code="+".join(codes),
            reason=f"{details}. Constraint origin: {rule.origin}",
            requirements_met=requirements,
            override_allowed="batch_not_identifiable" not in codes,
        )
    return MethodEligibility(
        name=rule.name, kind=rule.kind, eligible=True, reason_code="eligible",
        reason=f"all declared requirements met. Constraint origin: {rule.origin}",
        requirements_met=requirements,
    )


def audit_design(
    samples: list[SampleRecord],
    *,
    dataset_id: str,
    correction_factor: str = "batch",
    nested_within: str | None = None,
    protected: list[str] | None = None,
    donor_factor: str | None = None,
    minimum_batches: int = 3,
    minimum_batches_origin: str = "configuration default",
    require_shared_reference: bool = True,
    populations: list[str] | None = None,
    population_counts: dict[str, dict[str, int]] | None = None,
) -> DesignAudit:
    """Audit the design and emit the eligibility manifest.

    ``population_counts`` maps ``sample_id`` to per-population cell counts.  When it is
    supplied the pairwise classification is repeated within each population, because a
    globally valid design can be unsupported for a rare cell type.
    """

    protected = list(protected or [])
    batch_of = _levels(samples, correction_factor)
    batches = sorted({b for b in batch_of.values() if b is not None})
    n_batches = len(batches)

    warnings: list[str] = []
    missing_batch = [sid for sid, b in batch_of.items() if b is None]
    if missing_batch:
        warnings.append(
            f"{len(missing_batch)} sample(s) have no value for the correction factor "
            f"{correction_factor!r} and are excluded from every batch contrast"
        )
    if donor_factor and all(_factor_value(s, donor_factor) is None for s in samples):
        warnings.append(
            f"no sample declares {donor_factor!r}; descriptive diagnosis remains available "
            "but donor-matched metrics and donor-level uncertainty are disabled"
        )

    contingency: dict[str, dict[str, int]] = {}
    for factor in protected:
        table: dict[str, dict[str, int]] = {}
        for s in samples:
            b = _factor_value(s, correction_factor)
            v = _factor_value(s, factor)
            if b is None:
                continue
            table.setdefault(b, {}).setdefault(str(v), 0)
            table[b][str(v)] += 1
        contingency[factor] = {b: dict(v) for b, v in table.items()}

    residual: dict[str, int] = {}
    for s in samples:
        b = _factor_value(s, correction_factor)
        if b is None:
            continue
        key = b if not protected else (
            b + "|" + "|".join(f"{f}={_factor_value(s, f)}" for f in protected)
        )
        residual[key] = residual.get(key, 0) + 1

    factors = [correction_factor] + protected + ([nested_within] if nested_within else [])
    matrix, names = build_design_matrix(samples, factors)
    rank = int(np.linalg.matrix_rank(matrix)) if matrix.size else 0
    full_rank = rank == matrix.shape[1] if matrix.size else False

    confounded = [f for f in protected if _perfectly_confounded(samples, correction_factor, f)]
    batch_identifiable = not confounded

    coverage, components = reference_components(samples, correction_factor)
    batches_with_reference = sum(1 for b in batches if coverage.get(b))

    comparisons: list[ComparisonAssessment] = []
    for a, b in combinations(batches, 2):
        comparisons.append(
            _classify_pair(a, b, samples, correction_factor, protected, coverage)
        )
        for pop in populations or []:
            comparisons.append(
                _classify_pair(a, b, samples, correction_factor, protected, coverage,
                               population=pop, population_counts=population_counts)
            )

    below_floor = False
    floor_reasons: list[str] = []
    if n_batches < minimum_batches:
        below_floor = True
        floor_reasons.append(
            f"{n_batches} batch(es) present; the declared minimum is {minimum_batches}"
        )
    if require_shared_reference and not any(
        c.support is ComparisonSupport.REFERENCE_SUPPORTED
        for c in comparisons if c.population is None
    ):
        below_floor = True
        floor_reasons.append(
            "no measured reference sample is shared across any pair of the batches under "
            "comparison"
        )

    eligibility = [
        _evaluate_rule(rule, n_batches, batches_with_reference, batch_identifiable)
        for rule in ELIGIBILITY_RULES.values()
    ]

    if confounded:
        warnings.append(
            f"protected factor(s) {confounded} are perfectly confounded with "
            f"{correction_factor!r}; no automatic correction recommendation can be made for "
            "the affected contrast, and adding the factor as a covariate does not help"
        )
    if len(components) > 1:
        warnings.append(
            f"the measured-reference graph has {len(components)} disconnected component(s) "
            f"{components}; correction cannot be assumed to propagate between them"
        )

    return DesignAudit(
        dataset_id=dataset_id,
        correction_factor=correction_factor,
        nested_within=nested_within,
        protected_factors=protected,
        batches=batches,
        n_batches=n_batches,
        contingency=contingency,
        design_rank=rank,
        design_full_rank=full_rank,
        confounded_factors=confounded,
        residual_replication=residual,
        reference_coverage=coverage,
        reference_graph_components=components,
        below_design_floor=below_floor,
        design_floor_reason="; ".join(floor_reasons) or None,
        design_floor_origin=minimum_batches_origin,
        comparisons=comparisons,
        eligibility=eligibility,
        warnings=warnings + ([f"design matrix columns: {names}"] if not full_rank else []),
    )
