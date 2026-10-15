"""Train / tune / validate splits held out at the biospecimen level.

A fit-and-apply correction is fitted on one set of samples and then applied to
another.  If the two sets share a biospecimen, the model has already seen the
biology it is being scored on, and the held-out number stops being a held-out
number.  The partition here is therefore over *biospecimens*, never over samples
and never over cells: every sample of a biospecimen lands entirely in one split.

Four properties of this module are deliberate, and each is covered by a test:

* **The unit is declared, never chosen in silence.**  The manifest records which
  identity level the split is keyed on: ``donor_id`` when the design declares a
  donor factor, otherwise ``biospecimen_id``, otherwise ``sample_id``.  The last
  is an honest fallback and a weak one — with no biospecimen identity in the
  sample manifest there is nothing to hold out but the tube — and the manifest
  says so instead of presenting it as donor-level independence.
* **The assignment cannot depend on arrival order.**  Strata are sorted and the
  only randomness is a permutation of an already-sorted list.  Re-running, reading
  the samples in another order, or changing the candidate list cannot move a donor
  between splits.
* **A split that cannot hold both arms is refused.**  A split missing one arm of
  the contrast is not a weaker result, it is a meaningless one, and a summary
  computed over it would look like evidence.  :func:`build_splits` raises instead
  of returning one.
* **Split sizes do not depend on the seed; only membership does.**  Units are
  dealt to whichever split is furthest below its target share, so a stratum holding
  fewer biospecimens than there are splits still reaches the held-out split, and
  the seed decides which donor lands where rather than how many.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, Final

import numpy as np

from batchdx.cyto.schemas.contracts import ComparisonAssessment, SampleRecord
from batchdx.cyto.schemas.enums import ReferenceRole

#: The three splits, in the order every tie is broken.  A fixed order rather than a
#: set, because "which split wins a tie" must be a recorded fact and not an
#: artefact of set iteration.
SPLIT_NAMES: Final[tuple[str, str, str]] = ("train", "tune", "validate")

#: The only split whose metrics may support a recommendation.  A metric computed on
#: data the model was fitted on describes the fit; it is not evidence that the
#: correction generalises.
RECOMMENDATION_SPLIT: Final[str] = "validate"

#: Reserved key inside ``SamplingManifest.splits`` recording which identity level
#: the split is keyed on.  The field is ``dict[str, list[str]]`` and cannot carry a
#: nested metadata block without a contract change, so the unit is recorded as a
#: one-element list under a name that cannot collide with a split name.  Read it
#: through :func:`split_unit_field`, never by indexing ``splits`` directly.
SPLIT_UNIT_KEY: Final[str] = "split_unit"

#: Identity levels, most to least independent, in the order the brief fixes them.
BIOSPECIMEN_UNIT_ORDER: Final[tuple[str, str, str]] = (
    "donor_id", "biospecimen_id", "sample_id",
)

#: Target shares.  60/20/20 is a convention, not a result, and it is recorded with
#: its origin so a reader can disagree with it explicitly rather than inherit it.
#: Every share must be strictly positive: a split with no members is not held out,
#: so this build refuses to emit one instead of reporting an empty partition.
DEFAULT_SPLIT_PROPORTIONS: Final[dict[str, float]] = {
    "train": 0.6, "tune": 0.2, "validate": 0.2,
}
SPLIT_PROPORTION_ORIGIN: Final[str] = (
    "uncalibrated engineering default (60/20/20). It is a convention for reserving a "
    "tuning set and a held-out set, not a value derived from any property of the data; "
    "override it at the call site and record the override in the manifest"
)

#: Reference roles that are *measured* technical evidence.  A biological pseudo
#: control is deliberately absent: it is a sample, not a technical replicate, and
#: the design audit and docs/usage.md both forbid counting it as one.
_MEASURED_REFERENCE_ROLES: Final[frozenset] = frozenset({
    ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR,
})


class SplitError(ValueError):
    """Raised when an honest split cannot be formed, or cannot be read back.

    Refusing is the point.  A split that leaks a biospecimen, or that holds one arm
    of the contrast and not the other, is worse than no split: it is reported, and
    something downstream will treat it as a held-out result.
    """


# --------------------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Sample:
    """Per-sample facts the split needs, with the identity fields already resolved."""

    sample_id: str
    unit_value: str
    group: str
    batch: str
    is_measured_reference: bool


@dataclass(frozen=True)
class Biospecimen:
    """One partition key and every sample that carries it."""

    unit_value: str
    samples: list[_Sample]

    @property
    def is_measured_reference(self) -> bool:
        return all(s.is_measured_reference for s in self.samples)

    @property
    def stratum(self) -> tuple[bool, tuple[str, ...], tuple[str, ...]]:
        """Crossed stratification key: technical role, biological groups, batches.

        Defined on the *biospecimen*, not the sample, so a donor that appears in
        three batches is one unit in one stratum rather than three copies of
        itself.  Stratifying on the cross of group and batch rather than on group
        alone is what keeps every split able to see every batch: splitting on group
        alone can place all B2 donors in one split, and a contrast estimated inside
        a single batch is not a batch contrast.
        """

        return (
            self.is_measured_reference,
            tuple(sorted({s.group for s in self.samples if s.group})),
            tuple(sorted({s.batch for s in self.samples if s.batch})),
        )

    @property
    def arms(self) -> list[str]:
        """Biological group levels this biospecimen contributes."""

        if self.is_measured_reference:
            return []
        return sorted({s.group for s in self.samples if s.group})


@dataclass(frozen=True)
class SplitCoverage:
    """What one split actually contains, counted rather than asserted."""

    split: str
    units: list[str]
    samples: list[str]
    arm_counts: dict[str, int] = field(default_factory=dict)
    batch_counts: dict[str, int] = field(default_factory=dict)
    reference_batches: list[str] = field(default_factory=list)

    @property
    def n_units(self) -> int:
        return len(self.units)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def summary(self) -> str:
        arms = ", ".join(f"{a}={n}" for a, n in sorted(self.arm_counts.items()))
        refs = sorted(self.reference_batches)
        return (
            f"{self.split}: {self.n_units} biospecimen(s), {self.n_samples} sample(s); "
            f"arms {arms or 'none'}; batches {dict(sorted(self.batch_counts.items())) or 'none'}"
            f"; measured reference in batch(es) {refs or 'none'}"
        )


@dataclass(frozen=True)
class SplitAssignment:
    """A complete, checked partition of samples over biospecimen-level splits."""

    unit_field: str
    unit_reason: str
    assignments: dict[str, str]           # biospecimen value -> split
    coverage: dict[str, SplitCoverage]    # split -> what it holds
    arms: list[str]                       # biological group levels present overall
    strata: dict[str, int]                # stratum key -> number of biospecimens
    proportions: dict[str, float]
    seed: int
    warnings: list[str] = field(default_factory=list)

    def split_of_unit(self, unit_value: str) -> str:
        return self.assignments[unit_value]

    def samples_in(self, split: str) -> list[str]:
        return self.coverage[split].samples

    def units_in(self, split: str) -> list[str]:
        return self.coverage[split].units

    def split_of_sample(self, sample_id: str) -> str:
        for split, cov in self.coverage.items():
            if sample_id in cov.samples:
                return split
        raise KeyError(f"sample {sample_id!r} is in no split")

    def to_manifest(self) -> dict[str, list[str]]:
        """The value written into ``SamplingManifest.splits``."""

        out: dict[str, list[str]] = {SPLIT_UNIT_KEY: [self.unit_field]}
        for name in SPLIT_NAMES:
            out[name] = list(self.coverage[name].samples)
        return out

    def describe(self) -> str:
        return describe(self)


@dataclass(frozen=True)
class ComparisonInSplits:
    """Where a batch contrast can be estimated, and whether it is held out."""

    comparison_id: str
    batch_a: str
    batch_b: str
    evaluable_in: list[str]
    held_out: bool
    reason: str

    @property
    def can_support_recommendation(self) -> bool:
        return self.held_out


# --------------------------------------------------------------------------------------
# Unit resolution
# --------------------------------------------------------------------------------------


def _field_value(record: SampleRecord, field_name: str) -> str | None:
    if field_name == "sample_id":
        raw: Any = record.sample_id
    else:
        raw = getattr(record, field_name, None)
    if raw is None:
        return None
    text = str(raw).strip()
    return text or None


def resolve_biospecimen_unit(
    samples: Mapping[str, SampleRecord], donor_factor: str | None = None
) -> tuple[str, str]:
    """Return ``(unit_field, reason)``, preferring the most independent identity.

    A donor factor declared by the design wins, but only if every sample carries
    one.  A factor populated for some samples and missing for others is not a
    usable partition key, and mixing the two would silently put one donor's samples
    under two different units — precisely the leak this module exists to prevent.
    ``biospecimen_id`` is the next most independent declared identity, and
    ``sample_id`` is the last resort.
    """

    if not samples:
        raise SplitError("no samples to split; the sample manifest is empty")

    if donor_factor:
        missing = sorted(s for s, r in samples.items() if _field_value(r, donor_factor) is None)
        if not missing:
            return donor_factor, (
                f"the design declares a donor factor ({donor_factor!r}) and every sample "
                f"carries one, so the split is keyed on {donor_factor!r}, the most "
                "independent biospecimen identity available"
            )
        note = (
            f"{donor_factor!r} is declared as the donor factor but is absent for "
            f"{len(missing)} sample(s) ({missing[:3]}{'...' if len(missing) > 3 else ''}), so "
            "it cannot key the partition without leaking one donor's samples across splits"
        )
    else:
        note = "the design declares no donor factor, so the split cannot be keyed on a donor"

    missing_bio = sorted(
        s for s, r in samples.items() if _field_value(r, "biospecimen_id") is None
    )
    if not missing_bio:
        return "biospecimen_id", f"{note}; every sample declares biospecimen_id, so it is used"
    return "sample_id", (
        f"{note}; biospecimen_id is absent for {len(missing_bio)} sample(s) "
        f"({missing_bio[:3]}{'...' if len(missing_bio) > 3 else ''}), so the split falls "
        "back to sample_id"
    )


def _biospecimens(
    samples: Mapping[str, SampleRecord], unit_field: str
) -> dict[str, Biospecimen]:
    grouped: dict[str, list[_Sample]] = {}
    for sid in sorted(samples):
        rec = samples[sid]
        value = _field_value(rec, unit_field)
        if value is None:
            raise SplitError(
                f"sample {sid!r} has no value for the split unit {unit_field!r}; a sample "
                "without an identity cannot be assigned to a split"
            )
        grouped.setdefault(value, []).append(_Sample(
            sample_id=sid,
            unit_value=value,
            group=str(rec.group) if rec.group is not None else "",
            batch=str(rec.batch) if rec.batch is not None else "",
            is_measured_reference=rec.reference_role in _MEASURED_REFERENCE_ROLES,
        ))
    return {k: Biospecimen(unit_value=k, samples=v) for k, v in sorted(grouped.items())}


# --------------------------------------------------------------------------------------
# The partition
# --------------------------------------------------------------------------------------


def _normalise_proportions(proportions: Mapping[str, float] | None) -> dict[str, float]:
    raw = dict(DEFAULT_SPLIT_PROPORTIONS if proportions is None else proportions)
    unknown = sorted(set(raw) - set(SPLIT_NAMES))
    missing = sorted(set(SPLIT_NAMES) - set(raw))
    if unknown or missing:
        raise SplitError(
            f"split proportions must name exactly {list(SPLIT_NAMES)}; "
            f"unknown {unknown}, missing {missing}"
        )
    values = {k: float(v) for k, v in raw.items()}
    non_positive = sorted(k for k, v in values.items() if not v > 0)
    if non_positive:
        raise SplitError(
            f"split proportion(s) {non_positive} are not strictly positive. A split with no "
            "members is not held out, so this build refuses to emit one rather than "
            "reporting an empty split as a partition"
        )
    total = sum(values.values())
    if total <= 0:
        raise SplitError("split proportions sum to zero")
    return {k: v / total for k, v in values.items()}


def _deal(
    strata: dict[tuple[bool, tuple[str, ...], tuple[str, ...]], list[Biospecimen]],
    shares: dict[str, float],
    seed: int,
) -> dict[str, str]:
    """Assign each biospecimen to a split, deterministically from ``seed``.

    Biospecimens are dealt one at a time to whichever split is furthest below its
    target share, ties broken by the fixed :data:`SPLIT_NAMES` order.  Dealing
    rather than apportioning each stratum separately is what lets a stratum with
    fewer biospecimens than there be splits still reach the held-out split: a
    largest-remainder apportionment of two units over 60/20/20 yields 1/0/0, so
    ``validate`` would be empty on an entirely ordinary design.

    Because the deficit depends only on how many units have already been dealt and
    not on which ones, the *sizes* this function produces are seed-independent and
    the seed decides only which donor goes where.  The arm repair that follows can
    still move units, and does, so the final sizes need not match these exactly; the
    moves it makes are recorded as warnings.

    The randomness is a permutation of an already-sorted list, so the result
    depends on the seed and on the content of each stratum and on nothing else.
    """

    rng = np.random.default_rng(seed)
    filled = dict.fromkeys(SPLIT_NAMES, 0)
    out: dict[str, str] = {}

    for key in sorted(strata):
        units = sorted(strata[key], key=lambda b: b.unit_value)
        for i in rng.permutation(len(units)):
            total = sum(filled.values())
            best, best_deficit = SPLIT_NAMES[0], None
            for name in SPLIT_NAMES:
                share = filled[name] / total if total else 0.0
                deficit = shares[name] - share
                # Strict ``>`` over the fixed name order makes ties resolve to the
                # first name rather than to an arbitrary one.
                if best_deficit is None or deficit > best_deficit:
                    best, best_deficit = name, deficit
            filled[best] += 1
            out[units[i].unit_value] = best
    return out


def _repair_arms(
    biospecimens: dict[str, Biospecimen],
    assignments: dict[str, str],
    arms: list[str],
) -> tuple[dict[str, str], list[str]]:
    """Move whole biospecimens until every split holds every biological arm.

    Dealing by target share alone can leave a split without an arm — and then, on a
    design where biology is confounded with batch, the only splits that can hold the
    contrast are the ones that were given the arm.  That is not a weak partition, it
    is an unusable one, and a summary computed over it would still be reported.

    The repair moves a whole biospecimen (never a sample) from the split holding the
    *most* of the missing arm, and only when that split would keep at least one.
    Each pass therefore fills at least one (split, arm) cell and can never empty a
    filled one, so the loop terminates in at most ``len(arms) * len(SPLIT_NAMES)``
    passes.  If it stops without filling every cell, the arms are not jointly
    realisable and the caller refuses.

    Returns the repaired assignment and a list of the moves made, so the adjustment
    is a recorded fact rather than a silent one.
    """

    counts: dict[str, dict[str, int]] = {name: dict.fromkeys(arms, 0) for name in SPLIT_NAMES}
    arms_of: dict[str, list[str]] = {}
    for unit, bio in biospecimens.items():
        arms_of[unit] = bio.arms
        for a in bio.arms:
            counts[assignments[unit]][a] += 1

    order = {name: i for i, name in enumerate(SPLIT_NAMES)}
    current = dict(assignments)
    moves: list[str] = []

    for _ in range(len(arms) * len(SPLIT_NAMES) + 1):
        missing = [
            (name, a)
            for name in SPLIT_NAMES for a in arms
            if counts[name][a] == 0
        ]
        if not missing:
            return current, moves
        progressed = False
        for name, arm in missing:
            donors = [d for d in SPLIT_NAMES if d != name and counts[d][arm] > 1]
            if not donors:
                continue
            source = max(donors, key=lambda d: (counts[d][arm], -order[d]))
            candidates = sorted(
                u for u, a_list in arms_of.items()
                if arm in a_list and current[u] == source
            )
            if not candidates:  # pragma: no cover - implied by counts[source][arm] > 1
                continue
            # Prefer the unit whose batches are already represented in the split it
            # leaves, so the repair damages batch coverage as little as possible.
            source_batches = {
                s.batch for bio in biospecimens.values()
                if current[bio.unit_value] == source for s in bio.samples
            }
            move = min(candidates, key=lambda u: (
                not {s.batch for s in biospecimens[u].samples} <= source_batches,
                len(biospecimens[u].samples),
                u,
            ))
            current[move] = name
            for a in arms_of[move]:
                counts[source][a] -= 1
                counts[name][a] += 1
            moves.append(
                f"biospecimen {move!r} moved from split {source!r} to {name!r} so that "
                f"split {name!r} holds arm {arm!r}"
            )
            progressed = True
        if not progressed:
            break
    return current, moves


def build_splits(
    samples: Mapping[str, SampleRecord],
    *,
    donor_factor: str | None = None,
    seed: int,
    proportions: Mapping[str, float] | None = None,
    unit_field: str | None = None,
    require_all_arms: bool = True,
) -> SplitAssignment:
    """Partition samples into ``train`` / ``tune`` / ``validate`` by biospecimen.

    With ``require_all_arms`` (the default) this raises unless every split holds every
    biological arm, because a split that cannot be compared is not a weak partition but
    an unusable one.

    ``require_all_arms=False`` is for callers that must produce *something* whatever the
    design — the sampling manifest, which also serves detectors, metrics and the report
    and has no business refusing to run because one arm has fewer biospecimens than
    there are splits.  It never hides the deficiency: every missing (split, arm) cell
    becomes a warning naming the split, the arm, and the reason, and
    :func:`missing_arms` lets a consumer re-check the frozen manifest.  A fit-and-apply
    backend must use :func:`splits_hold_every_arm`, not this flag, because it may not
    fit on a partition that cannot be compared.
    """

    shares = _normalise_proportions(proportions)
    if unit_field is None:
        unit_field, reason = resolve_biospecimen_unit(samples, donor_factor)
    else:
        reason = f"the split unit was pinned to {unit_field!r} by the caller"

    biospecimens = _biospecimens(samples, unit_field)
    warnings: list[str] = []

    if unit_field == "sample_id":
        warnings.append(
            "the split unit is sample_id: no donor factor and no biospecimen_id is "
            "declared for every sample, so the only identity that can be held out is the "
            "tube. If several tubes come from one donor, that donor's data spans more "
            "than one split and this partition does not provide donor-level "
            "independence. Declare design.donor, or a biospecimen_id column in "
            "samples.tsv, to fix it"
        )

    strata: dict[tuple[bool, tuple[str, ...], tuple[str, ...]], list[Biospecimen]] = {}
    for bio in biospecimens.values():
        strata.setdefault(bio.stratum, []).append(bio)

    assignments = _deal(strata, shares, seed)

    arms = sorted({a for bio in biospecimens.values() for a in bio.arms})
    assignments, repairs = _repair_arms(biospecimens, assignments, arms)

    unassigned = sorted(set(biospecimens) - set(assignments))
    if unassigned:
        raise SplitError(  # pragma: no cover - defensive
            f"{len(unassigned)} biospecimen(s) were not assigned to any split: "
            f"{unassigned[:5]}"
        )

    # The leak check, stated as an invariant rather than assumed from the algorithm.
    by_split: dict[str, list[str]] = {name: [] for name in SPLIT_NAMES}
    unit_split: dict[str, str] = {}
    for bio in biospecimens.values():
        split = assignments[bio.unit_value]
        if split not in by_split:
            raise SplitError(f"internal assignment produced unknown split {split!r}")
        unit_split[bio.unit_value] = split
        by_split[split].extend(s.sample_id for s in bio.samples)
    for name in SPLIT_NAMES:
        by_split[name] = sorted(by_split[name])
    split_of = unit_split

    coverage: dict[str, SplitCoverage] = {}
    for name in SPLIT_NAMES:
        sids = by_split[name]
        if not sids:
            raise SplitError(
                f"split {name!r} ended up with no samples; the design has too few "
                f"{unit_field} values to fill {len(SPLIT_NAMES)} splits. Declare more "
                "biospecimens rather than accepting an empty split"
            )
        members = [s for bio in biospecimens.values()
                   if split_of[bio.unit_value] == name for s in bio.samples]
        arm_counts: dict[str, int] = {}
        for s in members:
            if not s.is_measured_reference and s.group:
                arm_counts[s.group] = arm_counts.get(s.group, 0) + 1
        batch_counts: dict[str, int] = {}
        for s in members:
            if s.batch:
                batch_counts[s.batch] = batch_counts.get(s.batch, 0) + 1
        coverage[name] = SplitCoverage(
            split=name,
            units=sorted(b.unit_value for b in biospecimens.values()
                         if split_of[b.unit_value] == name),
            samples=sids,
            arm_counts=dict(sorted(arm_counts.items())),
            batch_counts=dict(sorted(batch_counts.items())),
            reference_batches=sorted({s.batch for s in members
                                      if s.is_measured_reference and s.batch}),
        )

    if arms:
        gaps: list[tuple[str, list[str]]] = []
        for name in SPLIT_NAMES:
            absent = [a for a in arms if coverage[name].arm_counts.get(a, 0) == 0]
            if absent:
                gaps.append((name, absent))
        if gaps:
            detail = "; ".join(f"{name} lacks {absent}" for name, absent in gaps)
            per_arm = ", ".join(
                f"{a}: {sum(1 for bio in biospecimens.values() if a in bio.arms)} "
                f"biospecimen(s)"
                for a in arms
            )
            reason = (
                f"split(s) {detail}. This design holds {per_arm}, and an arm with fewer "
                f"biospecimens than there are splits cannot be spread across all of them"
            )
            if require_all_arms:
                raise SplitError(
                    f"the partition cannot hold every arm of the biological contrast: "
                    f"{reason}. Every split would be scored on a contrast that exists only "
                    "in the other splits, so no partition is emitted. Add biospecimens to "
                    "the thinner arm rather than accepting a split that cannot be compared"
                )
            warnings.append(
                f"the train/tune/validate split is NOT arm-balanced: {reason}. A batch or "
                "population contrast computed inside such a split is conditional on the "
                "other splits and must not be read as held-out evidence. A fit-and-apply "
                "backend must refuse such a partition"
            )
    else:
        warnings.append(
            "no biological group level is present outside the measured reference, so "
            "stratification by group is vacuous and this partition cannot support a "
            "biological preservation claim"
        )

    all_batches = sorted({s.batch for bio in biospecimens.values() for s in bio.samples
                          if s.batch})
    has_reference = any(b.is_measured_reference for b in biospecimens.values())
    for name in SPLIT_NAMES:
        absent = [b for b in all_batches if coverage[name].batch_counts.get(b, 0) == 0]
        if absent:
            warnings.append(
                f"split {name!r} contains no sample from batch(es) {absent}. Any batch "
                "contrast involving them is unsupported inside that split, and a contrast "
                "estimated from it is not held-out evidence"
            )
        if has_reference and not coverage[name].reference_batches:
            warnings.append(
                f"split {name!r} contains no measured technical reference, so evidence "
                "computed there is conditional on comparable biological samples behaving "
                "identically rather than reference-supported"
            )

    thin = sorted(str(k) for k, v in strata.items() if len(v) < len(SPLIT_NAMES))
    if thin:
        warnings.append(
            f"{len(thin)} stratification stratum/strata hold fewer biospecimens than there "
            f"are splits ({thin if len(thin) <= 5 else 'several; see SplitAssignment.strata'}"
            "); those biospecimens cannot be spread across every split, which is why "
            "batch and reference coverage is reported per split rather than assumed"
        )

    if repairs:
        warnings.append(
            f"{len(repairs)} biospecimen(s) were moved after allocation so that every split "
            "holds every arm of the biological contrast; the target shares are a preference "
            "and arm coverage is not: " + "; ".join(repairs)
        )

    return SplitAssignment(
        unit_field=unit_field,
        unit_reason=reason,
        assignments=assignments,
        coverage=coverage,
        arms=arms,
        strata={str(k): len(v) for k, v in sorted(strata.items())},
        proportions=shares,
        seed=seed,
        warnings=warnings,
    )


# --------------------------------------------------------------------------------------
# Reading the frozen manifest back
# --------------------------------------------------------------------------------------


def split_members(splits: Mapping[str, Any]) -> dict[str, list[str]]:
    """The split -> sample mapping recorded in the manifest, without the unit key.

    Strict on purpose: a key that is neither a split nor :data:`SPLIT_UNIT_KEY` is an
    error rather than something dropped, so a later writer cannot add metadata here
    and have it disappear at the reader.
    """

    unknown = sorted(set(splits) - set(SPLIT_NAMES) - {SPLIT_UNIT_KEY})
    if unknown:
        raise SplitError(
            f"SamplingManifest.splits carries unknown key(s) {unknown}; expected "
            f"{list(SPLIT_NAMES)} and optionally {SPLIT_UNIT_KEY!r}"
        )
    absent = sorted(set(SPLIT_NAMES) - set(splits))
    if absent:
        raise SplitError(
            f"SamplingManifest.splits is missing split(s) {absent}; a partial partition "
            "cannot be read as a train/tune/validate split"
        )
    out = {name: [str(s) for s in splits[name]] for name in SPLIT_NAMES}
    shared = sorted(set(out["train"]) & set(out["tune"])
                    | set(out["train"]) & set(out["validate"])
                    | set(out["tune"]) & set(out["validate"]))
    if shared:
        raise SplitError(
            f"sample(s) {shared} appear in more than one split; the recorded partition "
            "leaks and is not read as a held-out split"
        )
    return out


def split_unit_field(splits: Mapping[str, Any]) -> str | None:
    """The identity level the recorded split is keyed on, or None if unrecorded."""

    values = splits.get(SPLIT_UNIT_KEY)
    if values is None:
        return None
    if len(values) != 1:
        raise SplitError(
            f"{SPLIT_UNIT_KEY!r} must name exactly one identity level, got {list(values)}"
        )
    return str(values[0])


def has_splits(splits: Mapping[str, Any] | None) -> bool:
    """Whether a manifest carries a usable partition, without raising on an absent one."""

    if not splits:
        return False
    try:
        members = split_members(splits)
    except SplitError:
        return False
    return all(members[name] for name in SPLIT_NAMES)


def missing_arms(
    splits: Mapping[str, Any],
    samples: Mapping[str, SampleRecord],
) -> dict[str, list[str]]:
    """Split -> biological arms that split does not hold, for a frozen manifest.

    Recomputed from the manifest's own sample list rather than read from a stored
    count, so a consumer is checking the artefact it was handed and not a claim
    someone made about it.
    """

    members = split_members(splits)
    arms = sorted({
        str(samples[s].group)
        for s in members["train"] + members["tune"] + members["validate"]
        if s in samples and samples[s].group
        and samples[s].reference_role not in _MEASURED_REFERENCE_ROLES
    })
    out: dict[str, list[str]] = {}
    for name in SPLIT_NAMES:
        held = {
            str(samples[s].group) for s in members[name]
            if s in samples and samples[s].group
            and samples[s].reference_role not in _MEASURED_REFERENCE_ROLES
        }
        absent = [a for a in arms if a not in held]
        if absent:
            out[name] = absent
    return out


def splits_hold_every_arm(
    splits: Mapping[str, Any], samples: Mapping[str, SampleRecord]
) -> bool:
    """Whether the recorded partition can compare both arms inside every split.

    This is the check a fit-and-apply backend must make before fitting on it.  A
    model fitted on a partition whose validation split cannot hold one arm produces a
    held-out number for a contrast that does not exist inside the held-out data.
    """

    return not missing_arms(splits, samples)


# --------------------------------------------------------------------------------------
# The decision-layer hook
# --------------------------------------------------------------------------------------


def split_for_comparison(
    comparison: ComparisonAssessment,
    splits: Mapping[str, Any],
    samples: Mapping[str, SampleRecord],
) -> ComparisonInSplits:
    """Where a batch contrast can be estimated, and whether it is held out.

    A contrast between two batches can only be estimated from a split that holds at
    least one sample of each, because a contrast estimated inside a single batch is
    not a batch contrast.  A contrast estimable from the training split but not from
    the validation split is **not** held out, and presenting it as supporting
    evidence is the failure this function exists to make visible.

    ``decisions/policy.py`` does not consult this yet; wiring it in is a required
    central change.  This is the half that can be written and tested without it.
    """

    members = split_members(splits)
    evaluable: list[str] = []
    for name in SPLIT_NAMES:
        batches = {samples[s].batch for s in members[name] if s in samples}
        if comparison.batch_a in batches and comparison.batch_b in batches:
            evaluable.append(name)

    held_out = RECOMMENDATION_SPLIT in evaluable
    if held_out:
        reason = (
            f"the {RECOMMENDATION_SPLIT!r} split holds samples from both "
            f"{comparison.batch_a} and {comparison.batch_b}, so this contrast can be "
            "estimated on data held out from any fit"
        )
    elif evaluable:
        reason = (
            f"this contrast is estimable only in {evaluable}, never in "
            f"{RECOMMENDATION_SPLIT!r}: it was computed on data a fit could have seen, so "
            "it cannot by itself support a recommendation"
        )
    else:
        reason = (
            f"no split holds a sample of both {comparison.batch_a} and "
            f"{comparison.batch_b}; the contrast is not estimable in this partition at "
            "all, which is a coverage failure rather than a weak result"
        )
    return ComparisonInSplits(
        comparison_id=comparison.comparison_id,
        batch_a=comparison.batch_a,
        batch_b=comparison.batch_b,
        evaluable_in=evaluable,
        held_out=held_out,
        reason=reason,
    )


def comparisons_usable_for_recommendation(
    comparisons: Iterable[ComparisonAssessment],
    splits: Mapping[str, Any],
    samples: Mapping[str, SampleRecord],
) -> list[ComparisonInSplits]:
    """Every comparison in the frozen set, each with its held-out verdict."""

    return [split_for_comparison(c, splits, samples) for c in comparisons]


def held_out_batch_pairs(
    splits: Mapping[str, Any], samples: Mapping[str, SampleRecord]
) -> list[tuple[str, str]]:
    """Batch pairs the ``validate`` split can actually contrast."""

    members = split_members(splits).get(RECOMMENDATION_SPLIT, [])
    batches = sorted({
        str(samples[s].batch) for s in members if s in samples and samples[s].batch
    })
    return [(a, b) for i, a in enumerate(batches) for b in batches[i + 1:]]


def describe(assignment: SplitAssignment) -> str:
    """A short human-readable account, used in warnings and in failure messages."""

    return "\n".join([
        f"split unit: {assignment.unit_field} ({assignment.unit_reason}); target shares "
        f"{assignment.proportions} — {SPLIT_PROPORTION_ORIGIN}",
        *(assignment.coverage[name].summary() for name in SPLIT_NAMES),
    ])
