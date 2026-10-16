"""Evaluation labels, deterministic sampling, and the frozen common partition.

The labels, evaluation cells, marker subset and comparison pairs are resolved once,
from the uncorrected baseline, and are then reused unchanged by every candidate.  No
candidate may re-derive them, because a candidate that re-clusters its own output
would be scored against a partition it created.

The same freeze applies to the train/tune/validate split written to
``SamplingManifest.splits``: it is resolved once, from the manifest seed, and no
candidate sees a different one.  The split itself lives in
:mod:`batchdx.cyto.sampling.splits`, because a fit-and-apply backend needs it
independently of the evaluation draw.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from batchdx.cyto.config import Config
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.sampling.splits import SPLIT_PROPORTION_ORIGIN, build_splits
from batchdx.cyto.schemas.contracts import (
    DesignAudit,
    ProvenanceRef,
    RepresentationSettings,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import LABEL_TIER_CONFIDENCE, ComparisonSupport, LabelTier

UNGATED = "ungated"

#: The only partition scheme this build has ever had: labels, cells, markers and pairs
#: resolved once from the uncorrected baseline and reused by every candidate.  The
#: configuration key ``evaluation.partitions`` is pinned to this value, and the check
#: below is what makes that pin load-bearing rather than decorative.
FROZEN_COMMON = "frozen_common"

#: The unit the evaluation draw is balanced over, for the same reason.
SAMPLE_SAMPLING_UNIT = "sample_id"


@dataclass(frozen=True)
class GatingRule:
    """One threshold rule: ``feature`` is positive or negative about ``threshold``."""

    population: str
    feature_id: str
    sign: str  # "+" or "-"
    threshold: float
    priority: int = 0


class PartitionError(ValueError):
    """Raised when evaluation labels cannot be resolved as declared."""


def read_gating_thresholds(path: str | Path) -> list[GatingRule]:
    """Read the threshold table.

    Thresholds are fixed in the declared coordinate system *before* any candidate
    runs, and the same numbers are used before and after correction.  There is no
    per-candidate retuning, because five experienced analysts setting five different
    thresholds for the same marker is a documented failure mode.
    """

    p = Path(path)
    if not p.exists():
        raise PartitionError(f"gating threshold table not found: {p}")
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep)
    required = {"population", "feature_id", "sign", "threshold"}
    missing = required - set(df.columns)
    if missing:
        raise PartitionError(f"{p}: missing column(s) {sorted(missing)}")
    rules: list[GatingRule] = []
    for i, row in df.iterrows():
        sign = str(row["sign"]).strip()
        if sign not in {"+", "-"}:
            raise PartitionError(f"{p} row {i}: sign must be '+' or '-', got {sign!r}")
        rules.append(GatingRule(
            population=str(row["population"]).strip(),
            feature_id=str(row["feature_id"]).strip(),
            sign=sign,
            threshold=float(row["threshold"]),
            priority=int(row["priority"]) if "priority" in df.columns else i,
        ))
    return rules


def _labels_from_gating(dataset: Dataset, rules: list[GatingRule]) -> np.ndarray:
    populations: dict[str, list[GatingRule]] = {}
    order: list[str] = []
    for rule in sorted(rules, key=lambda r: (r.priority, r.population)):
        if rule.population not in populations:
            populations[rule.population] = []
            order.append(rule.population)
        populations[rule.population].append(rule)

    features = sorted({r.feature_id for r in rules})
    values = dataset.expression(features=features)
    col = {f: i for i, f in enumerate(features)}

    labels = np.full(dataset.n_cells, UNGATED, dtype=object)
    assigned = np.zeros(dataset.n_cells, dtype=bool)
    for pop in order:
        mask = np.ones(dataset.n_cells, dtype=bool)
        for rule in populations[pop]:
            v = values[:, col[rule.feature_id]]
            mask &= (v >= rule.threshold) if rule.sign == "+" else (v < rule.threshold)
        take = mask & ~assigned
        labels[take] = pop
        assigned |= take
    return labels


def _labels_from_clustering(dataset: Dataset, k: int, seed: int) -> np.ndarray:
    """Deterministic k-means on the clustering marker mask of the *uncorrected* baseline.

    These labels are the lowest-confidence tier: they are derived from the very data
    whose batch structure is under test, so they may themselves be batch-driven.
    """

    markers = dataset.marker_ids("clustering") or dataset.marker_ids()
    x = dataset.expression(features=markers)
    rng = np.random.default_rng(seed)
    n = x.shape[0]
    k = max(1, min(k, n))

    # k-means++ seeding, then Lloyd iterations; implemented here to keep the result
    # reproducible across environments rather than tied to a library's RNG.
    centres = [x[rng.integers(n)]]
    for _ in range(1, k):
        d2 = np.min(
            np.stack([np.sum((x - c) ** 2, axis=1) for c in centres], axis=0), axis=0
        )
        total = d2.sum()
        probs = d2 / total if total > 0 else np.full(n, 1.0 / n)
        centres.append(x[rng.choice(n, p=probs)])
    centres_arr = np.stack(centres)

    assign = np.zeros(n, dtype=int)
    for _ in range(50):
        dist = ((x[:, None, :] - centres_arr[None, :, :]) ** 2).sum(axis=2)
        new_assign = dist.argmin(axis=1)
        if np.array_equal(new_assign, assign):
            break
        assign = new_assign
        for j in range(k):
            members = x[assign == j]
            if len(members):
                centres_arr[j] = members.mean(axis=0)
    return np.array([f"cluster_{j}" for j in assign], dtype=object)


def resolve_labels(
    dataset: Dataset, config: Config, base_dir: str | Path = "."
) -> tuple[np.ndarray, LabelTier, str]:
    """Return frozen evaluation labels with their tier and confidence."""

    source = LabelTier(config.evaluation.label_source)
    if source is LabelTier.ANNOTATION:
        column = config.evaluation.labels
        if column not in dataset.adata.obs.columns:
            raise PartitionError(
                f"evaluation.labels column {column!r} is not present in the data; "
                "annotation labels are never invented"
            )
        labels = dataset.adata.obs[column].astype(str).to_numpy()
    elif source is LabelTier.GATING:
        rules = read_gating_thresholds(Path(base_dir) / config.evaluation.gating_thresholds)
        labels = _labels_from_gating(dataset, rules)
    else:
        labels = _labels_from_clustering(
            dataset, config.evaluation.clustering_k, config.project.seed
        )
    return labels, source, LABEL_TIER_CONFIDENCE[source]


SPLIT_FRACTIONS = {"train": 0.6, "tune": 0.2, "validation": 0.2}
SPLIT_ORDER = ("train", "tune", "validation")


def donor_splits(dataset: Dataset, seed: int) -> dict[str, list[str]]:
    """Deterministic donor-level train/tune/validation split of the biospecimens.

    Every sample of one donor lands in the same split: a donor split across splits
    would leak biological signal between tuning and validation. Donors are stratified
    by group so each split can contain both groups when the design allows it.
    Reference aliquots are technical replicates, not biospecimens, and are excluded.
    With fewer than three eligible donors no split is meaningful, so the single
    stratum is returned unsplit under the key "unsplit" and the empty keys are
    still present so consumers never key-error.
    """

    donors: dict[str, set[str]] = {}
    for sid, rec in sorted(dataset.samples.items()):
        if rec.reference_id:
            continue
        donor = rec.donor_id or sid  # without a donor factor, the sample is the unit
        donors.setdefault(donor, set()).add(rec.group or "unknown")

    splits: dict[str, list[str]] = {name: [] for name in SPLIT_ORDER}
    splits["unsplit"] = []
    if not donors:
        return splits

    rng = np.random.default_rng(seed)
    strata: dict[str, list[str]] = {}
    for donor, groups in donors.items():
        # a donor spanning groups (a data error) is stratified by its declared groups
        strata.setdefault("+".join(sorted(groups)), []).append(donor)

    total = sum(len(v) for v in strata.values())
    if total < 3:
        splits["unsplit"] = sorted(donors)
        return splits

    for _, members in sorted(strata.items()):
        ordered = sorted(members)
        rng.shuffle(ordered)
        n = len(ordered)
        boundary_train = int(round(SPLIT_FRACTIONS["train"] * n))
        boundary_tune = int(round((SPLIT_FRACTIONS["train"] + SPLIT_FRACTIONS["tune"]) * n))
        # with >= 3 donors every split gets at least one member
        boundary_train = min(max(boundary_train, 1), n - 2)
        boundary_tune = min(max(boundary_tune, boundary_train + 1), n - 1)
        splits["train"].extend(ordered[:boundary_train])
        splits["tune"].extend(ordered[boundary_train:boundary_tune])
        splits["validation"].extend(ordered[boundary_tune:])
    for name in SPLIT_ORDER:
        splits[name] = sorted(splits[name])
    return splits


def build_sampling_manifest(
    dataset: Dataset,
    config: Config,
    design: DesignAudit,
    labels: np.ndarray,
    label_tier: LabelTier,
    label_confidence: str,
) -> SamplingManifest:
    """Deterministic, sample-balanced evaluation cells plus the frozen comparison set."""

    if config.evaluation.partitions != FROZEN_COMMON:
        # Unreachable through parse_config, which pins the value.  Kept so that
        # widening the schema without widening this module fails loudly here instead
        # of quietly building the frozen partition under a name that promises another.
        raise PartitionError(
            f"evaluation.partitions={config.evaluation.partitions!r} is not implemented; "
            f"this build only builds the {FROZEN_COMMON!r} partition"
        )
    unit = config.evaluation.sampling.unit
    if unit != SAMPLE_SAMPLING_UNIT:
        # Likewise pinned by the schema.  The draw below is per sample, and the frozen
        # manifest records per-sample counts and weights, so honouring any other unit
        # would mean emitting a manifest that misdescribes the selection.
        raise PartitionError(
            f"evaluation.sampling.unit={unit!r} is not implemented; the evaluation draw is "
            f"balanced per {SAMPLE_SAMPLING_UNIT!r} and the frozen manifest records "
            f"per-sample counts, so no other unit can be recorded truthfully"
        )

    obs = dataset.obs
    obs["_label"] = labels
    rng = np.random.default_rng(config.project.seed)

    cap = config.evaluation.sampling.max_cells_per_sample
    selected: list[str] = []
    per_sample: dict[str, int] = {}
    original: dict[str, int] = {}

    for sid in sorted(obs[SAMPLE_ID].unique()):
        block = obs.loc[obs[SAMPLE_ID] == sid]
        original[sid] = int(len(block))
        ids = block[CELL_ID].to_numpy()
        if cap is not None and len(ids) > cap:
            # Sort first so the draw depends only on the seed and the cell ids, never
            # on row order in the input file.
            ids = np.sort(ids)
            ids = ids[np.sort(rng.choice(len(ids), size=cap, replace=False))]
        else:
            ids = np.sort(ids)
        per_sample[sid] = int(len(ids))
        selected.extend(ids.tolist())

    selected_set = set(selected)
    sel = obs.loc[obs[CELL_ID].isin(selected_set)]
    labels_by_cell = dict(zip(sel[CELL_ID].tolist(), sel["_label"].astype(str).tolist(),
                              strict=True))
    population_counts = sel["_label"].astype(str).value_counts().to_dict()

    warnings: list[str] = []
    min_cells = config.evaluation.sampling.min_cells_per_population
    for pop, n in sorted(population_counts.items()):
        if n < min_cells:
            warnings.append(
                f"population {pop!r} has {n} evaluation cells overall, below the declared "
                f"minimum of {min_cells}; metrics for it will report insufficient_data"
            )
    if UNGATED in population_counts:
        warnings.append(
            f"{population_counts[UNGATED]} evaluation cell(s) matched no gate and are "
            f"labelled {UNGATED!r}; they are reported, never silently dropped"
        )

    batch_of = {sid: rec.batch for sid, rec in dataset.samples.items()}
    per_batch_pop: dict[tuple[str, str], int] = {}
    for _, row in sel.iterrows():
        b = batch_of.get(str(row[SAMPLE_ID]))
        if b is not None:
            key = (b, str(row["_label"]))
            per_batch_pop[key] = per_batch_pop.get(key, 0) + 1
    for pop in sorted(population_counts):
        absent = [b for b in design.batches if per_batch_pop.get((b, pop), 0) == 0]
        if absent:
            warnings.append(
                f"population {pop!r} has no evaluation cells in batch(es) {absent}; "
                "comparisons involving those batches are unsupported for this population"
            )

    # A declared uncertainty unit is not a computed one.  Saying so here means the
    # report cannot be read as claiming donor-level uncertainty it never produced:
    # uncertainty_low/uncertainty_high/uncertainty_kind are null on every record.
    if config.evaluation.uncertainty_unit is not None:
        warnings.append(
            f"evaluation.uncertainty_unit={config.evaluation.uncertainty_unit!r} is recorded, "
            "but no unit-level uncertainty is computed in this build: the "
            "uncertainty_low/uncertainty_high/uncertainty_kind fields of every metric record "
            "are null, because donor-level uncertainty belongs to Stage 8 biological "
            "preservation, which is not implemented"
        )

    pairs = [
        (c.batch_a, c.batch_b)
        for c in design.comparisons
        if c.population is None and c.support is not ComparisonSupport.NON_IDENTIFIABLE
    ]

    # Biospecimen-level train/tune/validate.  Resolved from the same manifest seed
    # and the same frozen sample set, so a candidate cannot obtain a different
    # partition by being run in a different order.
    #
    # require_all_arms=False here, deliberately.  This manifest also serves the
    # detectors, the metrics and the report, and a design whose thinner arm holds
    # fewer biospecimens than there are splits is a reason to warn about the split,
    # not a reason to refuse to run: refusing would take the whole package down for a
    # limitation of one partition.  The deficiency is not hidden -- it becomes a
    # coverage warning rendered on the report -- and batchdx.cyto.sampling.splits
    # exposes missing_arms() so any consumer that must not use an unbalanced split
    # (notably the fit-and-apply CytoNorm backend) can refuse it.  A deficiency the
    # allocator could have fixed still raises, because that one is a defect.
    proportions = config.evaluation.sampling.split_proportions
    if proportions is None:
        proportion_origin = SPLIT_PROPORTION_ORIGIN
    else:
        proportion_origin = (
            "declared in the run configuration at evaluation.sampling.split_proportions; "
            "these user-configured shares replace the package default and are recorded "
            "here so a reader can disagree with them explicitly rather than inherit them"
        )
    assignment = build_splits(
        dataset.samples, donor_factor=config.design.donor, seed=config.project.seed,
        proportions=proportions, require_all_arms=False,
    )
    warnings.extend(assignment.warnings)

    rep = config.evaluation.representation
    representation = RepresentationSettings(
        space=rep.space, marker_mask=rep.marker_mask, scaling=rep.scaling,
        distance=rep.distance, k=rep.k,
    )

    # Equal weight per sample: a summary must not be dominated by whichever sample
    # happened to yield the most cells.
    weights = {sid: 1.0 for sid in per_sample}

    return SamplingManifest(
        dataset_id=dataset.manifest.dataset_id,
        seed=config.project.seed,
        unit=unit,
        max_cells_per_sample=cap,
        label_source=label_tier,
        label_confidence=label_confidence,
        label_column=config.evaluation.labels or "_batchdx_label",
        labels_by_cell=labels_by_cell,
        selected_cell_ids=sorted(selected_set),
        per_sample_counts=per_sample,
        per_sample_original_counts=original,
        population_counts={k: int(v) for k, v in population_counts.items()},
        population_coverage_warnings=warnings,
        comparison_pairs=pairs,
        evaluation_markers=dataset.marker_ids(rep.marker_mask),
        representation=representation,
        aggregation_weights=weights,
        splits=assignment.to_manifest(),
        split_proportions=dict(assignment.proportions),
        split_proportion_origin=proportion_origin,
        provenance=[ProvenanceRef(kind="dataset", identifier=dataset.manifest.dataset_id,
                                  content_hash=dataset.baseline_hash())],
    )
