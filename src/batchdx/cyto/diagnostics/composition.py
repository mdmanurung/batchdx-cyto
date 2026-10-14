"""Composition diagnostic, kept separate from intensity diagnostics.

A shift in how many cells of a population a batch contains is not evidence that marker
intensities need correcting.  This detector reports population frequencies per sample
and per batch, tests them at the sample level (the independent unit), and names
gating and viability differences as candidate explanations rather than asserting a
technical cause.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from batchdx.cyto.diagnostics.base import DetectorContext, eligibility_for, register_detector
from batchdx.cyto.io.dataset import SAMPLE_ID, Dataset
from batchdx.cyto.schemas.contracts import DetectorAssessment, DetectorResult
from batchdx.cyto.schemas.enums import MarkerRole, MetricStatus

MIN_SAMPLES_PER_BATCH = 2


def _bh_adjust(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, order-preserving and monotone."""

    n = len(pvalues)
    if n <= 1:
        return list(pvalues)
    order = np.argsort(pvalues)
    ranked = np.asarray(pvalues, dtype=float)[order]
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return [float(v) for v in out]


class CompositionDetector:
    """Population frequencies by sample and batch, with a sample-level test."""

    name = "composition"
    namespace = "common"
    eligibility_key = "composition"

    def run(self, dataset: Dataset, context: DetectorContext) -> DetectorResult:
        eligibility = eligibility_for(context.design, self.eligibility_key)
        if eligibility is not None and not eligibility.eligible:
            return DetectorResult(
                detector=self.name, namespace="common", run_id=context.run_id,
                status="ineligible", reason=eligibility.reason, eligibility=eligibility,
            )
        if context.sampling is None:
            return DetectorResult(
                detector=self.name, namespace="common", run_id=context.run_id,
                status="ineligible",
                reason="no frozen evaluation labels are available; composition needs populations",
                eligibility=eligibility,
            )

        labels = context.sampling.labels_by_cell
        obs = dataset.obs
        obs = obs.loc[obs["cell_id"].isin(labels.keys())].copy()
        obs["_label"] = obs["cell_id"].map(labels)

        batch_of = dataset.sample_batches()
        populations = sorted(set(labels.values()))
        samples = sorted(obs[SAMPLE_ID].unique())

        freq: dict[str, dict[str, float]] = {}
        empty_samples = []
        for sid in samples:
            block = obs.loc[obs[SAMPLE_ID] == sid]
            total = len(block)
            if total == 0:
                empty_samples.append(str(sid))
                freq[str(sid)] = {pop: float("nan") for pop in populations}
                continue
            counts = block["_label"].value_counts().to_dict()
            freq[str(sid)] = {
                pop: (counts.get(pop, 0) / total if total else float("nan"))
                for pop in populations
            }

        assessments: list[DetectorAssessment] = []
        for sid in samples:
            for pop in populations:
                value = freq[str(sid)][pop]
                if np.isfinite(value):
                    assessments.append(DetectorAssessment(
                        assessment_id=f"freq:{sid}:{pop}", detector=self.name,
                        namespace="common", scope="sample_population", sample_id=str(sid),
                        population=pop, batch=batch_of.get(str(sid)),
                        statistic="population_frequency", value=float(value),
                        status=MetricStatus.OK, evidence_class="descriptive",
                        reason="fraction of the sample's evaluation cells in this population",
                    ))
                else:
                    assessments.append(DetectorAssessment(
                        assessment_id=f"freq:{sid}:{pop}", detector=self.name,
                        namespace="common", scope="sample_population", sample_id=str(sid),
                        population=pop, batch=batch_of.get(str(sid)),
                        statistic="population_frequency", value=None,
                        status=MetricStatus.INSUFFICIENT_DATA,
                        reason="the sample has no evaluation cells, so a frequency is "
                               "undefined and is not reported as a number",
                    ))

        by_batch: dict[str, list[str]] = {}
        for sid in samples:
            b = batch_of.get(str(sid))
            if b is not None:
                by_batch.setdefault(b, []).append(str(sid))

        warnings: list[str] = []
        small = [b for b, sids in by_batch.items() if len(sids) < MIN_SAMPLES_PER_BATCH]
        tested: list[tuple[str, float, float]] = []  # (population, stat, pvalue)
        for pop in populations:
            groups = [
                [freq[sid][pop] for sid in sids]
                for b, sids in sorted(by_batch.items())
                if len(sids) >= MIN_SAMPLES_PER_BATCH
            ]
            if len(groups) < 2:
                assessments.append(DetectorAssessment(
                    assessment_id=f"freqtest:{pop}", detector=self.name, namespace="common",
                    scope="population", population=pop,
                    statistic="kruskal_batch_frequency_pvalue",
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=(
                        f"fewer than two batches have at least {MIN_SAMPLES_PER_BATCH} "
                        "samples; the sample is the independent unit, not the cell"
                    ),
                ))
                continue
            try:
                stat, pvalue = stats.kruskal(*groups)
            except ValueError as exc:
                assessments.append(DetectorAssessment(
                    assessment_id=f"freqtest:{pop}", detector=self.name, namespace="common",
                    scope="population", population=pop,
                    statistic="kruskal_batch_frequency_pvalue",
                    status=MetricStatus.UNDEFINED, reason=str(exc),
                ))
                continue
            tested.append((pop, float(stat), float(pvalue)))

        # Benjamini-Hochberg across the populations tested by this detector: each
        # population is one test, and per-population flags at raw alpha would emit
        # alpha*P false flags per run.
        adjusted = _bh_adjust([p for _, _, p in tested])
        for (pop, stat, pvalue), padj in zip(tested, adjusted, strict=True):
            assessments.append(DetectorAssessment(
                assessment_id=f"freqtest:{pop}", detector=self.name, namespace="common",
                scope="population", population=pop,
                statistic="kruskal_batch_frequency_pvalue", value=float(pvalue),
                flag=bool(padj < context.composition_alpha),
                threshold=context.composition_alpha, status=MetricStatus.OK,
                evidence_class="descriptive",
                reason=(
                    f"Kruskal-Wallis H={stat:.3g} on sample-level frequencies across "
                    f"batches; BH-adjusted p={padj:.3g} over {len(tested)} population(s) "
                    "tested. A frequency difference is not evidence that marker "
                    "intensities require correction"
                ),
            ))

        if small:
            warnings.append(
                f"batch(es) {sorted(small)} contain fewer than {MIN_SAMPLES_PER_BATCH} "
                "samples and are excluded from the frequency test"
            )

        viability = [m.feature_id for m in dataset.manifest.markers
                     if m.role in (MarkerRole.VIABILITY, MarkerRole.TECHNICAL)]
        explanations = [
            "differences in gating or in the declared threshold table",
            "differences in viability or cell recovery between batches",
            "genuine biological differences in composition between the donors in each batch",
        ]
        if viability:
            explanations.append(
                f"viability or technical channels are present ({viability}); inspect them "
                "before attributing a frequency shift to a technical batch effect"
            )

        flagged = [a.population for a in assessments
                   if a.statistic == "kruskal_batch_frequency_pvalue" and a.flag]
        if flagged:
            warnings.append(
                f"population frequencies differ across batches for {sorted(set(flagged))}. "
                "Candidate explanations: " + "; ".join(explanations)
            )

        return DetectorResult(
            detector=self.name, namespace="common", run_id=context.run_id, status="ok",
            eligibility=eligibility, assessments=assessments,
            requirements={"min_samples_per_batch": MIN_SAMPLES_PER_BATCH,
                          "independent_unit": "sample"},
            parameters={"alpha": context.composition_alpha,
                        "label_tier": str(context.sampling.label_source),
                        "candidate_explanations": explanations},
            warnings=warnings,
        )


register_detector(CompositionDetector())


def frequency_table(dataset: Dataset, labels: dict[str, str]) -> dict[str, dict[str, float]]:
    """Per-sample population frequencies; shared by the composition metric in Stage 7."""

    obs = dataset.obs
    obs = obs.loc[obs["cell_id"].isin(labels.keys())].copy()
    obs["_label"] = obs["cell_id"].map(labels)
    populations = sorted(set(labels.values()))
    out: dict[str, dict[str, float]] = {}
    for sid in sorted(obs[SAMPLE_ID].unique()):
        block = obs.loc[obs[SAMPLE_ID] == sid]
        total = len(block)
        counts = block["_label"].value_counts().to_dict()
        out[str(sid)] = {
            pop: (counts.get(pop, 0) / total if total else float(np.nan)) for pop in populations
        }
    return out
