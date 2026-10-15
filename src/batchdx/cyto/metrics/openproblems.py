"""OpenProblems metric components, extracted from the pinned Viash components.

The OpenProblems cytometry task ships each metric as a Viash component: a
``config.vsh.yaml`` declaring a container, plus a ``script.py`` / ``script.R`` and a
``helper.py`` / ``helper.R``.  B-007 recorded that extracting them without the
Nextflow runtime "was not attempted".  It is: the metric arithmetic lives in the
``helper`` file, and the three components extracted here are ordinary Python modules
with ordinary imports.  No Nextflow, no Viash and no container are involved in reading
or running them.

Six components are extracted, in the ``compatibility`` namespace, because their pinned
helper is importable here and can therefore be *run* and compared against:

``emd_mean_ct_horiz`` / ``emd_mean_ct_vert``
    ``src/metrics/emd/helper.py`` at commit ``37151479``.
``ratio_consistent_peaks``
    ``src/metrics/ratio_consistent_peaks/helper.py`` at the same commit.
``average_batch_r2_ct``
    ``src/metrics/average_batch_r2/helper.py``; needs scikit-learn, which is why this one
    was a ``common`` reimplementation until that dependency was declared.
``iLisi`` / ``cLisi``
    ``src/metrics/lisi/helper.py``; the LISI machinery is scib_metrics'.  The pinned
    calls are ``scib_metrics.nearest_neighbors.pynndescent`` and
    ``scib_metrics.lisi_knn``, both of which still exist in scib_metrics 0.6.1 with the
    same signatures, so the pinned helper runs unmodified against it.

What is verified is the arithmetic.  ``tests/scientific/test_stage7_openproblems.py``
imports the pinned helper file itself, runs it on a shared fixture, runs the code below
on the same fixture, and compares.  What is *not* verified by an upstream run is the
case-assembly loop in the components' ``script.py`` files: those read ``.h5ad``
resources (``resources_test/...``) that the pinned checkout does not carry, so the
loops here are documented transcriptions of that code rather than executed copies, and
the tests that cover them test the transcription.

Two departures from the pinned code, both deliberate:

*   ``emd/helper.py`` and ``ratio_consistent_peaks/script.py`` call ``pd.concat([])``
    when every group is too thin, which raises ``ValueError`` rather than reporting
    that the metric does not exist.  Here that case is ``insufficient_data`` with a
    reason, because a metric that could not be computed is neither a failure nor a
    zero.
*   The LISI components build a 100-nearest-neighbour graph and never say what happens
    below 101 cells.  Here that is ``insufficient_data`` with the cell count on the
    record.
*   The frozen common partition has no ``split`` and no biological ``group`` column,
    so the two grouping columns the components read are supplied differently.  Each
    substitution is named on the spec it affects.

No interval is attached to any record: the pinned components report a point statistic
and no uncertainty, and an interval invented here would not be the pinned component's.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, wasserstein_distance

from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.metrics.records import metric_record as _record
from batchdx.cyto.metrics.registry import (
    AVERAGE_BATCH_R2_CT,
    C_LISI,
    EMD_MEAN_CT_HORIZ,
    EMD_MEAN_CT_VERT,
    I_LISI,
    RATIO_CONSISTENT_PEAKS,
)
from batchdx.cyto.schemas.contracts import (
    DesignAudit,
    MetricRecord,
    MetricResult,
    SamplingManifest,
)
from batchdx.cyto.schemas.enums import MetricStatus

try:
    import scib_metrics
except ModuleNotFoundError:  # exercised by the wheel-smoke job's fresh venv
    # The LISI components are the only consumers.  The package imports and every
    # non-LISI metric works without it; the LISI rows report the obstacle instead of
    # a value, as the extras contract in pyproject.toml promises.
    scib_metrics = None


def _scib_metrics_obstacle() -> str | None:
    """Why the LISI components cannot run here, or ``None`` when they can."""

    if scib_metrics is not None:
        return None
    return (
        "scib-metrics is not installed in this environment; the LISI components are "
        "deferred until it is (install batchdx-cyto[metrics])"
    )

#: Cells required on *each* side of a comparison, as declared in the pinned components:
#: ``emd/helper.py`` skips a cell type below 50 cells, ``average_batch_r2/script.py``
#: skips below 20, and ``ratio_consistent_peaks/script.py`` skips below 100.
EMD_MIN_CELLS_PER_SIDE = 50
R2_MIN_CELLS_PER_SIDE = 20
PEAK_MIN_CELLS_PER_SIDE = 100

#: ``find_peaks`` thresholds from ``ratio_consistent_peaks/helper.py``.
PEAK_HEIGHT = 0.1
PEAK_PROMINENCE = 0.01

#: Neighbourhood size and seeds from ``lisi/helper.py``: 100 nearest neighbours for both
#: components, ``random_state=42`` for iLisi and ``random_state=0`` for cLisi.  The search
#: is approximate, so the seed is part of the metric rather than an implementation detail.
LISI_NEIGHBOURS = 100
ILISI_SEED = 42
CLISI_SEED = 0

#: The binning grid from ``emd/helper.py``: 2000 bins over [-100, 100.1), shifted by
#: 1e-7 so the left edge is excluded.  Binning with ``histogram`` rather than
#: digitising is what keeps distances physical, so moving mass one bin costs 0.1 and
#: moving it to the far end costs 200.
EMD_BIN_LOW = -100.0
EMD_BIN_HIGH = 100.1
EMD_BIN_WIDTH = 0.1
EMD_BIN_EDGE_SHIFT = 0.0000001


# --------------------------------------------------------------------------------------
# Extracted primitives.  Each is a port of the pinned helper function named in its
# docstring, with the AnnData the helper slices replaced by the array it would have
# taken out of ``layers["integrated"]``.
# --------------------------------------------------------------------------------------


def bin_array(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Port of ``emd/helper.py:bin_array``.

    Returns the bin centres, which are the explicit support of the distribution, and the
    bin probabilities, so two cell sets of different sizes can be compared.
    """

    bin_edges = np.arange(EMD_BIN_LOW, EMD_BIN_HIGH, EMD_BIN_WIDTH) + EMD_BIN_EDGE_SHIFT
    counts_per_bin, _ = np.histogram(np.asarray(values).flatten(), bins=bin_edges)
    bin_probabilities = counts_per_bin / np.sum(counts_per_bin)
    bin_indices = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_indices, bin_probabilities


def binned_emd(left_values: np.ndarray, right_values: np.ndarray) -> float:
    """Port of one marker's row of ``emd/helper.py:compute_emd``.

    Upstream returns a one-row DataFrame of every marker; for a single marker the whole
    of it is a weighted Wasserstein distance between two binned distributions.
    """

    left_values, left_weights = bin_array(left_values)
    right_values, right_weights = bin_array(right_values)
    return float(wasserstein_distance(
        left_values, right_values, left_weights, right_weights
    ))


def standardise_marker_expression(
    dist_1: np.ndarray, dist_2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Port of ``ratio_consistent_peaks/helper.py:standardise_marker_expression``.

    Both sides are standardised on the *pooled* mean and standard deviation, so that
    one pair of peak-calling thresholds means the same thing for every marker.
    """

    pooled = np.concatenate([dist_1, dist_2])
    mu, sd = pooled.mean(), pooled.std()
    return (dist_1 - mu) / sd, (dist_2 - mu) / sd


def kde_density(expression_array: np.ndarray) -> np.ndarray:
    """Port of ``ratio_consistent_peaks/helper.py:get_kde_density`` with no x grid.

    The pinned function's plotting branch is dropped: it is debugging output, it imports
    ``matplotlib`` and ``seaborn`` at call time, and nothing downstream reads a plot.  The
    bin-shift branch *is* kept, because it changes the peak count: when the density is
    highest in the first bin the grid shifts one step and a zero is prepended, so a mode
    at the edge is not silently lost.
    """

    expression_array = np.asarray(expression_array)
    min_val, max_val = expression_array.min(), expression_array.max()
    kde = gaussian_kde(np.reshape(expression_array, (1, -1)), bw_method="scott")
    x_grid = np.linspace(min_val, max_val, 100)
    density = kde(x_grid)

    if np.argmax(density) == 0 and density.size > 1:
        x_grid = np.linspace(min_val, max_val, 99)
        density = kde(x_grid)
        density = np.concatenate([[0.0], density])
        step = (max_val - min_val) / (len(x_grid)) if len(x_grid) > 1 else 0.0
        x_grid = np.concatenate([[min_val - step], x_grid])

    return density


def count_peaks(density: np.ndarray) -> int:
    """Port of ``ratio_consistent_peaks/helper.py:call_peaks``."""

    peaks, _ = find_peaks(density, prominence=PEAK_PROMINENCE, height=PEAK_HEIGHT)
    return int(len(peaks))


def peak_count(values: np.ndarray) -> int:
    """Peaks in one cell set's smoothed distribution, counted exactly as upstream does."""

    return count_peaks(kde_density(values))


# --------------------------------------------------------------------------------------
# The frozen evaluation cells, indexed the way the extracted components index them
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Panel:
    """The frozen evaluation cells of one matrix, in the columns the components read.

    The pinned components read ``adata.obs`` columns ``sample``, ``cell_type``,
    ``group`` and ``donor`` plus ``adata.layers["integrated"]``.  The frozen common
    partition has no ``split`` and its cell-type column is called ``population``, so
    ``population`` here is what the components call ``cell_type``, and ``group`` is
    supplied as the batch: upstream compares every pair of samples from the same
    *biological group*, whereas every pair of samples from the same *batch* is compared
    here.  That is a strictly harder test, because a correction which only aligns
    across batches leaves a within-batch pair untouched.
    """

    expression: np.ndarray
    markers: list[str]
    population: np.ndarray
    sample: np.ndarray
    donor: np.ndarray
    batch: np.ndarray
    #: where these cells sit in the dataset's own matrix, so the candidate matrix --
    #: which is indexed by the same cell order -- can be read the same way.
    rows: np.ndarray

    def n_cells(self) -> int:
        return int(self.expression.shape[0])


def _panel(dataset: Dataset, sampling: SamplingManifest) -> _Panel:
    """The frozen evaluation cells, with their positions in the dataset's own matrix."""

    obs = dataset.obs
    kept = obs.loc[obs[CELL_ID].isin(sampling.labels_by_cell)]
    positions = {cid: i for i, cid in enumerate(dataset.obs[CELL_ID])}
    rows = kept[CELL_ID].map(positions).to_numpy()
    features = dataset.feature_index(sampling.evaluation_markers)
    batch_of = dataset.sample_batches()
    donor_of = {sid: rec.donor_id for sid, rec in dataset.samples.items()}
    samples = kept[SAMPLE_ID].astype(str).to_numpy()

    return _Panel(
        expression=dataset.expression()[np.ix_(rows, features)],
        markers=list(sampling.evaluation_markers),
        population=kept[CELL_ID].map(sampling.labels_by_cell).to_numpy().astype(object),
        sample=samples,
        donor=np.array([str(donor_of.get(s, "")) for s in samples], dtype=object),
        batch=np.array([str(batch_of.get(s) or "") for s in samples], dtype=object),
        rows=rows,
    )


def _candidate_panel(
    panel: _Panel, dataset: Dataset, candidate_matrix: np.ndarray
) -> _Panel:
    """The same cells, read out of the candidate matrix instead of the baseline."""

    return _Panel(
        expression=candidate_matrix[np.ix_(panel.rows, dataset.feature_index(panel.markers))],
        markers=panel.markers,
        population=panel.population,
        sample=panel.sample,
        donor=panel.donor,
        batch=panel.batch,
        rows=panel.rows,
    )


def _pairs(values: list[str]) -> list[tuple[str, str]]:
    """Every unordered pair drawn from ``values`` without replacement, in sorted order."""

    unique = sorted(set(values))
    return [(a, b) for i, a in enumerate(unique) for b in unique[i + 1:]]


def _empirical_r2(values: np.ndarray, dummies: np.ndarray) -> float | None:
    """R^2 of the batch dummies for one marker, or ``None`` when it is undefined.

    The pinned ``average_batch_r2/helper.py`` fits ``sklearn``'s ``LinearRegression`` on
    an intercept plus ``pd.get_dummies(batch, drop_first=True)`` and reads
    ``model.score``, which is ``1 - SS_res / SS_tot``.  The dummies are the model's whole
    design: a donor measured in three batches contributes two of them, and dropping one
    would fit a different, wrong model.  A constant marker leaves ``SS_tot`` at zero and
    the R^2 does not exist; that is reported as undefined and never as ``0.0``.
    """

    n = values.shape[0]
    dummies = np.asarray(dummies, dtype=float)
    columns = dummies.reshape(n, -1) if dummies.ndim == 1 else dummies
    design = np.column_stack([np.ones(n), columns])
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    residuals = values - design @ beta
    ss_res = float(residuals @ residuals)
    ss_tot = float(((values - values.mean()) ** 2).sum())
    if ss_tot == 0.0:
        return None
    return 1.0 - ss_res / ss_tot


def _batch_dummies(batches: np.ndarray, levels: list[str]) -> np.ndarray:
    """``pd.get_dummies(batch, drop_first=True)`` for the batch labels in ``batches``.

    The dropped level is the first in sorted order, which is what pandas does, and it is
    the reference level the fitted intercept stands for.
    """

    return np.column_stack([
        (batches == level).astype(float) for level in levels[1:]
    ]) if len(levels) > 1 else np.zeros((batches.shape[0], 0))



def _pooled_record(
    spec, sampling: SamplingManifest, run_id: str, candidate_id: str, values: list[float],
    *, comparison_id: str | None, n_samples: int, n_donors: int,
    status_when_empty: MetricStatus, reason_when_empty: str, reason: str | None = None,
) -> MetricRecord:
    """The one row that carries the metric value itself.

    Upstream pools every marker into a single ``np.nanmean``, so the value is a mean
    over more numbers than there are rows above it.  The count is on the record, so a
    reader does not have to guess what went in.
    """

    pooled = [v for v in values if np.isfinite(v)]
    if pooled and reason is None:
        raise ValueError(f"{spec.name} pooled {len(pooled)} values but has no reason")
    return _record(
        spec, sampling, run_id=run_id, candidate_id=candidate_id,
        comparison_id=comparison_id, n_samples=n_samples, n_donors=n_donors,
        status=MetricStatus.OK if pooled else status_when_empty,
        value=float(np.mean(pooled)) if pooled else None,
        reason=reason if pooled else reason_when_empty,
        aggregation_weight=1.0,
    )


def compute_openproblems_metrics(
    dataset: Dataset,
    candidate_matrix: np.ndarray,
    sampling: SamplingManifest,
    design: DesignAudit,
    *,
    run_id: str,
    candidate_id: str,
) -> MetricResult:
    """The extracted OpenProblems components, on the frozen common partition.

    ``design`` is accepted for signature parity with :func:`compute_common_metrics`.  The
    extracted components group by donor, sample and batch and never by a declared batch
    pair, so nothing here reads it.
    """

    if not sampling.evaluation_markers:
        return MetricResult(
            run_id=run_id, candidate_id=candidate_id, metric_group="openproblems",
            status="ineligible",
            reason="no marker is flagged use_for_evaluation in the marker manifest",
        )

    panel = _panel(dataset, sampling)
    candidate = _candidate_panel(panel, dataset, candidate_matrix)
    records: list[MetricRecord] = []
    records += _horizontal_emd(panel, sampling, run_id, candidate_id)
    records += _vertical_emd(panel, sampling, run_id, candidate_id)
    records += _ratio_consistent_peaks(panel, candidate, sampling, run_id, candidate_id)
    records += _residual_batch_r2(panel, candidate, sampling, run_id, candidate_id)
    records += _ilisi(panel, candidate, sampling, run_id, candidate_id)
    records += _clisi(panel, candidate, sampling, run_id, candidate_id)

    n_failed = sum(1 for r in records if r.status is MetricStatus.FAILED)
    return MetricResult(
        run_id=run_id, candidate_id=candidate_id, metric_group="openproblems",
        records=records, status="ok" if not n_failed else "partial",
    )


# --------------------------------------------------------------------------------------
# emd_mean_ct_horiz and emd_mean_ct_vert
# --------------------------------------------------------------------------------------


def _emd_records(
    spec, sampling: SamplingManifest, run_id: str, candidate_id: str,
    comparisons: list[tuple[str, str, str, str, str]], cells: _Panel, min_cells: int,
    key_column: str,
) -> list[MetricRecord]:
    """Per-(group, side pair, population) EMD means, plus the pooled metric row.

    ``comparisons`` carries one tuple per pair: the group, the first side, the second
    side, the comparison id and a label for the group.  ``key_column`` says which panel
    column holds the group, and the two sides are the other two columns: a donor's two
    batches, or a batch's two samples.  A row is emitted even when a side is too thin,
    because a comparison that did not run has to say so.
    """

    key_values = cells.donor if key_column == "donor" else cells.batch
    side_values = cells.batch if key_column == "donor" else cells.sample
    records: list[MetricRecord] = []
    pooled: list[float] = []
    thin: list[str] = []

    for group, left, right, comparison_id, group_label in comparisons:
        left_rows = np.flatnonzero((side_values == left) & (key_values == group))
        right_rows = np.flatnonzero((side_values == right) & (key_values == group))
        for population in sorted(set(cells.population)):
            lp = left_rows[cells.population[left_rows] == population]
            rp = right_rows[cells.population[right_rows] == population]
            common = dict(
                run_id=run_id, candidate_id=candidate_id, population=population,
                comparison_id=comparison_id, n_cells=len(lp) + len(rp),
                aggregation_weight=1.0,
            )
            if len(lp) < min_cells or len(rp) < min_cells:
                thin.append(f"{group_label} {population!r} ({len(lp)}/{len(rp)} cells)")
                records.append(_record(
                    spec, sampling, status=MetricStatus.INSUFFICIENT_DATA,
                    reason=(
                        f"{group_label}, {left} versus {right}, population "
                        f"{population!r}: {len(lp)} and {len(rp)} frozen evaluation "
                        f"cells; {min_cells} are required on each side"
                    ),
                    **common,
                ))
                continue
            per_marker = [
                binned_emd(cells.expression[lp, k], cells.expression[rp, k])
                for k in range(len(cells.markers))
            ]
            pooled.extend(per_marker)
            records.append(_record(
                spec, sampling, value=float(np.mean(per_marker)), status=MetricStatus.OK,
                reason=(
                    f"mean over {len(per_marker)} markers of the binned Wasserstein "
                    f"distance between {left} and {right} cells of population "
                    f"{population!r} within {group_label}"
                ),
                **common,
            ))

    records.append(_pooled_record(
        spec, sampling, run_id, candidate_id, pooled,
        reason=(
            f"mean over all {len(pooled)} (comparison, population, marker) distances, "
            "pooled exactly as the pinned component pools them"
        ),
        comparison_id=None, n_samples=len(set(cells.sample)), n_donors=len(set(cells.donor)),
        status_when_empty=MetricStatus.INSUFFICIENT_DATA,
        reason_when_empty=(
            f"no comparison reached {min_cells} cells on both sides; thinnest: "
            + ("; ".join(thin[:4]) if thin else "none recorded")
        ),
    ))
    return records


def _batch_samples(cells: _Panel) -> dict[str, list[str]]:
    """Which samples each batch contributed frozen evaluation cells from."""

    out: dict[str, set[str]] = {}
    for sample, batch in zip(cells.sample, cells.batch, strict=True):
        if batch:
            out.setdefault(str(batch), set()).add(str(sample))
    return {batch: sorted(samples) for batch, samples in out.items()}


def _donors_batches(cells: _Panel) -> dict[str, list[str]]:
    """Which batches each donor has frozen evaluation cells in."""

    out: dict[str, set[str]] = {}
    for donor, batch in zip(cells.donor, cells.batch, strict=True):
        if batch:
            out.setdefault(str(donor), set()).add(str(batch))
    return {donor: sorted(batches) for donor, batches in out.items()}


def _donor_batch_comparisons(cells: _Panel) -> list[tuple[str, str, str, str, str]]:
    comparisons = []
    for donor, batches in sorted(_donors_batches(cells).items()):
        for batch_a, batch_b in _pairs(batches):
            comparisons.append((
                donor, batch_a, batch_b,
                f"donor:{donor}|{batch_a}__vs__{batch_b}", f"donor {donor!r}",
            ))
    return comparisons


def _sample_pair_comparisons(cells: _Panel) -> list[tuple[str, str, str, str, str]]:
    comparisons = []
    for batch, samples in sorted(_batch_samples(cells).items()):
        for sample_a, sample_b in _pairs(samples):
            comparisons.append((
                batch, sample_a, sample_b,
                f"batch:{batch}|{sample_a}__vs__{sample_b}", f"batch {batch!r}",
            ))
    return comparisons


def _no_comparisons(spec, sampling: SamplingManifest, run_id: str, candidate_id: str,
                    cells: _Panel, reason: str) -> list[MetricRecord]:
    return [_pooled_record(
        spec, sampling, run_id, candidate_id, [],
        comparison_id=None, n_samples=len(set(cells.sample)),
        n_donors=len(set(cells.donor)), status_when_empty=MetricStatus.INSUFFICIENT_DATA,
        reason_when_empty=reason,
    )]


def _horizontal_emd(
    cells: _Panel, sampling: SamplingManifest, run_id: str, candidate_id: str
) -> list[MetricRecord]:
    """``emd_mean_ct_horiz``: one donor, two batches.

    Upstream compares the two halves of a single donor, i.e. the same donor measured in
    two batches.  The frozen partition has no split, so a donor's cells are grouped by
    batch and every pair of batches that donor appears in is compared.  A donor that
    appears in only one batch -- the normal case when donors are nested inside batches,
    as they are in the paired-reference design -- has nothing to compare.
    """

    comparisons = _donor_batch_comparisons(cells)
    if not comparisons:
        return _no_comparisons(
            EMD_MEAN_CT_HORIZ, sampling, run_id, candidate_id, cells,
            "no donor in the frozen partition spans two batches, so there is no "
            "within-donor comparison to make; this component needs a design with "
            "technical replicates of one donor split across batches",
        )
    return _emd_records(
        EMD_MEAN_CT_HORIZ, sampling, run_id, candidate_id, comparisons, cells,
        EMD_MIN_CELLS_PER_SIDE, "donor",
    )


def _vertical_emd(
    cells: _Panel, sampling: SamplingManifest, run_id: str, candidate_id: str
) -> list[MetricRecord]:
    """``emd_mean_ct_vert``: every pair of samples inside one group.

    Upstream's group is a biological condition and every pair of its samples is
    compared; the partition's group column is the batch, so every pair of samples
    processed in the same batch is compared.
    """

    comparisons = _sample_pair_comparisons(cells)
    if not comparisons:
        return _no_comparisons(
            EMD_MEAN_CT_VERT, sampling, run_id, candidate_id, cells,
            "no batch in the frozen partition holds two samples, so this component has "
            "no within-group sample pair to compare",
        )
    return _emd_records(
        EMD_MEAN_CT_VERT, sampling, run_id, candidate_id, comparisons, cells,
        EMD_MIN_CELLS_PER_SIDE, "batch",
    )


# --------------------------------------------------------------------------------------
# ratio_consistent_peaks
# --------------------------------------------------------------------------------------


def _donor_population_sides(
    cells: _Panel, donor: str, population: str
) -> dict[str, np.ndarray]:
    """Cell positions of one (donor, population) triple, split by batch.

    Empty when the donor spans fewer than two batches: there is no within-donor
    comparison to make, and the caller reports that rather than inventing one.
    """

    batches = _donors_batches(cells).get(donor, [])
    if len(batches) < 2:
        return {}
    return {
        batch: np.flatnonzero(
            (cells.donor == donor) & (cells.batch == batch) & (cells.population == population)
        )
        for batch in batches
    }


def _standardised_peaks(blocks: dict[str, np.ndarray]) -> dict[str, int | None] | None:
    """Peak count per block, after standardising the blocks on their pooled mean and SD.

    This is the order upstream works in -- ``standardise_marker_expression`` pools the
    two sides and rescales both, then ``get_kde_density`` and ``call_peaks`` run on the
    rescaled values -- and the order matters: it is what makes one pair of absolute
    ``find_peaks`` thresholds mean the same thing for a marker whose spread is 0.2 and
    one whose spread is 40.  Upstream pools the two sides it is comparing; a donor here
    may have more than two, so all of its blocks are pooled together.

    ``None`` means the pooled values carry no variance at all, so there is no peak
    structure anywhere in the triple.
    """

    pooled = np.concatenate([np.asarray(block).flatten() for block in blocks.values()])
    sd = float(np.std(pooled))
    if sd == 0.0:
        return None
    mu = float(np.mean(pooled))
    return {
        name: _peaks_unless_constant((np.asarray(block).flatten() - mu) / sd)
        for name, block in blocks.items()
    }


def _ratio_consistent_peaks(
    cells: _Panel, candidate: _Panel, sampling: SamplingManifest,
    run_id: str, candidate_id: str,
) -> list[MetricRecord]:
    """``ratio_consistent_peaks``: how much peak structure the correction destroyed.

    A ``(donor, population, marker)`` triple is *preserved* when every batch of that
    donor shows the same number of peaks and that number is the one the uncorrected
    data showed -- upstream's case 1, which requires agreement with the uncorrected
    count and not merely agreement between the batches.  Upstream compares two batches;
    a donor here may span more than two, so every batch has to agree.

    Peaks are counted exactly as upstream counts them: standardise the two sides on
    their pooled mean and standard deviation, smooth with a Gaussian KDE, and call
    peaks at ``height=0.1`` and ``prominence=0.01``.  A triple whose uncorrected peaks
    disagree is upstream's case NGT and is excluded from the ratio; a triple the
    correction flattened to a constant is counted as broken rather than skipped.
    """

    records: list[MetricRecord] = []
    n_consistent = n_broken = n_excluded = 0
    thin: list[str] = []

    for donor in sorted(set(cells.donor)):
        for population in sorted(set(cells.population)):
            sides = _donor_population_sides(cells, donor, population)
            if len(sides) < 2:
                continue
            if any(len(rows) < PEAK_MIN_CELLS_PER_SIDE for rows in sides.values()):
                counts = ", ".join(f"{b}={len(r)}" for b, r in sorted(sides.items()))
                thin.append(f"donor {donor!r} population {population!r} ({counts})")
                records.append(_record(
                    RATIO_CONSISTENT_PEAKS, sampling, run_id=run_id, candidate_id=candidate_id,
                    population=population, comparison_id=f"donor:{donor}",
                    n_samples=len(sides), n_donors=1,
                    n_cells=sum(len(r) for r in sides.values()),
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=(
                        f"a peak count is not stable below {PEAK_MIN_CELLS_PER_SIDE} "
                        f"cells; this triple has {counts}"
                    ),
                ))
                continue
            for k, marker in enumerate(cells.markers):
                common = dict(
                    run_id=run_id, candidate_id=candidate_id, marker=marker,
                    population=population, comparison_id=f"donor:{donor}",
                    n_samples=len(sides), n_donors=1,
                    n_cells=sum(len(r) for r in sides.values()), aggregation_weight=1.0,
                )
                before = _standardised_peaks(
                    {b: cells.expression[rows, k] for b, rows in sides.items()}
                )
                after = _standardised_peaks(
                    {b: candidate.expression[rows, k] for b, rows in sides.items()}
                )

                if before is None:
                    n_excluded += 1
                    records.append(_record(
                        RATIO_CONSISTENT_PEAKS, sampling, status=MetricStatus.OK,
                        reason=(
                            "excluded: the uncorrected data of this triple carries no "
                            "variance at all, so there is no peak structure for the "
                            "correction to preserve or destroy"
                        ),
                        **common,
                    ))
                    continue
                if after is None or any(v is None for v in after.values()):
                    # upstream's zero-variance guard, which runs before the KDE because
                    # the pinned helper cannot peak-call a constant array at all
                    flat = (
                        sorted(sides) if after is None
                        else sorted(b for b, v in after.items() if v is None)
                    )
                    n_broken += 1
                    records.append(_record(
                        RATIO_CONSISTENT_PEAKS, sampling, status=MetricStatus.OK,
                        reason=(
                            "case 2 (flattened): the correction left no variance in "
                            f"batch(es) {flat} while the uncorrected data had "
                            f"{sorted(v for v in before.values() if v is not None)} peak(s); "
                            "a flat side is broken, not skipped"
                        ),
                        **common,
                    ))
                    continue
                if any(v is None for v in before.values()):
                    n_excluded += 1
                    records.append(_record(
                        RATIO_CONSISTENT_PEAKS, sampling, status=MetricStatus.OK,
                        reason=(
                            f"excluded: the uncorrected data is already constant in "
                            f"batch(es) "
                            f"{sorted(b for b, v in before.items() if v is None)}"
                        ),
                        **common,
                    ))
                    continue
                if len(set(before.values())) != 1:
                    n_excluded += 1
                    records.append(_record(
                        RATIO_CONSISTENT_PEAKS, sampling, status=MetricStatus.OK,
                        reason=(
                            f"case NGT (excluded): the uncorrected data already disagrees "
                            f"across the donor's batches ({sorted(before.values())} peaks)"
                        ),
                        **common,
                    ))
                    continue
                agreed = sorted(before.values())[0]
                # upstream's case 1 is not merely "the batches still agree with each
                # other" but "they agree with each other AND with the uncorrected
                # count"; a correction that turned two modes into one has changed the
                # marker even though its batches match
                if all(after[b] == agreed for b in sides):
                    n_consistent += 1
                    reason = (
                        f"case 1 (preserved): {agreed} peak(s) in every batch of this "
                        "donor, before and after the correction"
                    )
                else:
                    n_broken += 1
                    reason = (
                        f"case 2 (broken): the uncorrected data agreed at {agreed} peak(s) "
                        f"across the donor's batches and the correction leaves "
                        f"{sorted(after.values())}"
                    )
                records.append(_record(
                    RATIO_CONSISTENT_PEAKS, sampling, status=MetricStatus.OK,
                    reason=reason, **common,
                ))

    denominator = n_consistent + n_broken
    if not records:
        records.append(_record(
            RATIO_CONSISTENT_PEAKS, sampling, run_id=run_id, candidate_id=candidate_id,
            status=MetricStatus.INSUFFICIENT_DATA,
            reason=(
                "no donor in the frozen partition spans two batches, so there is no "
                "within-donor peak comparison to make"
            ),
        ))
    records.append(_record(
        RATIO_CONSISTENT_PEAKS, sampling, run_id=run_id, candidate_id=candidate_id,
        comparison_id=None, n_samples=len(set(cells.sample)), n_donors=len(set(cells.donor)),
        value=(n_consistent / denominator) if denominator else None,
        status=MetricStatus.OK if denominator else MetricStatus.UNDEFINED,
        reason=(
            f"{n_consistent} of {denominator} triples that agreed before correction still "
            f"agree; {n_excluded} triple(s) were excluded because the uncorrected data "
            f"already disagreed"
            if denominator else
            f"no triple was comparable: {n_excluded} excluded, and "
            f"{len(thin)} below {PEAK_MIN_CELLS_PER_SIDE} cells per side"
        ),
    ))
    return records


def _peaks_unless_constant(values: np.ndarray) -> int | None:
    """Peak count, or ``None`` when the values carry no variance at all.

    This is the guard ``ratio_consistent_peaks/script.py`` applies to a *corrected*
    side, which the correction flattened to a constant while the uncorrected side was
    not: such a triple is case 2 by inspection, without running the KDE.  The same test
    applied to the uncorrected side means there was never any peak structure there, so
    the triple is excluded from the ratio rather than counted on either side of it.
    """

    if float(np.std(values)) == 0.0:
        return None
    return peak_count(values)

# --------------------------------------------------------------------------------------
# average_batch_r2_ct
# --------------------------------------------------------------------------------------


def _residual_batch_r2(
    cells: _Panel, candidate: _Panel, sampling: SamplingManifest,
    run_id: str, candidate_id: str,
) -> list[MetricRecord]:
    """``average_batch_r2_ct``, extracted from the pinned helper.

    The pinned ``average_batch_r2/helper.py`` fits scikit-learn's ``LinearRegression``
    on the batch dummies and reads ``model.score``.  The arithmetic is reproduced in
    closed form by :func:`_empirical_r2` -- one dummy column per batch but the first, the
    same design pandas' ``get_dummies`` produces -- and the two are compared directly in
    ``tests/scientific/test_stage7_openproblems.py``.  What is transcribed rather than
    executed is the walk over donors and cell types in the component's ``script.py``,
    which needs ``.h5ad`` resources the pinned checkout does not carry.
    """

    records: list[MetricRecord] = []
    pooled: list[float] = []
    undefined = 0
    thin: list[str] = []

    for donor in sorted(set(cells.donor)):
        for population in sorted(set(cells.population)):
            sides = _donor_population_sides(cells, donor, population)
            if len(sides) < 2:
                continue
            if any(len(rows) < R2_MIN_CELLS_PER_SIDE for rows in sides.values()):
                counts = ", ".join(f"{b}={len(r)}" for b, r in sorted(sides.items()))
                thin.append(f"donor {donor!r} population {population!r} ({counts})")
                records.append(_record(
                    AVERAGE_BATCH_R2_CT, sampling, run_id=run_id, candidate_id=candidate_id,
                    population=population, comparison_id=f"donor:{donor}", n_donors=1,
                    n_samples=len(sides), n_cells=sum(len(r) for r in sides.values()),
                    status=MetricStatus.INSUFFICIENT_DATA,
                    reason=(
                        f"{R2_MIN_CELLS_PER_SIDE} cells are required on each side; this "
                        f"triple has {counts}"
                    ),
                ))
                continue
            batches = sorted(sides)
            rows = np.concatenate([sides[b] for b in batches])
            dummies = _batch_dummies(cells.batch[rows], batches)
            for k, marker in enumerate(cells.markers):
                common = dict(
                    run_id=run_id, candidate_id=candidate_id, marker=marker,
                    population=population, comparison_id=f"donor:{donor}", n_donors=1,
                    n_samples=len(batches), n_cells=len(rows), aggregation_weight=1.0,
                )
                r2 = _empirical_r2(candidate.expression[rows, k], dummies)
                if r2 is None:
                    undefined += 1
                    records.append(_record(
                        AVERAGE_BATCH_R2_CT, sampling, status=MetricStatus.UNDEFINED,
                        reason=(
                            f"marker {marker!r} is constant in this donor's cells, so the "
                            "batch covariate cannot explain any of a zero total variance; "
                            "the R^2 does not exist and is not reported as 0.0"
                        ),
                        **common,
                    ))
                    continue
                pooled.append(r2)
                records.append(_record(
                    AVERAGE_BATCH_R2_CT, sampling, value=r2, status=MetricStatus.OK,
                    reason=(
                        f"batch {batches[0]!r} as the reference level and {batches[1:]} as "
                        f"the contrast, over {len(rows)} cells of population {population!r} "
                        f"in donor {donor!r}"
                    ),
                    **common,
                ))

    if not records:
        records.append(_record(
            AVERAGE_BATCH_R2_CT, sampling, run_id=run_id, candidate_id=candidate_id,
            status=MetricStatus.INSUFFICIENT_DATA,
            reason=(
                "no donor in the frozen partition spans two batches, so no batch R^2 "
                "exists for this design"
            ),
        ))
    records.append(_record(
        AVERAGE_BATCH_R2_CT, sampling, run_id=run_id, candidate_id=candidate_id,
        comparison_id=None, n_samples=len(set(cells.sample)), n_donors=len(set(cells.donor)),
        value=float(np.mean(pooled)) if pooled else None,
        status=MetricStatus.OK if pooled else MetricStatus.UNDEFINED,
        reason=(
            f"mean over {len(pooled)} (donor, population, marker) fits; {undefined} "
            "triple(s) were undefined because a marker was constant"
            if pooled else
            f"no fit was defined: {len(thin)} triple(s) below the cell minimum, "
            f"{undefined} constant marker(s)"
        ),
    ))
    return records


# --------------------------------------------------------------------------------------
# iLisi and cLisi
# --------------------------------------------------------------------------------------


def confounded_groups(population: np.ndarray, batch: np.ndarray) -> list[str]:
    """Populations that appear in fewer than two batches, sorted.

    Port of the check in ``lisi/helper.py:_check_batch_group_confounding``, with the
    frozen population in the place of the component's ``group``.  A population measured
    in a single batch has no within-population batch mixing to measure, so a mixing
    score computed there would be an artefact of the design rather than a property of the
    correction.  The pinned helper skips the whole metric in that case; here the metric is
    reported ``undefined`` with the confounded populations named, which is a refusal to
    answer, not a score of zero.
    """

    counts: dict[str, set[str]] = {}
    for pop, bat in zip(population, batch, strict=True):
        counts.setdefault(str(pop), set()).add(str(bat))
    return sorted(pop for pop, batches in counts.items() if len(batches) < 2)


def local_inverse_simpson_index(
    expression: np.ndarray, labels: np.ndarray, *, seed: int
) -> np.ndarray:
    """Per-cell LISI over a 100-nearest-neighbour graph.  Port of the pinned call pair.

    ``lisi/helper.py`` builds the graph with
    ``scib_metrics.nearest_neighbors.pynndescent(layer, n_neighbors=100,
    random_state=seed)`` and scores it with ``scib_metrics.lisi_knn(graph, labels)``.  The
    seed is part of the metric, not an implementation detail: the search is approximate,
    so without it the same data would give a different neighbourhood structure each run.
    """

    graph = scib_metrics.nearest_neighbors.pynndescent(
        np.asarray(expression, dtype=np.float64), n_neighbors=LISI_NEIGHBOURS,
        random_state=seed,
    )
    return scib_metrics.lisi_knn(graph, np.asarray(labels))


def _lisi_record(
    spec, sampling: SamplingManifest, run_id: str, candidate_id: str, value: float | None,
    *, per_cell: np.ndarray, n_cells: int, n_levels: int, n_donors: int,
    status: MetricStatus, reason: str,
) -> MetricRecord:
    """One row for a whole-candidate mixing score, with the per-cell summary on it.

    There is one score per candidate and one row says so.  A record per cell would be
    6000 rows of the same number, and the reader wants the aggregate; the per-cell
    distribution is described on the record instead of enumerated.
    """

    finite = int(np.isfinite(per_cell).sum())
    return _record(
        spec, sampling, run_id=run_id, candidate_id=candidate_id,
        comparison_id=None, n_cells=n_cells, n_donors=n_donors, n_samples=n_donors,
        value=value, status=status, reason=reason, aggregation_weight=1.0,
        coverage=float(finite) / n_cells if n_cells else None,
    )


def _ilisi(
    cells: _Panel, candidate: _Panel, sampling: SamplingManifest,
    run_id: str, candidate_id: str,
) -> list[MetricRecord]:
    """``iLisi``: how well the batches mixed, over every frozen evaluation cell.

    Upstream computes it globally rather than per group, and the component's script
    averages the two splits; the frozen partition holds one candidate, so the value below
    is the single-split computation.
    """

    n_cells = candidate.n_cells()
    batches = sorted({str(b) for b in candidate.batch if b})
    if n_cells <= LISI_NEIGHBOURS:
        return [_lisi_record(
            I_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(batches), n_donors=len(set(cells.donor)),
            status=MetricStatus.INSUFFICIENT_DATA,
            reason=(
                f"a {LISI_NEIGHBOURS}-nearest-neighbour graph needs more than "
                f"{LISI_NEIGHBOURS} cells; the frozen partition has {n_cells}"
            ),
        )]
    if len(batches) < 2:
        return [_lisi_record(
            I_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(batches), n_donors=len(set(cells.donor)),
            status=MetricStatus.UNDEFINED,
            reason=(
                f"the frozen partition has one batch ({batches}), so there is nothing to "
                "mix; the normalisation 1/(n_batches - 1) has no value"
            ),
        )]
    confounded = confounded_groups(cells.population, cells.batch)
    if confounded:
        return [_lisi_record(
            I_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(batches), n_donors=len(set(cells.donor)),
            status=MetricStatus.UNDEFINED,
            reason=(
                f"batch is confounded with population in {confounded}, which appear in a "
                "single batch each; a mixing score over them would measure the design "
                "rather than the correction, so none is reported"
            ),
        )]

    obstacle = _scib_metrics_obstacle()
    if obstacle:
        return [_lisi_record(
            I_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(batches), n_donors=len(set(cells.donor)),
            status=MetricStatus.FAILED, reason=obstacle,
        )]
    per_cell = local_inverse_simpson_index(
        candidate.expression, cells.batch, seed=ILISI_SEED
    )
    raw = float(np.nanmedian(per_cell))
    value = (raw - 1.0) / (len(batches) - 1)
    return [_lisi_record(
        I_LISI, sampling, run_id, candidate_id, value, per_cell=per_cell,
        n_cells=n_cells, n_levels=len(batches), n_donors=len(set(cells.donor)),
        status=MetricStatus.OK,
        reason=(
            f"median per-cell LISI {raw:.4f} over {int(np.isfinite(per_cell).sum())} of "
            f"{n_cells} cells, normalised by (1, {len(batches) - 1}) with the pinned "
            f"seed {ILISI_SEED}; the frozen population is the component's group and the "
            "whole partition is its single split"
        ),
    )]


def _clisi(
    cells: _Panel, candidate: _Panel, sampling: SamplingManifest,
    run_id: str, candidate_id: str,
) -> list[MetricRecord]:
    """``cLisi``: whether the populations stayed apart, over the same graph."""

    n_cells = candidate.n_cells()
    populations = sorted({str(p) for p in cells.population if p})
    n_donors = len(set(cells.donor))
    if n_cells <= LISI_NEIGHBOURS:
        return [_lisi_record(
            C_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(populations), n_donors=n_donors,
            status=MetricStatus.INSUFFICIENT_DATA,
            reason=(
                f"a {LISI_NEIGHBOURS}-nearest-neighbour graph needs more than "
                f"{LISI_NEIGHBOURS} cells; the frozen partition has {n_cells}"
            ),
        )]
    if len(populations) < 2:
        return [_lisi_record(
            C_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(populations), n_donors=n_donors,
            status=MetricStatus.UNDEFINED,
            reason=(
                f"the frozen partition has {len(populations)} population, so the "
                "separation normalisation (n - 1) has no value"
            ),
        )]

    obstacle = _scib_metrics_obstacle()
    if obstacle:
        return [_lisi_record(
            C_LISI, sampling, run_id, candidate_id, None, per_cell=np.array([]),
            n_cells=n_cells, n_levels=len(populations), n_donors=n_donors,
            status=MetricStatus.FAILED, reason=obstacle,
        )]
    per_cell = local_inverse_simpson_index(
        candidate.expression, cells.population, seed=CLISI_SEED
    )
    raw = float(np.nanmedian(per_cell))
    value = (len(populations) - raw) / (len(populations) - 1)
    return [_lisi_record(
        C_LISI, sampling, run_id, candidate_id, value, per_cell=per_cell,
        n_cells=n_cells, n_levels=len(populations), n_donors=n_donors,
        status=MetricStatus.OK,
        reason=(
            f"median per-cell LISI {raw:.4f} over {int(np.isfinite(per_cell).sum())} of "
            f"{n_cells} cells against {len(populations)} populations, normalised by "
            f"({len(populations)}, {len(populations) - 1}) with the pinned seed "
            f"{CLISI_SEED}; 1 means every neighbourhood is pure and 0 that they are fully "
            "mixed"
        ),
    )]
