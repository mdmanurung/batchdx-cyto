"""Deterministic synthetic cytometry datasets used by examples and by the test suite.

The generator is a *fixture*, not a scientific claim.  It produces a mixture of
lineage populations with declared marker means, then applies separately controllable
technical and biological perturbations so that a test can state exactly which effect
it expects a method or metric to see.

The seed defaults here are :data:`batchdx.cyto.config.DEFAULT_SEED`, the same value
:attr:`batchdx.cyto.config.ProjectConfig.seed` carries.  A run's seed therefore has one
source, and a fixture built with an explicit seed is the same seed the run records.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from batchdx.cyto.config import DEFAULT_SEED

DEFAULT_POPULATIONS: dict[str, dict[str, float]] = {
    # population -> marker -> mean on the declared (asinh-like) scale
    "T_CD4": {"CD3": 3.4, "CD4": 3.2, "CD8": 0.3, "CD19": 0.2, "CD14": 0.2, "CD56": 0.2,
              "HLADR": 0.6, "CD69": 0.5},
    "T_CD8": {"CD3": 3.3, "CD4": 0.3, "CD8": 3.1, "CD19": 0.2, "CD14": 0.2, "CD56": 0.3,
              "HLADR": 0.6, "CD69": 0.5},
    "B": {"CD3": 0.2, "CD4": 0.2, "CD8": 0.2, "CD19": 3.3, "CD14": 0.2, "CD56": 0.2,
          "HLADR": 2.9, "CD69": 0.4},
    "Mono": {"CD3": 0.2, "CD4": 1.1, "CD8": 0.2, "CD19": 0.2, "CD14": 3.5, "CD56": 0.2,
             "HLADR": 2.6, "CD69": 0.4},
    "NK": {"CD3": 0.3, "CD4": 0.2, "CD8": 0.9, "CD19": 0.2, "CD14": 0.2, "CD56": 3.0,
           "HLADR": 0.5, "CD69": 0.8},
}

DEFAULT_FREQUENCIES = {"T_CD4": 0.34, "T_CD8": 0.20, "B": 0.14, "Mono": 0.22, "NK": 0.10}

LINEAGE_MARKERS = ["CD3", "CD4", "CD8", "CD19", "CD14", "CD56"]
FUNCTIONAL_MARKERS = ["HLADR", "CD69"]


@dataclass
class SampleSpec:
    """One simulated sample."""

    sample_id: str
    donor_id: str
    batch: str
    group: str
    reference_id: str | None = None
    reference_role: str = "none"
    n_cells: int = 2000
    stimulation: str = "none"
    timepoint: str = "T0"
    frequencies: dict[str, float] | None = None


@dataclass
class SimulationSpec:
    """Full description of a synthetic dataset."""

    samples: list[SampleSpec]
    seed: int = DEFAULT_SEED
    sd: float = 0.35
    #: Strength of a shared per-cell latent factor within each population.  Real
    #: cytometry markers co-vary within a cell type (cell size, staining index,
    #: activation state); without that the joint-structure metric would have nothing
    #: to measure and a joint-destroying corruption would be undetectable by
    #: construction rather than by the metric working.
    within_population_correlation: float = 0.6
    #: batch -> marker -> additive shift on the declared scale (a technical offset)
    batch_shift: dict[str, dict[str, float]] = field(default_factory=dict)
    #: batch -> marker -> multiplicative spread factor (a technical dispersion change)
    batch_spread: dict[str, dict[str, float]] = field(default_factory=dict)
    #: group -> marker -> additive shift (a genuine biological difference)
    group_shift: dict[str, dict[str, float]] = field(default_factory=dict)
    #: sample -> marker -> shift accumulated linearly across acquisition order
    drift: dict[str, dict[str, float]] = field(default_factory=dict)
    populations: dict[str, dict[str, float]] = field(
        default_factory=lambda: {k: dict(v) for k, v in DEFAULT_POPULATIONS.items()}
    )
    frequencies: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_FREQUENCIES))
    scale: str = "asinh_cofactor5"


def simulate(spec: SimulationSpec) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict]:
    """Generate an AnnData plus the sample, marker and processing manifests.

    The real expression matrix is written to ``adata.X`` *and* to the declared
    ``preprocessed`` layer.  An earlier build wrote zeros to ``X`` on the grounds
    that the declared layer is the only source of truth; that zero was a trap for
    any generic AnnData reader of the fixture, so it is gone.  Within this package
    the declared layer remains the source of truth
    (:meth:`batchdx.cyto.io.dataset.Dataset.expression` reads it and nothing reads
    ``X``).
    """

    markers = LINEAGE_MARKERS + FUNCTIONAL_MARKERS
    rng = np.random.default_rng(spec.seed)

    blocks, obs_rows = [], []
    for sample in spec.samples:
        freqs = sample.frequencies or spec.frequencies
        pops = list(freqs)
        probs = np.array([freqs[p] for p in pops], dtype=float)
        probs = probs / probs.sum()
        counts = rng.multinomial(sample.n_cells, probs)

        values = np.zeros((sample.n_cells, len(markers)), dtype=np.float64)
        labels = np.empty(sample.n_cells, dtype=object)

        start = 0
        for pop, n in zip(pops, counts, strict=True):
            if n == 0:
                continue
            means = np.array([spec.populations[pop][m] for m in markers], dtype=float)
            rho = float(spec.within_population_correlation)
            shared = rng.normal(size=(n, 1))
            independent = rng.normal(size=(n, len(markers)))
            noise = rho * shared + np.sqrt(max(0.0, 1.0 - rho * rho)) * independent
            block = means + spec.sd * noise
            values[start:start + n, :] = block
            labels[start:start + n] = pop
            start += n

        order = rng.permutation(sample.n_cells)
        values, labels = values[order], labels[order]

        for m_i, m in enumerate(markers):
            values[:, m_i] += spec.batch_shift.get(sample.batch, {}).get(m, 0.0)
            values[:, m_i] += spec.group_shift.get(sample.group, {}).get(m, 0.0)
            spread = spec.batch_spread.get(sample.batch, {}).get(m)
            if spread is not None:
                centre = values[:, m_i].mean()
                values[:, m_i] = centre + (values[:, m_i] - centre) * spread
            slope = spec.drift.get(sample.sample_id, {}).get(m)
            if slope is not None:
                values[:, m_i] += slope * np.linspace(0.0, 1.0, sample.n_cells)

        blocks.append(values)
        obs_rows.append(pd.DataFrame({
            "sample_id": np.repeat(sample.sample_id, sample.n_cells),
            "original_event_index": np.arange(sample.n_cells, dtype=int),
            "acquisition_order": np.arange(sample.n_cells, dtype=int),
            "true_population": labels,
        }))

    X = np.vstack(blocks)
    obs = pd.concat(obs_rows, ignore_index=True)
    # X carries the same values as the declared "preprocessed" layer.  The declared
    # layer stays the source of truth for this package, but a zero X in a fixture is
    # not a convention, it is a footgun for every external reader.
    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(markers, dtype=object)),
    )
    adata.layers["preprocessed"] = X

    sample_table = pd.DataFrame([{
        "sample_id": s.sample_id,
        "donor_id": s.donor_id,
        "biospecimen_id": f"{s.donor_id}_{s.timepoint}",
        "batch": s.batch,
        "group": s.group,
        "stimulation": s.stimulation,
        "timepoint": s.timepoint,
        "reference_id": s.reference_id or "",
        "reference_role": s.reference_role,
    } for s in spec.samples])

    marker_table = pd.DataFrame([{
        "feature_id": m,
        "channel": f"Ch_{m}",
        "marker": m,
        "panel_id": "panelA",
        "role": "lineage" if m in LINEAGE_MARKERS else "functional",
        "use_for_clustering": m in LINEAGE_MARKERS,
        "use_for_correction": True,
        "use_for_evaluation": True,
    } for m in markers])

    processing = {
        "compensated": True,
        "unmixed": False,
        "compensation_source": "simulated; no instrument spillover applied",
        "qc_history": ["simulated data: no QC performed"],
        "transformation": "asinh",
        "transformation_parameters": {"cofactor": 5.0},
        "scale": spec.scale,
        "units": "arcsinh units",
        "corrected": False,
        "notes": "synthetic fixture produced by batchdx.cyto.simulate",
    }
    return adata, sample_table, marker_table, processing


def write_simulation(spec: SimulationSpec, out_dir: str | Path) -> dict[str, Path]:
    """Write a simulated dataset and its manifests into ``out_dir``."""

    out = Path(out_dir)
    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / "config").mkdir(parents=True, exist_ok=True)

    adata, sample_table, marker_table, processing = simulate(spec)
    paths = {
        "h5ad": out / "data" / "preprocessed.h5ad",
        "samples": out / "config" / "samples.tsv",
        "markers": out / "config" / "markers.tsv",
        "processing": out / "config" / "processing.json",
    }
    adata.write_h5ad(paths["h5ad"])
    sample_table.to_csv(paths["samples"], sep="\t", index=False)
    marker_table.to_csv(paths["markers"], sep="\t", index=False)
    paths["processing"].write_text(json.dumps(processing, indent=2) + "\n", encoding="utf-8")
    return paths


# --------------------------------------------------------------------------------------
# Named scenarios used by the examples and the scientific suite
# --------------------------------------------------------------------------------------


def paired_reference_spec(seed: int = DEFAULT_SEED, n_cells: int = 1500) -> SimulationSpec:
    """Three batches, a measured reference aliquot in each, balanced biology.

    Batch B and C carry a technical offset on CD4 and HLADR; the biological
    difference between groups sits on CD69 and is present in every batch.
    """

    samples: list[SampleSpec] = []
    for b_i, batch in enumerate(["B1", "B2", "B3"]):
        for ref in (("ref", "REFDONOR", "REF1"), ("ref2", "REFDONOR2", "REF2")):
            sid, donor, ref_id = ref
            samples.append(SampleSpec(
                sample_id=f"{sid}_{batch}", donor_id=donor, batch=batch, group="reference",
                reference_id=ref_id, reference_role="technical_replicate", n_cells=n_cells,
            ))
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                donor = f"D{b_i}{k}{tag}"
                samples.append(SampleSpec(
                    sample_id=f"{donor}_{batch}", donor_id=donor, batch=batch,
                    group=group, n_cells=n_cells,
                ))
    return SimulationSpec(
        samples=samples,
        seed=seed,
        batch_shift={"B2": {"CD4": 0.8, "HLADR": -0.5}, "B3": {"CD4": -0.6, "HLADR": 0.7}},
        group_shift={"case": {"CD69": 0.9}},
    )


def subtle_effect_spec(seed: int = 473, n_cells: int = 1500) -> SimulationSpec:
    """Batch shifts of 0.1-0.2 SD: 0.05-0.1 on the declared scale (sd = 0.35).

    The layout matches :func:`paired_reference_spec`, including the 0.9 biological
    case/control difference on CD69, but the technical offsets are an order of
    magnitude smaller.  A correction that removes these small offsets must leave the
    biological contrast in place, so a test can state a tolerance that is small
    relative to the biology yet large relative to the residual technical error.
    """

    samples: list[SampleSpec] = []
    for b_i, batch in enumerate(["B1", "B2", "B3"]):
        samples.append(SampleSpec(
            sample_id=f"ref_{batch}", donor_id="REFDONOR", batch=batch, group="reference",
            reference_id="REF1", reference_role="technical_replicate", n_cells=n_cells,
        ))
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                donor = f"D{b_i}{k}{tag}"
                samples.append(SampleSpec(
                    sample_id=f"{donor}_{batch}", donor_id=donor, batch=batch,
                    group=group, n_cells=n_cells,
                ))
    return SimulationSpec(
        samples=samples,
        seed=seed,
        batch_shift={
            "B2": {"CD4": 0.08, "HLADR": -0.06},
            "B3": {"CD4": -0.07, "HLADR": 0.09},
        },
        group_shift={"case": {"CD69": 0.9}},
    )


def rare_population_spec(seed: int = 473, n_cells: int = 1500) -> SimulationSpec:
    """A fifth population below the 30-cell reporting floor in some samples.

    Every sample devotes either 0.05% or 0.2% of its cells to the rare population
    (roughly one to three cells per 1500-cell sample), so every (population, batch)
    group and every group-level contrast sits far below the 30-cell minimum.  Any
    metric that needs a minimum number of cells must report ``insufficient_data``
    for this population, never a zero.
    """

    populations = {k: dict(v) for k, v in DEFAULT_POPULATIONS.items()}
    populations["Rare"] = {
        "CD3": 0.4, "CD4": 0.3, "CD8": 0.3, "CD19": 1.8, "CD14": 0.2, "CD56": 1.2,
        "HLADR": 2.4, "CD69": 0.9,
    }
    samples = []
    for b_i, batch in enumerate(["B1", "B2", "B3"]):
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                donor = f"D{b_i}{k}{tag}"
                freqs = dict(DEFAULT_FREQUENCIES)
                freqs["Rare"] = 0.0005 if (b_i + k) % 2 == 0 else 0.002
                samples.append(SampleSpec(
                    sample_id=f"{donor}_{batch}", donor_id=donor, batch=batch,
                    group=group, n_cells=n_cells, frequencies=freqs,
                ))
    return SimulationSpec(
        samples=samples,
        seed=seed,
        populations=populations,
        batch_shift={"B2": {"CD4": 0.8}, "B3": {"CD4": -0.6}},
        group_shift={"case": {"CD69": 0.9}},
    )


def confounded_spec(seed: int = 473, n_cells: int = 1500) -> SimulationSpec:
    """Two batches, perfectly confounded with biological group, no reference."""

    samples = []
    for k in range(3):
        samples.append(SampleSpec(sample_id=f"C{k}_B1", donor_id=f"C{k}", batch="B1",
                                  group="case", n_cells=n_cells))
        samples.append(SampleSpec(sample_id=f"H{k}_B2", donor_id=f"H{k}", batch="B2",
                                  group="control", n_cells=n_cells))
    return SimulationSpec(
        samples=samples,
        seed=seed,
        batch_shift={"B2": {"CD4": 0.7}},
        group_shift={"case": {"CD69": 0.9}},
    )


def overlapping_no_reference_spec(seed: int = DEFAULT_SEED, n_cells: int = 1500) -> SimulationSpec:
    """Three batches, biologically overlapping, but no measured technical replicate."""

    samples = []
    for b_i, batch in enumerate(["B1", "B2", "B3"]):
        for k in range(2):
            for group, tag in (("case", "CA"), ("control", "CO")):
                donor = f"D{b_i}{k}{tag}"
                samples.append(SampleSpec(sample_id=f"{donor}_{batch}", donor_id=donor,
                                          batch=batch, group=group, n_cells=n_cells))
    return SimulationSpec(
        samples=samples,
        seed=seed,
        batch_shift={"B2": {"CD4": 0.8}, "B3": {"CD4": -0.6}},
        group_shift={"case": {"CD69": 0.9}},
    )
