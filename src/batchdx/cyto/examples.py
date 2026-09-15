"""Shipped examples that initialise every file they reference.

``batchdx-cyto init`` writes a complete, runnable project into a directory the user
chooses: synthetic data, the three manifests, the gating threshold table and a
configuration whose paths are all relative to that directory.  No example depends on a
path that exists only on a developer's machine.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from batchdx.cyto.simulate import (
    SimulationSpec,
    confounded_spec,
    overlapping_no_reference_spec,
    paired_reference_spec,
    write_simulation,
)

#: Thresholds in the declared coordinate system of the synthetic fixture
#: (asinh, cofactor 5).  Populations are matched in ``priority`` order and the first
#: match wins, so more specific gates come first.
GATING_THRESHOLDS = [
    ("Mono", "CD14", "+", 1.75, 0),
    ("B", "CD19", "+", 1.75, 1),
    ("NK", "CD3", "-", 1.75, 2),
    ("NK", "CD56", "+", 1.65, 2),
    ("T_CD4", "CD3", "+", 1.75, 3),
    ("T_CD4", "CD4", "+", 1.75, 3),
    ("T_CD8", "CD3", "+", 1.75, 4),
    ("T_CD8", "CD8", "+", 1.65, 4),
]

EXAMPLES: dict[str, dict] = {
    "paired_reference": {
        "spec": paired_reference_spec,
        "description": (
            "three batches, a measured reference aliquot in each, balanced case/control "
            "biology, and a technical offset on CD4 and HLA-DR"
        ),
        "design": {
            "reference_identity": "reference_id",
            "reference_role": "reference_role",
            "require_shared_reference": True,
            "minimum_batches": 3,
        },
        "candidates": ["uncorrected", "cycombine_default", "cytonorm_reference"],
    },
    "no_reference": {
        "spec": overlapping_no_reference_spec,
        "description": (
            "three batches with biologically overlapping groups but no measured technical "
            "replicate; evidence is conditional, not reference-supported"
        ),
        "design": {
            "reference_identity": None,
            "reference_role": None,
            "require_shared_reference": False,
            "minimum_batches": 3,
        },
        "candidates": ["uncorrected", "cycombine_default"],
    },
    "confounded": {
        "spec": confounded_spec,
        "description": (
            "two batches perfectly confounded with biological group and no reference; a "
            "diagnostic-only design that must not yield a correction recommendation"
        ),
        "design": {
            "reference_identity": None,
            "reference_role": None,
            "require_shared_reference": True,
            "minimum_batches": 3,
        },
        "candidates": ["uncorrected", "cycombine_default"],
    },
}

CANDIDATE_SPECS: dict[str, dict] = {
    "uncorrected": {"id": "uncorrected", "backend": "none", "parameters": {}},
    "cycombine_default": {"id": "cycombine_default", "backend": "cycombine", "parameters": {}},
    "cytonorm_reference": {"id": "cytonorm_reference", "backend": "cytonorm", "parameters": {}},
}


def _config_dict(name: str, example: dict) -> dict:
    return {
        "schema_version": 2,
        "project": {"id": f"{name}_example", "seed": 473},
        "input": {
            "format": "h5ad",
            "path": "data/preprocessed.h5ad",
            "expression_layer": "preprocessed",
            "samples": "config/samples.tsv",
            "markers": "config/markers.tsv",
            "processing_manifest": "config/processing.json",
        },
        "preprocessing": {"enabled": False},
        "design": {
            "correction_batch": "batch",
            "protect": ["group"],
            "donor": "donor_id",
            "reference_identity": example["design"]["reference_identity"],
            "reference_role": example["design"]["reference_role"],
            "allow_exploratory_override": False,
            "minimum_batches": example["design"]["minimum_batches"],
            "require_shared_reference": example["design"]["require_shared_reference"],
            "minimum_batches_origin": (
                "shipped example default: below three batches a batch effect cannot be "
                "separated from a single outlying batch, and the IQR and cluster detectors "
                "are degenerate"
            ),
        },
        "evaluation": {
            "label_source": "gating",
            "gating_thresholds": "config/thresholds.tsv",
            "representation": {
                "space": "marker", "marker_mask": "evaluation",
                "scaling": "none", "distance": "euclidean", "k": 30,
            },
            "sampling": {
                "unit": "sample_id", "max_cells_per_sample": 1200,
                "min_cells_per_population": 50,
            },
            "partitions": "frozen_common",
            "uncertainty_unit": "donor_id",
            "policy": "conservative",
        },
        "correction": {
            "control_policy": "measured_only",
            "candidates": [CANDIDATE_SPECS[c] for c in example["candidates"]],
        },
        "detectors": {
            "acquisition_drift": True,
            "composition": True,
            "cycombine_express": False,
            "cycombine_cluster": False,
            "cytobatchflagr": False,
        },
        "upstream": {"lock": "config/upstream_lock.yaml", "checkout_dir": None,
                     "rscript": "Rscript"},
        "report": {"offline": True, "redact_identifiers": True,
                   "title": f"batchdx-cyto — {name} example"},
        "exports": {"h5ad": True, "fcs": False},
    }


def materialise_example(
    name: str, target: str | Path, cells_per_sample: int = 1500, seed: int = 473
) -> dict[str, Path]:
    """Write a complete example project and return the paths created."""

    if name not in EXAMPLES:
        raise KeyError(f"unknown example {name!r}; available: {sorted(EXAMPLES)}")
    example = EXAMPLES[name]
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)

    spec: SimulationSpec = example["spec"](seed=seed, n_cells=cells_per_sample)
    paths = write_simulation(spec, target)

    thresholds = pd.DataFrame(
        GATING_THRESHOLDS, columns=["population", "feature_id", "sign", "threshold", "priority"]
    )
    thresholds_path = target / "config" / "thresholds.tsv"
    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    paths["thresholds"] = thresholds_path

    config_path = target / "config" / "config.yaml"
    config_path.write_text(
        "# " + example["description"] + "\n"
        + yaml.safe_dump(_config_dict(name, example), sort_keys=False),
        encoding="utf-8",
    )
    paths["config"] = config_path

    readme = target / "README.md"
    readme.write_text(
        f"# batchdx-cyto example: {name}\n\n{example['description']}\n\n"
        "```\n"
        "batchdx-cyto validate --config config/config.yaml\n"
        "batchdx-cyto diagnose --config config/config.yaml --outdir results/run1\n"
        "batchdx-cyto benchmark --config config/config.yaml --outdir results/run1\n"
        "```\n\n"
        "Equivalent native Snakemake invocation:\n\n"
        "```\n"
        "snakemake --snakefile $(batchdx-cyto info --json | python -c "
        "'import json,sys; print(json.load(sys.stdin)[\"workflow\"])') \\\n"
        "    --cores 4 --config config=config/config.yaml outdir=results/run1 mode=benchmark\n"
        "```\n\n"
        "The data here are synthetic fixtures produced by `batchdx.cyto.simulate`. They are\n"
        "useful for exercising the workflow and are not a scientific claim about real\n"
        "cytometry.\n",
        encoding="utf-8",
    )
    paths["readme"] = readme

    manifest = target / "config" / "example_manifest.json"
    manifest.write_text(json.dumps({
        "example": name,
        "description": example["description"],
        "seed": seed,
        "cells_per_sample": cells_per_sample,
        "n_samples": len(spec.samples),
        "generated_by": "batchdx.cyto.examples.materialise_example",
    }, indent=2) + "\n", encoding="utf-8")
    paths["example_manifest"] = manifest
    return paths
