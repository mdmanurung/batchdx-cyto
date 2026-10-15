"""cyCombine correction backend.

The adapter targets the pinned upstream API recorded in ``config/upstream_lock.yaml``
and audited by ``scripts/audit_upstream.R``: ``normalize`` then ``create_som`` then
``correct_data``.

Two upstream properties drive the design:

* ``correct_data`` takes ``method = c("ComBat", "ComBat_seq")`` and the package imports
  both from *sva*.  ``ComBat_seq`` is a negative-binomial model for counts and is
  wrong for transformed cytometry, so this adapter pins ``method = "ComBat"`` and
  refuses any configuration that would reach the count model.
* ``correct_data`` runs ComBat *within each SOM node*.  Identifiability therefore has
  to hold inside the correction partitions and not only globally, so the R side
  reports per-node batch and covariate coverage back to Python.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from batchdx.cyto.adapters.base import (
    BackendCapabilities,
    BackendContext,
    BackendError,
    BackendNotAvailable,
    make_spec,
    register_backend,
    resolve_parameters,
)
from batchdx.cyto.adapters.r_runner import (
    RNotAvailable,
    RScriptFailed,
    r_package_version,
    r_version,
    run_r_script,
)
from batchdx.cyto.adapters.validate import changed_samples, correction_magnitude, validate_output
from batchdx.cyto.io.binary import read_matrix, write_matrix
from batchdx.cyto.io.dataset import CELL_ID, SAMPLE_ID, Dataset
from batchdx.cyto.io.hashing import matrix_hash
from batchdx.cyto.schemas.contracts import CorrectionResult, ProvenanceRef, SoftwareVersion
from batchdx.cyto.schemas.enums import CorrectionStatus, ReferenceRole

R_SCRIPT = "cycombine_correct.R"


class CyCombineBackend:
    """Wrap ``cyCombine::correct_data`` on the continuous-expression route."""

    capabilities = BackendCapabilities(
        name="cycombine",
        required_sample_fields=("batch",),
        reference_layouts=("none", "shared_reference_optional", "anchor"),
        supports_covariates=True,
        supports_marker_mask=True,
        input_scale="transformed (asinh or logicle); not counts",
        output_scale="same transformed scale as the input",
        supports_fit_apply=False,
        fit_apply_verified=False,
        resource_profile=(
            "memory: all cells of all batches are held jointly; ComBat is fitted per SOM "
            "node across the whole dataset, so per-file chunking is not equivalent"
        ),
        accepted_parameters={
            "norm_method": "scale",
            "cluster_method": "kohonen",
            "distf": "euclidean",
            "xdim": 8,
            "ydim": 8,
            "rlen": 10,
            "nClus": None,
            "parametric": True,
            "covar": None,
            "anchor": None,
            "ref_batch": None,
        },
        rejected_parameters={
            "method": (
                "the correction route is pinned to ComBat. cyCombine also exposes "
                "ComBat_seq, a negative-binomial model for count data, which is not a "
                "valid model for transformed cytometry expression."
            ),
            "mc.cores": "parallelism is owned by Snakemake resources, not by adapter parameters",
        },
        notes=(
            "fit and apply are not separable in the audited API: correct_data returns a "
            "corrected table, not a reusable model, so frozen application to unseen "
            "batches is not advertised. The anchor parameter is explicitly experimental "
            "upstream (correct_data docstring: 'Experimental ... may be confounded with "
            "the condition'); this adapter builds it from reference identity and never "
            "presents an anchor-corrected run as ordinary ComBat output."
        ),
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [
            r_version(context.rscript),
            r_package_version("cyCombine", context.rscript),
            r_package_version("sva", context.rscript),
            r_package_version("kohonen", context.rscript),
        ]

    def run(
        self,
        dataset: Dataset,
        candidate_id: str,
        parameters: dict[str, Any],
        context: BackendContext,
    ) -> tuple[np.ndarray, CorrectionResult]:
        resolved = resolve_parameters(self.capabilities, parameters)
        software = self.availability(context)

        cycombine = next(s for s in software if s.name == "cyCombine")
        if not cycombine.available:
            raise BackendNotAvailable(
                "the cyCombine R package is not installed in this environment "
                f"({cycombine.detail}). Install it from the pinned commit recorded in "
                "config/upstream_lock.yaml; this adapter never substitutes another algorithm."
            )

        correction_markers = dataset.marker_ids("correction")
        if not correction_markers:
            raise BackendError(
                "no marker is flagged use_for_correction in the marker manifest"
            )

        obs = dataset.obs
        batches = obs[SAMPLE_ID].map(dataset.sample_batches()).astype(str).tolist()
        if any(b in {"None", "nan", ""} for b in batches):
            raise BackendError("at least one sample has no batch; cyCombine requires a batch column")

        covar = resolved.get("covar")
        covar_values = None
        if covar is not None:
            sample_covar = {
                sid: (getattr(rec, covar, None) or rec.technical_factors.get(covar))
                for sid, rec in dataset.samples.items()
            }
            if any(v is None for v in sample_covar.values()):
                raise BackendError(
                    f"covariate {covar!r} is requested but is missing for at least one "
                    "sample; a requested protected covariate is never quietly dropped"
                )
            covar_values = obs[SAMPLE_ID].map(sample_covar).astype(str).tolist()

        work = Path(context.work_dir)
        work.mkdir(parents=True, exist_ok=True)
        matrix_path = work / "input.f64"
        write_matrix(
            matrix_path,
            dataset.expression(features=correction_markers),
            row_ids=obs[CELL_ID].astype(str).tolist(),
            col_names=correction_markers,
        )

        # The anchor covariate is built on the R side from reference identity, so the
        # per-cell sample ids and the per-cell reference identity cross the boundary
        # here.  Only measured reference roles may anchor the correction; a biological
        # pseudo-control is never technical evidence.
        measured_roles = (ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR)
        reference_of = {
            sid: (rec.reference_id or "") if rec.reference_role in measured_roles else ""
            for sid, rec in dataset.samples.items()
        }
        reference_per_cell = obs[SAMPLE_ID].map(reference_of).fillna("").astype(str).tolist()

        params = {
            "input_matrix": str(matrix_path),
            "output_matrix": str(work / "corrected.f64"),
            "markers": correction_markers,
            "batch": batches,
            "sample": obs[SAMPLE_ID].astype(str).tolist(),
            "reference": reference_per_cell,
            "covar": covar_values,
            "covar_name": covar,
            "anchor": resolved.get("anchor"),
            "ref_batch": resolved.get("ref_batch"),
            "parametric": bool(resolved.get("parametric", True)),
            "norm_method": resolved.get("norm_method", "scale"),
            "cluster_method": resolved.get("cluster_method", "kohonen"),
            "distf": resolved.get("distf", "euclidean"),
            "xdim": int(resolved.get("xdim", 8)),
            "ydim": int(resolved.get("ydim", 8)),
            "rlen": int(resolved.get("rlen", 10)),
            "nClus": resolved.get("nClus"),
            "seed": int(context.seed),
            # Pinned by the adapter, never taken from user configuration.
            "method": "ComBat",
        }

        try:
            r_result = run_r_script(
                R_SCRIPT, params, work, rscript=context.rscript,
                log_path=work / "cycombine.log",
            )
        except RNotAvailable as exc:
            # No interpreter at all is an environment gap, not a method failure.
            raise BackendNotAvailable(str(exc)) from exc
        except RScriptFailed as exc:
            raise BackendError(
                f"cyCombine correction failed: {exc}. "
                f"Last lines of the R log:\n{exc.log[-2000:]}"
            ) from exc

        if r_result.payload.get("status") != "ok":
            raise BackendError(
                f"the cyCombine script did not report success: {r_result.payload.get('status')!r}"
            )

        corrected_block, row_ids, col_names, _ = read_matrix(params["output_matrix"])
        if row_ids != obs[CELL_ID].astype(str).tolist():
            raise BackendError(
                "cyCombine returned rows in a different order or with different cell ids; "
                "event identity could not be verified"
            )
        if col_names != correction_markers:
            raise BackendError(
                f"cyCombine returned columns {col_names}, expected {correction_markers}"
            )

        matrix = dataset.expression()
        matrix[:, dataset.feature_index(correction_markers)] = corrected_block

        checks = validate_output(dataset, matrix, correction_markers, candidate_id=candidate_id)
        labels = context.sampling.labels_by_cell if context.sampling else None
        by_marker, by_sample, by_population, changed, unchanged = correction_magnitude(
            dataset, matrix, labels=labels
        )

        spec = make_spec(candidate_id, "cycombine", cycombine.version, resolved,
                         context.seed, dataset, correction_markers)
        warnings = list(r_result.payload.get("warnings", []))
        if resolved.get("anchor") is not None:
            warnings.append(
                "the anchor route is explicitly experimental upstream: the "
                "correct_data documentation marks anchor 'Experimental' and warns it "
                "may be confounded with the condition. The anchor covariate was built "
                "from reference identity (cells of measured reference samples carry "
                "their reference id; every other cell carries its own sample id), and "
                "upstream silently ignores the anchor inside a label group that is "
                "confounded with batch or dominated by a single anchor."
            )
        for node in r_result.payload.get("confounded_nodes", []):
            warnings.append(
                f"SOM node {node} contains a single batch or a covariate level perfectly "
                "confounded with batch; the node contributes no estimable batch adjustment"
            )

        result = CorrectionResult(
            candidate=spec,
            status=CorrectionStatus.SUCCESS,
            output_scale=dataset.processing.scale,
            output_hash=matrix_hash(matrix),
            n_cells=dataset.n_cells,
            markers_changed=changed,
            markers_unchanged=unchanged,
            samples_changed=changed_samples(dataset, matrix),
            populations_changed=sorted(by_population),
            magnitude_by_marker=by_marker,
            magnitude_by_sample=by_sample,
            magnitude_by_population=by_population,
            validation=checks,
            warnings=warnings,
            software=software,
            provenance=[
                ProvenanceRef(kind="baseline", identifier=dataset.manifest.dataset_id,
                              content_hash=dataset.baseline_hash()),
                ProvenanceRef(kind="r_log", identifier=R_SCRIPT,
                              path=str(work / "cycombine.log")),
            ],
        )
        return matrix, result


register_backend(CyCombineBackend())
