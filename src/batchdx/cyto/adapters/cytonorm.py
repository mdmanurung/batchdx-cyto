"""CytoNorm correction backend — fit on one split, applied to the rest.

The audited API recorded in ``config/upstream_lock.yaml`` is
``CytoNorm.train(files, labels, channels, transformList, ...)`` followed by
``CytoNorm.normalize(model, files, ...)``, and both take **FCS file paths**.  That
shapes the whole adapter:

* the baseline is written out as one FCS file per sample, because that is the only
  thing CytoNorm will accept;
* the ``transformList`` is the **identity**, because the baseline is already on the
  declared transformed scale.  CytoNorm's own vignettes pass a ``flowCore``
  transform because they start from raw files; passing one here would transform the
  values twice, and the R script refuses any declaration other than ``identity``
  rather than ignoring it;
* the model is fitted on the **training split's** measured reference aliquots only,
  which is what makes the split a fit/apply split rather than a label;
* the normalised files are read back and every event is matched to the baseline
  before any value is accepted.

Event identity is verified with a dedicated channel.  ``CytoNorm.normalize``
reassigns each FlowSOM metacluster back into the slots it came from, so row order
should survive; that is exactly the kind of thing that has to be checked rather than
assumed, so a small integer identity channel is carried through both writes and
compared cell by cell against the baseline.  If the ordering cannot be verified the
run fails.  Returning the input unchanged and calling it a success is the one
outcome this adapter will not produce.

CytoNorm is GPL (>= 2) and is **called, never copied or vendored**
(``docs/decisions/0003-cytobatchflagr-integration.md``).  ``fit_apply_verified`` below
is a claim about a test that actually ran against the installed package; see the
evidence note on the constant.
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
from batchdx.cyto.io.dataset import EVENT_INDEX, SAMPLE_ID, Dataset
from batchdx.cyto.io.fcs import (
    FCSError,
    assert_fcs_round_trip,
    read_fcs,
    round_trip_tolerance,
    write_fcs,
)
from batchdx.cyto.io.hashing import matrix_hash
from batchdx.cyto.sampling.splits import (
    RECOMMENDATION_SPLIT,
    SplitError,
    missing_arms,
    split_members,
    split_unit_field,
)
from batchdx.cyto.schemas.contracts import CorrectionResult, ProvenanceRef, SoftwareVersion
from batchdx.cyto.schemas.enums import CorrectionStatus, ReferenceRole

R_SCRIPT = "cytonorm_correct.R"

#: Channel carrying the per-sample event index so that a corrected file can be
#: matched back to the baseline cell by cell.  It is written to disk but is never
#: named as a normalisation channel, so CytoNorm leaves it alone and it survives
#: byte-identically; that is what makes the identity check a real check.
IDENTITY_CHANNEL = "BATCHDX_EVENT_IDX"

#: The transform list handed to CytoNorm.  Not a free choice: the input is the
#: already-transformed preprocessed baseline, so the honest list is the identity.
#: Anything else would apply a second transform to values that already carry one.
TRANSFORM_LIST = "identity"

#: Reference roles that count as measured technical evidence.  A biological pseudo
#: control is not one, so it is never fitted on.
_MEASURED_REFERENCE_ROLES = {ReferenceRole.TECHNICAL_REPLICATE, ReferenceRole.ANCHOR}

#: Fit-and-apply verification.
#:
#: Set to True only after a real correction was executed against the installed
#: CytoNorm 2.0.12 at the pinned commit and the event-identity check passed on the
#: files it wrote.  The check is `tests/integration/test_cytonorm_backend.py::
#: test_the_cytonorm_backend_corrects_a_held_out_split_and_verifies_event_identity`,
#: which skips (visibly, and the Stage 6 gate records it) when the package is
#: absent.  If that test cannot run in an environment, this stays False and the
#: capability reads as declared-but-unverified rather than as working.
FIT_APPLY_VERIFIED = True


def _as_list(value: Any) -> list[str]:
    """Normalise a JSON scalar-or-array into a list of strings.

    ``jsonlite`` with ``auto_unbox = TRUE`` collapses a length-one vector to a bare
    string, so a payload field that is normally a list is sometimes a str.  Iterating
    that str would iterate its *characters*, and the batch name ``"B3"`` would become
    ``["3", "B"]`` -- a wrong answer that reads as right.  Every list-valued field from
    R goes through here.
    """

    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(v) for v in value]


def _unavailable_here(reason: str) -> BackendNotAvailable:
    """Type a design-coverage refusal as *not available*, not as a failure.

    These refusals are not defects: the algorithm is fine, the data simply does not
    contain the layout this backend declares it needs
    (``reference_layouts = measured_reference_in_every_batch``) inside the split it must
    fit on.  ``docs/decisions/0006-candidate-results-and-failure-handling.md`` requires
    "this method did not run here" to stay distinguishable from "this method ran and
    broke", and recording a coverage limitation as ``failed`` would tell a reader of a
    report that CytoNorm is broken on designs where it simply has nothing to fit on.

    ``BackendNotAvailable`` is used for that, with a reason that names the coverage
    property and not a missing package, so nothing here claims a package is absent when
    it is present.  ``pipeline.py`` has no third branch for an ineligibility exception
    type, so adding one is a required central change; until then this is the closest
    honest status the pipeline can record.
    """

    return BackendNotAvailable(
        "the layout this backend requires is not available in the split it must fit on: "
        + reason
    )


class CytoNormBackend:
    """Quantile normalisation fitted on held-in technical references, applied to all."""

    capabilities = BackendCapabilities(
        name="cytonorm",
        required_sample_fields=("batch", "reference_id", "reference_role"),
        reference_layouts=("measured_reference_in_every_batch",),
        supports_covariates=False,
        supports_marker_mask=True,
        input_scale="transformed, with the transform declared to CytoNorm as a transformList",
        output_scale="same transformed scale as the input",
        supports_fit_apply=True,
        fit_apply_verified=FIT_APPLY_VERIFIED,
        resource_profile=(
            "memory: per-file; the model is fitted from control files only"
        ),
        accepted_parameters={
            "nCells": 1_000_000,
            "xdim": 15,
            "ydim": 15,
            "nClus": 10,
            "nQ": 99,
            "scale": False,
        },
        rejected_parameters={
            "aggregate_control": (
                "an artificial aggregate control is outside the validated path and is "
                "not offered as a standard option"
            ),
            "transformList": (
                "the input is the already-transformed preprocessed baseline, so the "
                "transform list handed to CytoNorm is the identity and is pinned by the "
                "adapter. Declaring a flowCore transform here would transform the values "
                "a second time"
            ),
        },
        notes=(
            "fit and apply are separable in the audited API: the model is fitted on the "
            "training split's measured reference aliquots and applied to every sample, "
            "including the tuning and validation splits. The fit requires a measured "
            "reference in at least two batches inside the training split, and the backend "
            "reports ineligibility rather than fitting on one. CytoNorm leaves a cluster "
            "on the identity function when its per-batch quantiles are degenerate, and the "
            "count of such clusters is reported rather than hidden"
        ),
    )

    def availability(self, context: BackendContext) -> list[SoftwareVersion]:
        return [
            r_version(context.rscript),
            r_package_version("CytoNorm", context.rscript),
            r_package_version("FlowSOM", context.rscript),
            r_package_version("flowCore", context.rscript),
        ]

    # -- the split and the references ------------------------------------------------

    def _split(
        self, dataset: Dataset, context: BackendContext
    ) -> tuple[dict[str, list[str]], str | None]:
        if context.sampling is None or not context.sampling.splits:
            raise _unavailable_here(
                "CytoNorm is a fit-and-apply backend and needs a train/tune/validate "
                "partition to fit on, but no SamplingManifest with splits was supplied. "
                "Fit and apply cannot be separated without one"
            )
        splits = context.sampling.splits
        try:
            members = split_members(splits)
            unit = split_unit_field(splits)
        except SplitError as exc:
            raise BackendError(f"the frozen split cannot be read: {exc}") from exc
        # A fit-and-apply backend may not fit on a partition that cannot be compared:
        # the sampling manifest tolerates an arm with fewer biospecimens than there
        # are splits, because it also serves detectors and metrics, but fitting on such
        # a partition would produce a held-out number for a contrast that does not
        # exist inside the held-out data.
        gaps = missing_arms(splits, dataset.samples)
        if gaps:
            raise _unavailable_here(
                "the frozen train/tune/validate partition is not arm-balanced: "
                + "; ".join(f"split {name!r} lacks arm(s) {arms}" for name, arms in gaps.items())
                + ". CytoNorm may not fit on a partition whose validation split cannot hold "
                "one arm of the biological contrast, because the held-out comparison it "
                "would produce does not exist inside the held-out data. Declare more "
                "biospecimens in the thinner arm"
            )
        return members, unit

    def _training_controls(
        self, dataset: Dataset, split: dict[str, list[str]], unit_field: str | None
    ) -> tuple[list[str], dict[str, str]]:
        """Measured reference samples in the training split, by sample id.

        The references are what CytoNorm learns the per-batch quantiles from, so they
        must be on the *training* side.  A reference that ended up in the validation
        split cannot be used: the point of the split is that the fit never sees it.
        """

        train = set(split["train"])
        controls = {
            sid for sid, rec in dataset.samples.items()
            if sid in train and rec.reference_role in _MEASURED_REFERENCE_ROLES
        }
        if not controls:
            raise _unavailable_here(
                "no measured technical reference falls in the training split, so there is "
                "nothing for CytoNorm to fit on. CytoNorm learns per-batch quantiles from "
                "control aliquots; fitting on validation or tuning samples would make the "
                "held-out comparison meaningless, and fitting on a pseudo control is "
                "forbidden. Declare more measured reference aliquots, or accept that this "
                "candidate is ineligible"
            )
        batches = sorted({
            str(dataset.samples[sid].batch) for sid in controls if dataset.samples[sid].batch
        })
        if len(batches) < 2:
            raise _unavailable_here(
                f"the training split holds a measured reference in {len(batches)} batch(es) "
                f"({batches}); CytoNorm fits a per-batch quantile goal and needs at least "
                "two. This is a coverage failure of the split, not a defect in CytoNorm, and "
                "it is not worked around by fitting on a single batch"
            )
        if unit_field in {"donor_id", "biospecimen_id"}:
            # A control whose donor also contributes biological samples to the
            # validation split would leak the donor's biology into the fit.  Say so.
            donors = {
                str(getattr(dataset.samples[s], unit_field, "") or "") for s in controls
            }
            validation = set(split[RECOMMENDATION_SPLIT])
            leaked = sorted(
                sid for sid in validation
                if str(getattr(dataset.samples[sid], unit_field, "") or "") in donors
            )
            if leaked:
                raise _unavailable_here(
                    f"the training split's reference aliquot(s) {sorted(controls)} share a "
                    f"{unit_field} with validation-split sample(s) {leaked}; fitting on them "
                    "would put held-out biospecimens into the model. The split is keyed on "
                    f"{unit_field} for this reason and this check is not optional"
                )
        return sorted(controls), {sid: str(dataset.samples[sid].batch) for sid in controls}

    # -- the run ---------------------------------------------------------------------

    def run(
        self,
        dataset: Dataset,
        candidate_id: str,
        parameters: dict[str, Any],
        context: BackendContext,
    ) -> tuple[np.ndarray, CorrectionResult]:
        resolved = resolve_parameters(self.capabilities, parameters)
        software = self.availability(context)

        cytonorm = next(s for s in software if s.name == "CytoNorm")
        if not cytonorm.available:
            raise BackendNotAvailable(
                "the CytoNorm R package is not installed in this environment "
                f"({cytonorm.detail}). Install it at the commit pinned in "
                "config/upstream_lock.yaml; this adapter never substitutes another "
                "algorithm, and CytoNorm is GPL (>= 2) so it is called, never vendored"
            )

        markers = dataset.marker_ids("correction")
        if not markers:
            raise BackendError("no marker is flagged use_for_correction in the marker manifest")

        split, unit_field = self._split(dataset, context)
        controls, control_batches = self._training_controls(dataset, split, unit_field)

        obs = dataset.obs
        sample_ids = sorted(dataset.samples)
        work = Path(context.work_dir)
        in_dir = work / "fcs_in"
        work.mkdir(parents=True, exist_ok=True)

        # The identity channel must not collide with a real marker.
        clash = [m for m in markers if m.upper() == IDENTITY_CHANNEL.upper()]
        if clash:
            raise BackendError(
                f"the marker {clash[0]!r} collides with the {IDENTITY_CHANNEL!r} identity "
                "channel CytoNorm's event-identity check requires; rename the marker"
            )

        # One FCS file per sample, in the sample's own acquisition order, carrying its
        # original event index so the corrected file can be matched back cell by cell.
        # ``rows_of`` fixes that ordering once; the write, the identity check and the
        # scatter all read it, so they cannot drift apart.
        values = dataset.expression()
        feature_pos = dataset.feature_index(markers)
        rows_of = {sid: _rows_in_acquisition_order(dataset, sid) for sid in sample_ids}
        paths: dict[str, Path] = {}
        write_error = 0.0
        for sid in sample_ids:
            rows = rows_of[sid]
            index = dataset.obs[EVENT_INDEX].to_numpy().astype(np.int64)[rows]
            events = np.column_stack([index.astype(np.float64), values[rows][:, feature_pos]])
            try:
                paths[sid] = write_fcs(
                    in_dir / f"{sid}.fcs",
                    events,
                    [IDENTITY_CHANNEL, *markers],
                    [IDENTITY_CHANNEL, *markers],
                    scale=dataset.processing.scale,
                    processing_note=(
                        f"batchdx-cyto baseline for CytoNorm; already transformed "
                        f"({dataset.processing.scale}); the transform list handed to "
                        f"CytoNorm is the identity"
                    ),
                    source_sample=sid,
                    candidate_id=candidate_id,
                    overwrite=True,
                )
                # Prove the write itself lost nothing beyond what the format cannot
                # hold, so any later difference is CytoNorm's and not ours.  Checked
                # on the *input* file: the corrected file is supposed to differ.
                write_error = max(write_error, assert_fcs_round_trip(
                    events, read_fcs(paths[sid]), [IDENTITY_CHANNEL, *markers],
                    context=f"sample {sid} written for CytoNorm",
                ))
            except FCSError as exc:
                raise BackendError(
                    f"the per-sample FCS CytoNorm requires is not faithful to the baseline "
                    f"for sample {sid!r}: {exc}"
                ) from exc

        train_files = [str(paths[sid]) for sid in controls]
        train_labels = [control_batches[sid] for sid in controls]
        apply_files = [str(paths[sid]) for sid in sample_ids]
        apply_labels = [str(dataset.samples[sid].batch) for sid in sample_ids]

        params = {
            "train_files": train_files,
            "train_labels": train_labels,
            "apply_files": apply_files,
            "apply_labels": apply_labels,
            "channels": markers,
            "identity_channel": IDENTITY_CHANNEL,
            "transform_list": TRANSFORM_LIST,
            "train_dir": str(work / "cytonorm_train"),
            "output_dir": str(work / "cytonorm_out"),
            "prefix": "Norm_",
            "seed": int(context.seed),
            "nCells": int(resolved["nCells"]),
            "xdim": int(resolved["xdim"]),
            "ydim": int(resolved["ydim"]),
            "nClus": int(resolved["nClus"]),
            "nQ": int(resolved["nQ"]),
            "scale": bool(resolved["scale"]),
        }

        try:
            r_result = run_r_script(
                R_SCRIPT, params, work, rscript=context.rscript,
                log_path=work / "cytonorm.log",
            )
        except RNotAvailable as exc:
            raise BackendNotAvailable(str(exc)) from exc
        except RScriptFailed as exc:
            raise BackendError(
                f"CytoNorm correction failed (exit {exc.returncode}). "
                f"Last lines of the R log:\n{exc.log[-2000:]}"
            ) from exc

        # Read the corrected files back and check every event against the baseline
        # before a single value is accepted.  A missing output is a failure, never a
        # silent pass-through of the input.
        payload = r_result.payload
        if payload.get("status") != "ok":
            raise BackendError(
                f"the CytoNorm script did not report success: {payload.get('status')!r}"
            )

        out_dir = Path(params["output_dir"])
        cache: dict[str, Any] = {}

        def contents_of(sid: str):
            """Read one corrected file, once, and check its channel layout."""

            if sid in cache:
                return cache[sid]
            path = out_dir / f"{params['prefix']}{sid}.fcs"
            if not path.exists():
                raise BackendError(
                    f"CytoNorm wrote no corrected file for sample {sid!r}; expected {path}. "
                    "A missing file is not an unchanged file"
                )
            try:
                contents = read_fcs(path)
            except FCSError as exc:
                raise BackendError(
                    f"the file CytoNorm wrote for {sid!r} could not be read back: {exc}"
                ) from exc
            if contents.channel_names != [IDENTITY_CHANNEL, *markers]:
                raise BackendError(
                    f"CytoNorm's output for {sid!r} has channels {contents.channel_names}, "
                    f"expected {[IDENTITY_CHANNEL, *markers]}; the channel to feature mapping "
                    "could not be verified"
                )
            cache[sid] = contents
            return contents

        matrix = dataset.expression()
        max_abs_change = 0.0
        uncorrected: list[str] = []
        for sid in sample_ids:
            rows = rows_of[sid]
            index = obs[EVENT_INDEX].to_numpy().astype(np.int64)[rows]
            contents = contents_of(sid)
            verify_event_identity(contents, index, sid)
            corrected = contents.events[:, 1:]
            if not np.isfinite(corrected).all():
                bad = int(np.argmax(~np.isfinite(corrected).all(axis=1)))
                raise BackendError(
                    f"CytoNorm's output for sample {sid!r} contains a non-finite value at "
                    f"event {int(index[bad])}; non-finite values are not repaired or replaced"
                )
            delta = float(np.max(np.abs(corrected - values[rows][:, feature_pos])))
            if delta <= round_trip_tolerance(values[rows][:, feature_pos]):
                # CytoNorm writes a file back untouched when it has no model for that
                # batch, so "wrote a file" and "corrected the cells" are not the same
                # claim.  The comparison is against the single-precision bound rather
                # than exact equality, because a file the backend never touched still
                # comes back through float32 FCS storage and is not bit-identical to a
                # float64 baseline.  A change smaller than what the format can represent
                # is not a correction.
                uncorrected.append(sid)
            max_abs_change = max(max_abs_change, delta)
            # np.ix_, not matrix[rows, cols]: the latter takes the *outer* product of
            # the two index arrays, so it would demand an (n_rows, n_markers) index grid
            # and silently misbehave for anything but a full block.
            matrix[np.ix_(rows, feature_pos)] = corrected
            dataset.assert_baseline_unchanged()

        checks = validate_output(dataset, matrix, markers, candidate_id=candidate_id)
        labels = context.sampling.labels_by_cell if context.sampling else None
        by_marker, by_sample, by_population, changed, unchanged = correction_magnitude(
            dataset, matrix, labels=labels
        )

        spec = make_spec(candidate_id, "cytonorm", cytonorm.version, resolved,
                         context.seed, dataset, markers)

        warnings_out: list[str] = []
        unfitted = sorted(_as_list(payload.get("unfitted_labels")))
        if unfitted:
            warnings_out.append(
                f"PARTIAL CORRECTION: no quantile model was fitted for batch(es) {unfitted}, "
                f"because the training split held no measured reference there. CytoNorm "
                f"writes those files back untouched, so {len(uncorrected)} of "
                f"{len(sample_ids)} sample(s) came back identical to the baseline. A batch "
                "contrast involving those batches is not corrected by this candidate, and "
                "the candidate must not be read as if it were"
            )
        if uncorrected:
            batches_left = sorted({str(dataset.samples[s].batch) for s in uncorrected})
            warnings_out.append(
                f"sample(s) {uncorrected} were returned by CytoNorm bit-identical to the "
                f"baseline (batch(es) {batches_left}); they are counted as uncorrected, not as "
                "a successful correction of unchanged data"
            )
        identity_clusters = int(payload.get("n_identity_function_clusters") or 0)
        if identity_clusters:
            warnings_out.append(
                f"CytoNorm left {identity_clusters} FlowSOM metacluster(s) on the identity "
                "function because their per-batch quantiles were degenerate; events in those "
                "clusters come back unchanged. This is a partial correction and is reported "
                "as one, not as a complete one"
            )
        seen: set[str] = set()
        for message in payload.get("warnings", []) or []:
            if message in seen:
                continue
            seen.add(str(message))
            warnings_out.append(f"CytoNorm: {message}")
        warnings_out.append(
            "the model was fitted on the training split's measured reference aliquots only: "
            + ", ".join(f"{sid} ({control_batches[sid]})" for sid in controls)
            + f"; it was then applied to every sample, including the {RECOMMENDATION_SPLIT!r} "
            "split, so no validation sample contributed to the fit"
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
            validation={**checks, "all_samples_corrected": not uncorrected},
            warnings=warnings_out,
            software=software,
            resources={
                "r_elapsed_seconds": float(payload.get("elapsed_seconds", 0.0)),
                # The observed FCS single-precision loss on the *input* write, and the
                # observed size of the correction.  Both are reported rather than
                # assumed to be zero.
                "fcs_write_max_abs_error": write_error,
                "max_abs_change": max_abs_change,
                "n_metaclusters": float(payload.get("n_metaclusters", 0) or 0),
                "n_fitted_batches": float(payload.get("n_fitted_batches", 0) or 0),
                "n_uncorrected_samples": float(len(uncorrected)),
                "n_housekeeping_warnings": float(payload.get("n_housekeeping_warnings", 0) or 0),
            },
            provenance=[
                ProvenanceRef(kind="baseline", identifier=dataset.manifest.dataset_id,
                              content_hash=dataset.baseline_hash()),
                ProvenanceRef(kind="r_log", identifier=R_SCRIPT,
                              path=str(work / "cytonorm.log")),
            ],
        )
        return matrix, result


class EventIdentityError(BackendError):
    """The corrected file could not be matched to the baseline it corrects.

    Distinct from a plain :class:`BackendError` only in the message it carries, but
    the distinction is load-bearing: this is not "CytoNorm ran and produced
    something", it is "what CytoNorm produced is not a correction *of these cells*",
    and the package will not score the first as if it were the second.
    """


def verify_event_identity(
    contents, expected_index: np.ndarray, sample_id: str
) -> np.ndarray:
    """Require a corrected FCS file to be the baseline, event for event and in order.

    ``CytoNorm.normalize`` splits each file by FlowSOM metacluster and writes the
    corrected values back into the slots the cluster came from, so row order ought
    to survive.  "Ought to" is not a verification: this function compares the
    identity channel against the baseline's own event indices and requires the two
    ordered sequences to be *equal*.  Anything else — a shuffle, a dropped event, a
    duplicated one, an extra channel ahead of the identity channel — raises, because
    scattering the corrected values onto the wrong cells is a silent, total failure
    that no downstream check would catch.

    Returns the identity channel as integer event indices, so the caller can also
    report what it saw.
    """

    expected = np.asarray(expected_index, dtype=np.int64)
    if contents.n_events != expected.size:
        raise EventIdentityError(
            f"CytoNorm's output for sample {sample_id!r} holds {contents.n_events} events, "
            f"expected {expected.size}; event identity could not be verified and the file is "
            "refused rather than padded or truncated to fit"
        )
    written = np.rint(np.asarray(contents.events[:, 0], dtype=np.float64)).astype(np.int64)
    if not np.array_equal(written, expected):
        first = int(np.argmax(written != expected)) if written.size else 0
        raise EventIdentityError(
            f"CytoNorm's output for sample {sample_id!r} does not correspond event by event "
            f"to the baseline: the identity channel first diverges at row {first} "
            f"(got event index {int(written[first]) if written.size else 'n/a'}, expected "
            f"{int(expected[first]) if expected.size else 'n/a'}). The corrected file is "
            "refused rather than matched approximately, because a wrong cell order is a "
            "silent total failure rather than a small numerical one"
        )
    return written


def _rows_in_acquisition_order(dataset: Dataset, sample_id: str) -> np.ndarray:
    """Baseline row positions of one sample's events, in acquisition order.

    The FCS write, the identity check and the scatter back onto the baseline must
    all agree on which row is which event.  Defining the ordering once is the only
    way to stop them from drifting apart; three inline copies of the same sort is
    three chances to disagree about cell identity.
    """

    obs = dataset.obs
    # Positional, not label-based: ``Dataset.expression`` addresses the matrix by
    # position, and obs_names are cell ids whose dtype is not guaranteed to be an
    # integer index.  ``flatnonzero`` keeps the two addressing schemes identical.
    rows = np.flatnonzero(obs[SAMPLE_ID].to_numpy() == sample_id)
    index = obs[EVENT_INDEX].to_numpy().astype(np.int64)[rows]
    return rows[np.argsort(index, kind="stable")]


register_backend(CytoNormBackend())
