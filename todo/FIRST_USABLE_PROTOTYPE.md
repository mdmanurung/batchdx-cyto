# First usable batchdx-cyto prototype

Date: 2026-10-10. Status: P01/P02 implemented and verified; P03–P20 pending.

The detailed, source-grounded checklist is
[PROTOTYPE_TASKS.md](../review/2026-10-09/PROTOTYPE_TASKS.md).
It contains 20 prototype tasks, 20 conditional extension tasks, and five broader
capability choices. Preserve the existing mixed worktree when implementing.

1. P01–P05: contain failures, validate R payloads, fix DAG dependencies, escape and
   redact reports, refuse unsafe init and incompatible records.
2. P06–P10: implement and measure all nine pinned OpenProblems scores, including
   FlowSOM and both functional-marker metrics; report the complete roster and
   prevent missing required evidence from improving an assessment. P08 includes
   OpenProblems/scIB-style rank visuals, benchmark controls/scaling, normalization
   QC, and offline HTML/vector/data exports; see
   [BENCHMARK_VISUALIZATION_SPEC.md](../review/2026-10-09/BENCHMARK_VISUALIZATION_SPEC.md).
3. P11–P15: honor protection and marker masks, execute both benchmark splits, expose
   correction applicability, and disclose diagnostic support and threshold limits.
4. P16–P20: prepare and run one measured dataset, test real methods through the
   workflow, verify installation, and reconcile docs/status with complete receipts.

Proposed scope: local, preprocessed H5AD, one panel, identity plus cyCombine,
exploratory evidence, with the full nine-score OpenProblems suite required by
the user on 2026-10-10. Custom/common/native metrics are supplemental.
The user also requires the OpenProblems/scIB visualization style: method rows,
nine grouped score columns, category/overall bars, explicit ranks and failures.
Use the OpenProblems-style control-scaled mean, not implicit scIB 40/60 weights.
CytoNorm/raw-FCS/stronger science/HPC have separate extension
gates. No new automatic recommendation or publication claim is authorized here.

Acceptance: a clean install reproduces a full-suite eligible measured benchmark run and
inspectable corrected H5AD/report; failures and unsupported evidence stay visible;
the measured run includes rank visuals and normalization QC with reproducible
CSV/JSON/SVG/PDF/PNG exports;
no-op/input-change/recovery behavior is verified; scoped tests and documentation
agree. All nine scores must execute on eligible inputs; the control-only spectral
subset cannot alone satisfy this. Estimated 13–20 engineer-days, excluding
access/preparation/infrastructure delays.

Completed P01: 23 new candidate/config regressions plus a simulated missing-output
Snakemake run pass. Final config/candidate subset: 39 passed. Broader checks:
307 passed, 1 skipped, 1 failed, 5 deselected; the CytoNorm availability failure is
reproduced on untouched HEAD. See
[P01 evidence](../review/2026-10-09/P01_FAILURE_CONTAINMENT.md).

Completed P02: checked R payload/string-array boundary, shared list preservation,
bounded execution/probes, partial timeout logs and caller failure reasons. Final
focused tests: 40 passed. R4_51 exchange: 110 assertions passed, no skips. Related
checks: 93 passed, 3 failed, 25 deselected; all three CytoNorm failures reproduce on
untouched HEAD under the same default PATH. See
[P02 evidence](../review/2026-10-09/P02_R_BOUNDARY.md).
The updated simulated missing-matrix Snakemake regression also passed: 1 passed,
14 deselected. These checks do not execute real correction algorithms.

Next bounded task: P03, workflow input dependencies and conditional candidate
output bundles. Start with changed-input and missing-H5AD regressions. Full-suite
metrics and rank visuals remain required; the full package gate is not green.
