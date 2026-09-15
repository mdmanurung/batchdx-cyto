# Implementation status

**This file is generated. Do not edit it by hand.**

Derived by `scripts/generate_status.py` from `report.json` on 2026-09-15 06:48 UTC.
Every line below corresponds to test identifiers that were executed; a claim
without an executed test cannot appear here.

Regenerate with:

```
pytest --json-report --json-report-file=build/report.json
python scripts/generate_status.py --report build/report.json
```

## Legend

| Badge | Meaning |
|---|---|
| ✅ PASSED | every test for the gate executed and passed |
| 🟡 PARTIAL | the gate is met in part; the omission is stated |
| ⛔ BLOCKED | an external obstacle prevents part of the gate from being tested at all |
| ❌ FAILED | at least one test for the gate failed |
| ⬜ DEFERRED | not implemented in this build |
| ❔ NOT RUN | a named test was not present in this report |

## Summary by tier

| Tier | Stage | Title | Status | Tests passed |
|---|---|---|---|---|
| A | X | Cross-stage contracts (typed, versioned, round-tripped across the R boundary) | ✅ PASSED | 10/10 |
| A | 0 | Audit upstream resources and establish contracts | ✅ PASSED | 8/8 |
| A | 1 | Runnable package and minimal Snakemake workflow | ✅ PASSED | 6/6 |
| A | 2 | Ingestion, immutable data, optional preprocessing | 🟡 PARTIAL | 16/16 |
| A | 3 | Design audit and correction eligibility | ✅ PASSED | 12/12 |
| A | 5 | Detection adapters and diagnosis report | ⛔ BLOCKED | 10/11 |
| A | 6 | Correction backends and output validation | ⛔ BLOCKED | 14/15 |
| A | 7 | Native and common comparative metrics | ⛔ BLOCKED | 14/14 |
| A | 10 | Reporting and safe exports | 🟡 PARTIAL | 10/10 |
| B | 4 | Freeze evaluation populations, sampling and splits | 🟡 PARTIAL | 10/10 |
| B | 8 | Biological preservation and uncertainty | 🟡 PARTIAL | 5/5 |
| B | 9 | Candidate assessment and recommendation rules | ✅ PASSED | 12/12 |
| C | 11 | Snakemake execution, performance, deployment | 🟡 PARTIAL | 11/11 |
| C | 12 | Scientific validation, documentation, release readiness | ⬜ DEFERRED | 0/0 |

## Gate detail

### Stage X (Tier A) — Cross-stage contracts (typed, versioned, round-tripped across the R boundary)

**✅ PASSED** — 10 passed, 0 failed, 0 skipped, 0 not executed.

<details><summary>All test identifiers for this gate</summary>

* `tests/unit/test_contracts.py::test_all_seven_contracts_exist_and_are_versioned` — passed
* `tests/unit/test_contracts.py::test_contracts_reject_unknown_keys` — passed
* `tests/unit/test_contracts.py::test_contracts_round_trip_through_json` — passed
* `tests/unit/test_contracts.py::test_metric_record_carries_every_field_the_contract_requires` — passed
* `tests/unit/test_contracts.py::test_correction_record_distinguishes_success_ineligible_and_failed` — passed
* `tests/unit/test_contracts.py::test_json_schemas_are_exported_for_the_r_side` — passed
* `tests/unit/test_contracts.py::test_r_reads_every_contract_and_the_binary_matrix_format` — passed
* `tests/unit/test_contracts.py::test_binary_matrix_rejects_a_shape_mismatch` — passed
* `tests/unit/test_contracts.py::test_downstream_never_infers_success_from_a_directory` — passed
* `tests/unit/test_contracts.py::test_contract_schema_version_is_stamped_on_every_record` — passed

</details>

### Stage 0 (Tier A) — Audit upstream resources and establish contracts

**✅ PASSED** — 8 passed, 0 failed, 0 skipped, 0 not executed.

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage0_audit.py::test_lock_lists_a_commit_licence_and_strategy_for_every_resource` — passed
* `tests/integration/test_stage0_audit.py::test_fetch_script_pins_the_same_commits_as_the_lock` — passed
* `tests/integration/test_stage0_audit.py::test_python_audit_runs_and_reconciles_commits_and_metrics` — passed
* `tests/integration/test_stage0_audit.py::test_python_audit_fails_on_a_commit_mismatch` — passed
* `tests/integration/test_stage0_audit.py::test_r_audit_runs_and_finds_every_adapter_referenced_function` — passed
* `tests/integration/test_stage0_audit.py::test_audit_output_documents_the_facts_the_adapters_rely_on` — passed
* `tests/integration/test_stage0_audit.py::test_cytobatchflagr_licence_forbids_vendoring_and_the_adapter_says_so` — passed
* `tests/integration/test_stage0_audit.py::test_every_adapter_referenced_r_function_appears_in_the_lock` — passed

</details>

### Stage 1 (Tier A) — Runnable package and minimal Snakemake workflow

**✅ PASSED** — 6 passed, 0 failed, 0 skipped, 0 not executed.

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage1_slice.py::test_cli_module_entry_point_is_importable_and_reports_what_it_ships` — passed
* `tests/integration/test_stage1_slice.py::test_example_initialises_every_file_it_references` — passed
* `tests/integration/test_stage1_slice.py::test_snakemake_runs_the_slice_and_a_second_invocation_does_no_analytical_work` — passed
* `tests/integration/test_stage1_slice.py::test_a_bad_configuration_propagates_a_useful_error_and_exit_code` — passed
* `tests/integration/test_stage1_slice.py::test_correct_requires_an_explicitly_selected_candidate` — passed
* `tests/integration/test_stage1_slice.py::test_packaged_workflow_is_found_from_outside_the_source_checkout` — passed

</details>

### Stage 2 (Tier A) — Ingestion, immutable data, optional preprocessing

**🟡 PARTIAL** — 16 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** SPiSCy-style preprocessing (declared gating, PeacoQC, required technical channel validation, per-removed-event QC tracking) is NOT implemented. The configuration schema declares it and refuses an undeclared cofactor, and already-QCed data bypasses preprocessing, which is the path Tier A uses.

<details><summary>All test identifiers for this gate</summary>

* `tests/unit/test_stage2_io.py::test_h5ad_roundtrip_preserves_marker_and_channel_order` — passed
* `tests/unit/test_stage2_io.py::test_cell_identity_is_stable_and_unique` — passed
* `tests/unit/test_stage2_io.py::test_baseline_is_immutable_under_read_operations` — passed
* `tests/unit/test_stage2_io.py::test_duplicate_sample_ids_are_rejected` — passed
* `tests/unit/test_stage2_io.py::test_ambiguous_channel_mapping_is_rejected` — passed
* `tests/unit/test_stage2_io.py::test_undeclared_scale_is_rejected` — passed
* `tests/unit/test_stage2_io.py::test_double_transformation_is_rejected` — passed
* `tests/unit/test_stage2_io.py::test_already_corrected_input_is_refused` — passed
* `tests/unit/test_stage2_io.py::test_missing_marker_is_not_treated_as_zero` — passed
* `tests/unit/test_stage2_io.py::test_nonfinite_values_are_rejected` — passed
* `tests/unit/test_stage2_io.py::test_missing_expression_layer_names_what_is_available` — passed
* `tests/unit/test_stage2_io.py::test_multiple_panels_must_be_run_separately` — passed
* `tests/unit/test_stage2_io.py::test_qc_history_and_transformation_provenance_are_carried` — passed
* `tests/unit/test_stage2_io.py::test_fcs_write_read_roundtrip_and_keyword_policy` — passed
* `tests/unit/test_stage2_io.py::test_fcs_export_refuses_to_overwrite_a_source_file` — passed
* `tests/unit/test_stage2_io.py::test_uninvertible_transformation_keeps_transformed_values_and_says_so` — passed

</details>

### Stage 3 (Tier A) — Design audit and correction eligibility

**✅ PASSED** — 12 passed, 0 failed, 0 skipped, 0 not executed.

<details><summary>All test identifiers for this gate</summary>

* `tests/scientific/test_stage3_design.py::test_balanced_design_with_references_is_reference_supported` — passed
* `tests/scientific/test_stage3_design.py::test_partial_reference_overlap_is_reported_per_pair` — passed
* `tests/scientific/test_stage3_design.py::test_perfect_confounding_is_non_identifiable_and_blocks_correction` — passed
* `tests/scientific/test_stage3_design.py::test_two_batch_design_makes_the_cluster_detector_ineligible` — passed
* `tests/scientific/test_stage3_design.py::test_thirty_one_batches_exceeds_the_cluster_detector_range` — passed
* `tests/scientific/test_stage3_design.py::test_cytobatchflagr_needs_controls_in_every_batch_and_four_batches` — passed
* `tests/scientific/test_stage3_design.py::test_iqr_fence_cannot_flag_with_three_values` — passed
* `tests/scientific/test_stage3_design.py::test_repeated_stimulation_and_timepoints_do_not_become_technical_replicates` — passed
* `tests/scientific/test_stage3_design.py::test_missing_donor_ids_permit_descriptive_diagnosis_only` — passed
* `tests/scientific/test_stage3_design.py::test_population_specific_lack_of_support_is_reported` — passed
* `tests/scientific/test_stage3_design.py::test_nested_technical_structure_is_represented_and_the_target_named` — passed
* `tests/scientific/test_stage3_design.py::test_below_floor_is_a_configuration_value_with_a_declared_origin` — passed

</details>

### Stage 4 (Tier B) — Freeze evaluation populations, sampling and splits

**🟡 PARTIAL** — 10 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** Training / tuning / validation splits at the biospecimen or donor level are NOT implemented. They matter only for a backend with a verified fit-and-apply interface, and no such backend is executable in this build; the manifest carries an empty `splits` field rather than a fabricated one.

<details><summary>All test identifiers for this gate</summary>

* `tests/unit/test_stage4_partitions.py::test_sampling_manifest_is_reproducible` — passed
* `tests/unit/test_stage4_partitions.py::test_sampling_is_balanced_at_the_sample_level_and_reports_coverage` — passed
* `tests/unit/test_stage4_partitions.py::test_gating_labels_are_medium_confidence_and_recorded` — passed
* `tests/unit/test_stage4_partitions.py::test_clustering_labels_are_flagged_lowest_confidence` — passed
* `tests/unit/test_stage4_partitions.py::test_unmatched_cells_are_labelled_not_dropped` — passed
* `tests/unit/test_stage4_partitions.py::test_rare_population_generates_a_coverage_warning` — passed
* `tests/unit/test_stage4_partitions.py::test_candidate_order_does_not_change_labels_or_comparisons` — passed
* `tests/unit/test_stage4_partitions.py::test_every_metric_record_carries_label_tier_and_representation` — passed
* `tests/unit/test_stage4_partitions.py::test_non_identifiable_comparisons_are_excluded_from_the_frozen_pairs` — passed
* `tests/unit/test_stage4_partitions.py::test_sampling_does_not_redefine_population_frequency` — passed

</details>

### Stage 5 (Tier A) — Detection adapters and diagnosis report

**⛔ BLOCKED** — 10 passed, 0 failed, 1 skipped, 0 not executed.

> **Blocked:** Neither native detection system can be EXECUTED here: cyCombine needs sva (Bioconductor) and kohonen (CRAN), and CytoBatchFlagR needs a dozen CRAN packages; the development environment's egress policy blocks CRAN, Bioconductor and conda mirrors (docs/blocked.md). The adapters, their eligibility rules and their unavailable paths are tested; the "outputs match a direct upstream call within a documented tolerance" clause is NOT satisfied.

Tests needing attention:

* `tests/integration/test_stage5_detectors.py::test_native_cycombine_detector_matches_a_direct_upstream_call — skipped`

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage5_detectors.py::test_drift_detector_finds_an_injected_within_file_trend` — passed
* `tests/integration/test_stage5_detectors.py::test_drift_detector_is_ineligible_without_an_acquisition_order` — passed
* `tests/integration/test_stage5_detectors.py::test_composition_detector_separates_frequency_from_intensity` — passed
* `tests/integration/test_stage5_detectors.py::test_composition_detector_reports_insufficient_data_rather_than_a_p_value` — passed
* `tests/integration/test_stage5_detectors.py::test_two_batch_fixture_makes_the_cluster_detector_ineligible_not_uninformative` — passed
* `tests/integration/test_stage5_detectors.py::test_absent_r_package_yields_not_available_never_an_absence_of_batch_effects` — passed
* `tests/integration/test_stage5_detectors.py::test_cytobatchflagr_without_a_checkout_is_not_available_with_a_stated_remedy` — passed
* `tests/integration/test_stage5_detectors.py::test_cytobatchflagr_refuses_to_choose_thresholds_itself` — passed
* `tests/integration/test_stage5_detectors.py::test_native_cycombine_detector_matches_a_direct_upstream_call` — skipped
* `tests/integration/test_stage5_detectors.py::test_detector_results_round_trip_through_json` — passed
* `tests/integration/test_stage5_detectors.py::test_flag_counts_are_not_presented_as_probabilities` — passed

</details>

### Stage 6 (Tier A) — Correction backends and output validation

**⛔ BLOCKED** — 14 passed, 0 failed, 1 skipped, 0 not executed.

> **Blocked:** The cycombine backend cannot be smoke-tested here: the package and its Bioconductor dependencies are not installable in this environment (docs/blocked.md). The identity backend, the parameter contract, the count-model refusal, output validation and the failure path are all tested. cytonorm is declared and deliberately deferred.

Tests needing attention:

* `tests/integration/test_stage6_backends.py::test_cycombine_smoke_test — skipped`

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage6_backends.py::test_identity_backend_is_numerically_unchanged` — passed
* `tests/integration/test_stage6_backends.py::test_every_backend_declares_its_capabilities_explicitly` — passed
* `tests/integration/test_stage6_backends.py::test_cycombine_refuses_the_count_data_route` — passed
* `tests/integration/test_stage6_backends.py::test_the_r_script_would_refuse_a_non_combat_method` — passed
* `tests/integration/test_stage6_backends.py::test_unsupported_arguments_are_rejected_not_ignored` — passed
* `tests/integration/test_stage6_backends.py::test_a_requested_protected_covariate_is_never_quietly_dropped` — passed
* `tests/integration/test_stage6_backends.py::test_absent_backend_reports_not_available_with_an_actionable_reason` — passed
* `tests/integration/test_stage6_backends.py::test_cytonorm_is_declared_but_explicitly_deferred` — passed
* `tests/integration/test_stage6_backends.py::test_cytonorm_refuses_an_artificial_aggregate_control` — passed
* `tests/integration/test_stage6_backends.py::test_output_validation_catches_a_changed_excluded_marker` — passed
* `tests/integration/test_stage6_backends.py::test_output_validation_catches_shape_and_finiteness_problems` — passed
* `tests/integration/test_stage6_backends.py::test_a_runtime_failure_never_returns_unchanged_data_as_success` — passed
* `tests/integration/test_stage6_backends.py::test_correction_magnitude_is_reported_by_marker_sample_and_population` — passed
* `tests/integration/test_stage6_backends.py::test_candidate_identity_is_stable_under_parameter_reordering` — passed
* `tests/integration/test_stage6_backends.py::test_cycombine_smoke_test` — skipped

</details>

### Stage 7 (Tier A) — Native and common comparative metrics

**⛔ BLOCKED** — 14 passed, 0 failed, 0 skipped, 0 not executed.

> **Blocked:** The NATIVE namespace (cyCombine::evaluate_emd and evaluate_mad) is registered and documented but cannot be executed here, so the "compatibility components match pinned implementations within documented tolerances" clause is NOT satisfied. The common-namespace metrics, including the required joint-structure row, are implemented and tested. Every OpenProblems metric row remains deferred.

<details><summary>All test identifiers for this gate</summary>

* `tests/scientific/test_stage7_metrics.py::test_identical_distributions_give_near_zero_emd` — passed
* `tests/scientific/test_stage7_metrics.py::test_a_known_shift_moves_emd_by_the_size_of_the_shift` — passed
* `tests/scientific/test_stage7_metrics.py::test_scale_change_alone_is_visible_in_dispersion_not_only_in_emd` — passed
* `tests/scientific/test_stage7_metrics.py::test_constant_marker_yields_undefined_not_zero` — passed
* `tests/scientific/test_stage7_metrics.py::test_rare_population_returns_insufficient_data_with_a_reason` — passed
* `tests/scientific/test_stage7_metrics.py::test_unequal_sample_sizes_do_not_change_the_weighting` — passed
* `tests/scientific/test_stage7_metrics.py::test_zero_denominator_makes_the_reduction_undefined_not_infinite` — passed
* `tests/scientific/test_stage7_metrics.py::test_composition_only_change_moves_composition_and_not_within_population_emd` — passed
* `tests/scientific/test_stage7_metrics.py::test_mad_uses_the_normal_consistency_constant` — passed
* `tests/scientific/test_stage7_metrics.py::test_empirical_and_binned_emd_are_named_differently` — passed
* `tests/scientific/test_stage7_metrics.py::test_no_local_metric_claims_to_be_an_upstream_component` — passed
* `tests/scientific/test_stage7_metrics.py::test_unsupported_metrics_report_reasons_rather_than_zero_scores` — passed
* `tests/scientific/test_stage7_metrics.py::test_missing_metrics_cannot_improve_a_candidate` — passed
* `tests/scientific/test_stage7_metrics.py::test_sampling_manifest_round_trips_through_json` — passed

</details>

### Stage 8 (Tier B) — Biological preservation and uncertainty

**🟡 PARTIAL** — 5 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** Prespecified biological contrasts, paired and repeated-measure models, donor-block bootstrap and restricted permutation, and the OpenProblems functional-marker compatibility component are NOT implemented. Three of the four clauses of the Stage 8 gate are met by the joint-structure metric and the distortion guardrails; the "synthetic biological effect survives correction within predeclared tolerances" clause is only partially met, because no biological contrast model exists to state a tolerance against.

<details><summary>All test identifiers for this gate</summary>

* `tests/scientific/test_stage8_preservation.py::test_a_correct_offset_removal_improves_alignment_without_distortion` — passed
* `tests/scientific/test_stage8_preservation.py::test_a_collapsing_correction_is_rejected_despite_perfect_alignment` — passed
* `tests/scientific/test_stage8_preservation.py::test_marginal_preserving_joint_destroying_corruption_is_detected` — passed
* `tests/scientific/test_stage8_preservation.py::test_biology_only_differences_are_not_labelled_technical` — passed
* `tests/scientific/test_stage8_preservation.py::test_stage8_is_marked_unimplemented_rather_than_quietly_missing` — passed

</details>

### Stage 9 (Tier B) — Candidate assessment and recommendation rules

**✅ PASSED** — 12 passed, 0 failed, 0 skipped, 0 not executed.

<details><summary>All test identifiers for this gate</summary>

* `tests/scientific/test_stage9_decisions.py::test_absent_batch_effect_yields_not_needed` — passed
* `tests/scientific/test_stage9_decisions.py::test_successful_correction_reaches_review_but_not_automatic_approval` — passed
* `tests/scientific/test_stage9_decisions.py::test_a_calibrated_policy_can_recommend` — passed
* `tests/scientific/test_stage9_decisions.py::test_information_destroying_correction_is_rejected` — passed
* `tests/scientific/test_stage9_decisions.py::test_sparse_validation_coverage_yields_insufficient_evidence` — passed
* `tests/scientific/test_stage9_decisions.py::test_confounded_design_is_non_identifiable_for_every_candidate` — passed
* `tests/scientific/test_stage9_decisions.py::test_below_floor_design_blocks_every_recommendation` — passed
* `tests/scientific/test_stage9_decisions.py::test_failed_and_unavailable_backends_get_distinct_states` — passed
* `tests/scientific/test_stage9_decisions.py::test_dropping_a_failed_metric_cannot_improve_a_recommendation` — passed
* `tests/scientific/test_stage9_decisions.py::test_clustering_labels_alone_cannot_support_a_recommendation` — passed
* `tests/scientific/test_stage9_decisions.py::test_policy_settings_and_threshold_origins_are_published` — passed
* `tests/scientific/test_stage9_decisions.py::test_overrides_are_recorded_beside_the_original_assessment` — passed

</details>

### Stage 10 (Tier A) — Reporting and safe exports

**🟡 PARTIAL** — 10 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** Minimal form, as Tier A specifies. Embeddings and independently fitted layout panels are not rendered, and the narrative summary is the rule-derived text in the template rather than a generated narrative section.

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage10_report.py::test_report_is_self_contained_and_offline` — passed
* `tests/integration/test_stage10_report.py::test_every_headline_number_is_traceable_to_a_record` — passed
* `tests/integration/test_stage10_report.py::test_failed_and_unavailable_candidates_stay_visible` — passed
* `tests/integration/test_stage10_report.py::test_missing_preservation_evidence_is_not_presented_as_preservation` — passed
* `tests/integration/test_stage10_report.py::test_identifiers_are_redacted_and_the_mapping_is_stored_separately` — passed
* `tests/integration/test_stage10_report.py::test_regeneration_triggers_no_correction` — passed
* `tests/integration/test_stage10_report.py::test_diagnosis_only_report_renders_before_any_correction` — passed
* `tests/integration/test_stage10_report.py::test_exports_round_trip_with_preserved_ids_and_declared_scale` — passed
* `tests/integration/test_stage10_report.py::test_status_file_reports_every_planned_unit_not_just_the_ones_that_ran` — passed
* `tests/integration/test_stage10_report.py::test_report_settings_do_not_change_the_analytical_fingerprint` — passed

</details>

### Stage 11 (Tier C) — Snakemake execution, performance, deployment

**🟡 PARTIAL** — 11 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** Deferred as Tier C, with the execution behaviour the Tier A slice depends on implemented and tested. NOT done: interruption-recovery testing, checkpoint-driven data-dependent branches (none are needed, because expansion is static), a built and digest-pinned Apptainer image, environment prefetching for offline nodes, and any cluster execution. The SLURM profile is shipped labelled UNTESTED.

<details><summary>All test identifiers for this gate</summary>

* `tests/integration/test_stage11_workflow.py::test_clean_execution_then_a_no_op_rerun` — passed
* `tests/integration/test_stage11_workflow.py::test_deleting_a_real_artifact_forces_recomputation` — passed
* `tests/integration/test_stage11_workflow.py::test_changing_the_configuration_invalidates_downstream_results` — passed
* `tests/integration/test_stage11_workflow.py::test_diagnose_mode_stops_before_any_correction` — passed
* `tests/integration/test_stage11_workflow.py::test_correct_mode_runs_only_the_selected_candidate` — passed
* `tests/integration/test_stage11_workflow.py::test_correct_mode_without_a_candidate_fails_at_parse_time` — passed
* `tests/integration/test_stage11_workflow.py::test_an_unavailable_method_produces_a_partial_report_not_a_failed_run` — passed
* `tests/integration/test_stage11_workflow.py::test_no_rule_hides_a_failure_behind_true_or_blanket_suppression` — passed
* `tests/integration/test_stage11_workflow.py::test_rules_declare_logs_threads_resources_and_benchmarks` — passed
* `tests/integration/test_stage11_workflow.py::test_the_candidate_roster_comes_from_the_configuration_not_the_results_tree` — passed
* `tests/integration/test_stage11_workflow.py::test_slurm_profile_is_present_but_labelled_untested` — passed

</details>

### Stage 12 (Tier C) — Scientific validation, documentation, release readiness

**⬜ DEFERRED** — 0 passed, 0 failed, 0 skipped, 0 not executed.

> **Partial:** Deferred as Tier C. Deterministic synthetic fixtures, a scientific suite covering technical-only shifts, biology-only differences, confounding, rare populations, constant markers, dispersion collapse, joint-structure corruption, within-batch drift and composition-only change DO exist and run. NOT done: a permissibly usable real cytometry example with provenance and checksums, multi-seed stochastic checks, R-side unit tests, a Snakemake lint step in CI, a compatibility matrix, and a release checklist. CITATION.cff and third-party notices are present.

## Report totals

* collected: 141
* passed: 139
* failed: 0
* skipped: 2
* errors: 0

Tests outside the gate lists (helpers and additional checks) are counted in the
totals but do not decide a gate.

