// ============================================================================
// MS-Standalone LFO: Compact cutoff-cohort mapping for endpoint aggregation
// ============================================================================
// Built after _ms_standalone_lfo_transformed_data.stan so that
// `cutoff_last_visit_idx` is in scope.
//
// A patient is "cutoff-observed" iff they had at least one visit on or before
// the cutoff. For those patients only, we build:
//   - a compact index array `cutoff_observed_patients` (1..n_cutoff_observed_patients)
//   - a reverse lookup `patient_to_cutoff_idx`
//   - compact per-trial position boundaries `cutoff_trial_patient_pos`
//   - a compact cond_group membership (`cutoff_cond_group` + `cutoff_cond_group_pos`)
//
// These feed `calculate_all_patients_endpoints_rng` + aggregate_trial_metrics +
// aggregate_conditional_group_metrics in `_ms_standalone_lfo_os_km_generated_quantities.stan`,
// matching the pattern used by sclc's `_lfo_endpoints_generated_quantities.stan`.

// --- Compact cohort mask ----------------------------------------------------
array[n_patients] int cutoff_observed_mask;
int n_cutoff_observed_patients = 0;
for (i in 1:n_patients) {
  cutoff_observed_mask[i] = cutoff_last_visit_idx[i] > 0;
  if (cutoff_observed_mask[i]) n_cutoff_observed_patients += 1;
}
print("n_cutoff_observed_patients = ", n_cutoff_observed_patients);

// --- Compact index + reverse lookup ----------------------------------------
array[n_cutoff_observed_patients] int cutoff_observed_patients;
array[n_patients] int patient_to_cutoff_idx;
(cutoff_observed_patients, patient_to_cutoff_idx) =
  create_compact_patient_mapping(cutoff_observed_mask, n_cutoff_observed_patients);

// --- Compact per-trial position boundaries ---------------------------------
array[n_trials + 1] int cutoff_trial_patient_pos =
  create_compact_group_pos(cutoff_observed_patients, patient_trial, n_trials);

// --- Compact cond_group membership -----------------------------------------
// Filters cond_group membership to cutoff-observed patients and remaps the
// entries to compact indices (expected by aggregate_conditional_group_metrics).
array[n_cond_group + 1] int cutoff_cond_group_pos;
int n_cutoff_cond_group_entries = 0;
for (g_idx in 1:size(cond_group)) {
  if (cutoff_observed_mask[cond_group[g_idx]]) n_cutoff_cond_group_entries += 1;
}
array[n_cutoff_cond_group_entries] int cutoff_cond_group;
(cutoff_cond_group, cutoff_cond_group_pos) = remap_group_to_compact(
  cond_group, cond_group_pos, cutoff_observed_mask, patient_to_cutoff_idx, n_cond_group
);
