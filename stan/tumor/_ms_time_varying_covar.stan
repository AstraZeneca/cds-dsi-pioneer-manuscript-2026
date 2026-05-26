// ============================================================================
// BUILD TIME-VARYING COVARIATE MATRIX FOR MULTISTATE
// ============================================================================
// This file is MODEL-SPECIFIC: it knows about SLD, decrease_rate, growth_rate.
// The multistate module is GENERIC: it just multiplies coefficients by matrices.
//
// MUST be included AFTER: modules/state_space/transformed_parameters.stan
// (needs states_full_grid, patient_log_decrease_rate, patient_log_growth_rate)
//
// MUST be included BEFORE: modules/multistate/transformed_parameters.stan
// (provides ms_time_varying_covar_01)

// Declare covariate matrix (sized by n_time_varying_covar from data)
// Structure: array[n_covar] of matrix[n_patients, max_all_t]
// Only allocated when a downstream consumer of the dense grid actually exists:
//   - continuous-time 0->1 hazard (visit-gating off)
//   - visit-gated 0->1 in latent mode (reads dense grid at visit times)
//   - 0->2 time-varying covariate path (always uses dense grid)
// Visit-gated 0->1 in observed mode reads ms_obs_visit_covar_flat instead, so
// the dense grid is unnecessary then.
array[enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
      ((enable_ms_01 && !enable_ms_visit_gated_01) ||
       (enable_ms_01 && enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01) ||
       (enable_ms_02 && enable_ms_02_time_varying_cov))
      ? n_time_varying_covar : 0]
  matrix[n_patients, max_all_t] ms_time_varying_covar_01;

if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    ((enable_ms_01 && !enable_ms_visit_gated_01) ||
     (enable_ms_01 && enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01) ||
     (enable_ms_02 && enable_ms_02_time_varying_cov))) {
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Patient's first visit time in absolute weeks
    int first_visit = t_patient_visits[visit_start];

    // Mapping: absolute time t -> patient-relative column = t - first_visit + 1
    // Data invariant: first_visit <= 0, so states_start_col = 2 - first_visit >= 2
    int states_start_col = 2 - first_visit;  // Column for absolute time 1
    int states_end_col = states_start_col + max_all_t - 1;  // Column for absolute time max_all_t

    // Feature 1: Standardized log(SLD)
    // states_full_grid gives log(normalized_SLD) where normalized = ratio to baseline
    row_vector[max_all_t] log_sld_normalized = log_sum_exp(
      states_full_grid[1][i, states_start_col:states_end_col],
      states_full_grid[2][i, states_start_col:states_end_col]
    );
    // Add baseline to get absolute SLD value, then standardize
    row_vector[max_all_t] log_sld_absolute = log_baseline_sld[i] + log_sld_normalized;
    ms_time_varying_covar_01[1][i] = (log_sld_absolute - median_log_sld_obs) / iqr_log_sld_obs;

    // Feature 2: log(decrease rate) - already on log scale
    if (n_time_varying_covar >= 2) {
      if (enable_any_process_noise_tr) {
        // Time-varying rates: extract contiguous slice
        ms_time_varying_covar_01[2][i] = patient_log_decrease_rate[i, states_start_col:states_end_col];
      } else {
        // Constant rates: broadcast scalar
        ms_time_varying_covar_01[2][i] = rep_row_vector(patient_log_decrease_rate[i, 1], max_all_t);
      }
    }

    // Feature 3: log(growth rate) - already on log scale
    if (n_time_varying_covar >= 3) {
      if (enable_any_process_noise_tr) {
        ms_time_varying_covar_01[3][i] = patient_log_growth_rate[i, states_start_col:states_end_col];
      } else {
        ms_time_varying_covar_01[3][i] = rep_row_vector(patient_log_growth_rate[i, 1], max_all_t);
      }
    }
  }
}
