// ============================================================================
// BUILD TIME-VARYING COVARIATE MATRIX FOR MULTISTATE (PSA MODEL)
// ============================================================================
// PSA-specific version of tumor/_ms_time_varying_covar.stan.
// Feature 1 uses log(absolute PSA) instead of log(SLD).
//
// MUST be included AFTER: modules/state_space/transformed_parameters.stan
// (needs states_full_grid, patient_log_decrease_rate, patient_log_growth_rate)
//
// MUST be included BEFORE: modules/multistate/transformed_parameters.stan
// (provides ms_time_varying_covar_01)

// Declare covariate matrix (sized by n_time_varying_covar from data)
// Structure: array[n_covar] of matrix[n_patients, max_all_t]
// Only allocated when both multistate and time-varying covariates are enabled
array[enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 ? n_time_varying_covar : 0]
  matrix[n_patients, max_all_t] ms_time_varying_covar_01;

if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Patient's first visit time in absolute weeks
    int first_visit = t_patient_visits[visit_start];

    // Mapping: absolute time t -> patient-relative column = t - first_visit + 1
    // Data invariant: first_visit <= 0, so states_start_col = 2 - first_visit >= 2
    int states_start_col = 2 - first_visit;  // Column for absolute time 1
    int states_end_col = states_start_col + max_all_t - 1;  // Column for absolute time max_all_t

    // Feature 1: Standardized log(PSA)
    // states_full_grid gives log(normalized_PSA) where normalized = ratio to baseline
    row_vector[max_all_t] log_psa_normalized = log_sum_exp(
      states_full_grid[1][i, states_start_col:states_end_col],
      states_full_grid[2][i, states_start_col:states_end_col]
    );
    // Add baseline to get absolute PSA value, then standardize
    row_vector[max_all_t] log_psa_absolute = log_baseline_psa[i] + log_psa_normalized;
    ms_time_varying_covar_01[1][i] = (log_psa_absolute - median_log_psa_obs) / iqr_log_psa_obs;

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
