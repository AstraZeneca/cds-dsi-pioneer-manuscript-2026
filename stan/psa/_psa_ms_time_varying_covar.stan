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
// Structure: array[n_covar] of matrix[n_hmc_patients, max_all_t]
// Only allocated when both multistate and time-varying covariates are enabled.
// Rows correspond to HMC patients (j = 1..n_hmc_patients); data arrays are
// accessed via hmc_patient_idx[j] to map to the unified patient arrays.
// Build modeled-PSA covariate matrix only when needed:
// - 0->1 continuous mode (not visit-gated), OR
// - 0->2 with time-varying covariate enabled
array[enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    (enable_ms_01 && !enable_ms_visit_gated_01 || enable_ms_02_time_varying_cov) ? n_time_varying_covar : 0]
  matrix[n_hmc_patients, max_all_t] ms_time_varying_covar_01;

if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    (enable_ms_01 && !enable_ms_visit_gated_01 || enable_ms_02_time_varying_cov)) {
  for (j in 1:n_hmc_patients) {
    int p = hmc_patient_idx[j];  // Unified patient index
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, p);

    // Patient's first visit time in absolute weeks
    int first_visit = t_patient_visits[visit_start];

    // Mapping: absolute time t -> patient-relative column = t - first_visit + 1
    // Data invariant: first_visit <= 0, so states_start_col = 2 - first_visit >= 2
    int states_start_col = 2 - first_visit;  // Column for absolute time 1
    int states_end_col = states_start_col + max_all_t - 1;  // Column for absolute time max_all_t

    // Feature 1: Standardized log(PSA)
    // states_full_grid gives log(normalized_PSA) where normalized = ratio to baseline.
    // states_full_grid is indexed by j (HMC-local), not p (unified).
    row_vector[max_all_t] log_psa_normalized = log_sum_exp(
      states_full_grid[1][j, states_start_col:states_end_col],
      states_full_grid[2][j, states_start_col:states_end_col]
    );
    // Cap log_psa_normalized so the log-hazard cap in transformed_parameters.stan
    // never binds during HMC warmup. With cap=10 here:
    //   log_psa_absolute ≤ log_bpsa + 10 ≈ 14
    //   standardized     ≤ (14 - 2.5) / 1 ≈ 11.5
    //   log_hazard       ≤ -4.5 + beta_max × 11.5 ≈ 5  (well below the cap=10)
    // When log-hazard hits its ceiling, fmin makes the gradient zero → HMC stuck.
    // This cap prevents that region from being reached; the hazard cap is then a
    // backstop only. Cap of 10 allows PSA up to exp(10)× baseline (22000×), which
    // is physically impossible over any clinical follow-up period.
    log_psa_normalized = fmin(log_psa_normalized, 10.0);
    // Add baseline to get absolute PSA value, then standardize.
    // log_baseline_psa is a unified array indexed by p.
    row_vector[max_all_t] log_psa_absolute = log_baseline_psa[p] + log_psa_normalized;
    ms_time_varying_covar_01[1][j] = (log_psa_absolute - median_log_psa_obs) / iqr_log_psa_obs;

    // Feature 2: log(decrease rate) - already on log scale.
    // patient_log_decrease_rate is HMC-local (indexed by j).
    if (n_time_varying_covar >= 2) {
      if (enable_any_process_noise_tr) {
        // Time-varying rates: extract contiguous slice
        ms_time_varying_covar_01[2][j] = patient_log_decrease_rate[j, states_start_col:states_end_col];
      } else {
        // Constant rates: broadcast scalar
        ms_time_varying_covar_01[2][j] = rep_row_vector(patient_log_decrease_rate[j, 1], max_all_t);
      }
    }

    // Feature 3: log(growth rate) - already on log scale.
    // patient_log_growth_rate is HMC-local (indexed by j).
    if (n_time_varying_covar >= 3) {
      if (enable_any_process_noise_tr) {
        ms_time_varying_covar_01[3][j] = patient_log_growth_rate[j, states_start_col:states_end_col];
      } else {
        ms_time_varying_covar_01[3][j] = rep_row_vector(patient_log_growth_rate[j, 1], max_all_t);
      }
    }
  }
}
