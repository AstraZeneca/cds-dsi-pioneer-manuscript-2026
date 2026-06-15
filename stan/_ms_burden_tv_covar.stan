// ============================================================================
// SHARED MULTISTATE TIME-VARYING COVARIATE MATRIX BUILDER
// ============================================================================
// Biomarker-agnostic TV covariate matrix for the multistate hazard model.
// Used by both tumor (SLD) and PSA models — they share the same shrinkage-
// fraction state representation, so the only biomarker-specific information
// is the (log_baseline_burden, median_log_burden_obs, iqr_log_burden_obs)
// triple, which each model declares as aliases over its own (SLD/PSA) names
// before including this file.
//
// Required aliases (declared in each top-level model's transformed parameters):
//   vector[n_patients] log_baseline_burden
//   real median_log_burden_obs
//   real iqr_log_burden_obs
//
// MUST be included AFTER:  modules/state_space/transformed_parameters.stan
//   (needs states_full_grid, patient_log_decrease_rate, patient_log_growth_rate;
//    all of these are sized by n_forecast_patients)
// MUST be included BEFORE: modules/multistate/transformed_parameters.stan
//   (provides ms_time_varying_covar_01)
//
// Indexing convention: rows are forecast-local (j = 1..n_forecast_patients).
// Patient-data lookups go through forecast_patient_idx[j] to reach unified
// arrays such as t_patient_visits, patient_visit_pos, log_baseline_burden.

// Allocate covariate matrix only when the multistate machinery actually
// consumes it. When ms_needs_inline_burden, the burden trajectory is computed
// inline in _ms_burden_inline_tv_covar.stan instead, so this matrix is sized
// to 0. Otherwise it is needed when any of:
//   - 0->1 continuous mode (not visit-gated)
//   - 0->1 visit-gated with latent burden trajectory
//   - 0->2 with TV covariate
array[!ms_needs_inline_burden && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    ((enable_ms_01 && !enable_ms_visit_gated_01) || enable_ms_02_time_varying_cov ||
     (enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01))
    ? n_time_varying_covar : 0]
  matrix[n_forecast_patients, max_all_t] ms_time_varying_covar_01;

if (!ms_needs_inline_burden && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    ((enable_ms_01 && !enable_ms_visit_gated_01) || enable_ms_02_time_varying_cov ||
     (enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01))) {
  profile("tv covariate") {
    for (j in 1:n_forecast_patients) {
      int p = forecast_patient_idx[j];
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, p);

      int first_visit = t_patient_visits[visit_start];
      int states_start_col = 2 - first_visit;
      int states_end_col = states_start_col + max_all_t - 1;

      // Feature 1: standardized log(absolute burden).
      // states_full_grid carries log(burden / baseline). Cap before adding
      // baseline so the downstream log-hazard cap never binds during HMC
      // warmup (untreated growth excursions can push the linear predictor
      // into the saturating region and stall adaptation). Cap of 10 allows
      // burden up to exp(10)x baseline (~22000x) — physically impossible
      // over any clinical follow-up.
      row_vector[max_all_t] log_burden_normalized = fmin(log_sum_exp(
        states_full_grid[1][j, states_start_col:states_end_col],
        states_full_grid[2][j, states_start_col:states_end_col]
      ), 10.0);
      row_vector[max_all_t] log_burden_absolute = log_baseline_burden[p] + log_burden_normalized;
      ms_time_varying_covar_01[1][j] = (log_burden_absolute - median_log_burden_obs) / iqr_log_burden_obs;

      if (enable_ms_velocity_basis) {
        // Feature 2 (velocity mode): central-difference velocity of the
        // standardized log-burden trajectory. We difference the SAME capped
        // level values used in feature 1 (log_burden_normalized is already
        // fmin(·,10)-capped), then standardize with the velocity constants.
        // No additional cap on velocity (spec §3.1).
        row_vector[max_all_t] vel_raw = central_difference_row(log_burden_normalized);
        row_vector[max_all_t] vel_std;
        for (w in 1:max_all_t)
          vel_std[w] = standardize_velocity(vel_raw[w], median_velocity_burden_obs, iqr_velocity_burden_obs);
        ms_time_varying_covar_01[2][j] = vel_std;
      } else {
        // Feature 2 (legacy): log(decrease rate)
        if (n_time_varying_covar >= 2) {
          if (enable_any_process_noise_tr) {
            ms_time_varying_covar_01[2][j] = patient_log_decrease_rate[j, states_start_col:states_end_col];
          } else {
            ms_time_varying_covar_01[2][j] = rep_row_vector(patient_log_decrease_rate[j, 1], max_all_t);
          }
        }
        // Feature 3 (legacy): log(growth rate)
        if (n_time_varying_covar >= 3) {
          if (enable_any_process_noise_tr) {
            ms_time_varying_covar_01[3][j] = patient_log_growth_rate[j, states_start_col:states_end_col];
          } else {
            ms_time_varying_covar_01[3][j] = rep_row_vector(patient_log_growth_rate[j, 1], max_all_t);
          }
        }
      }
    }
  }
}
