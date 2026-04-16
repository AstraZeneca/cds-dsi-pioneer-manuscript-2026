// ============================================================================
// INLINE PSA TIME-VARYING COVARIATE (no states_full_grid)
// ============================================================================
// PSA-specific optimization: when ms_needs_inline_psa = TRUE (process noise OFF),
// compute the PSA trajectory analytically per-patient and add the TV covariate
// directly into log_cond_surv_02 and log_cond_surv_01. This avoids materializing
// the dense states_full_grid matrix (783K autodiff vars for 1773 patients × 221 weeks).
//
// MUST be included AFTER: modules/multistate/transformed_parameters.stan
//   (needs log_cond_surv_01, log_cond_surv_02 constructed but NOT yet transformed)
// MUST be included BEFORE: modules/multistate/cond_surv_transform.stan
//   (which applies -exp() to convert log-hazard to log conditional survival)

if (ms_needs_inline_psa) {
  // ── 0→2: Dense inline computation up to per-patient event time ──────────
  if (enable_ms_02_time_varying_cov && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    profile("tv covariate inline 02") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, p);
        int first_visit = t_patient_visits[visit_start];
        int T_max = ms_max_time_02[p];

        if (T_max > 0) {
          real rate_d = patient_decrease_rate[j, 1];
          real rate_g = patient_growth_rate[j, 1];
          real log_bpsa = log_baseline_psa[p];

          // Precompute base offsets: state(t) = base + rate * t
          // Original: init - rate * (t - first_visit) = (init + rate * first_visit) - rate * t
          real base_d = init_log_decrease_patient[j] + rate_d * first_visit;
          real base_g = init_log_growth_patient[j] - rate_g * first_visit;

          // Feature 1: Standardized log PSA trajectory
          for (t in 1:T_max) {
            real s_d = base_d - rate_d * t;
            real s_g = base_g + rate_g * t;
            real std_psa = standardize_log_psa(
              log_sum_exp(s_d, s_g), log_bpsa, median_log_psa_obs, iqr_log_psa_obs);
            log_cond_surv_02[j, t] += time_varying_coef_02[1] * std_psa;
          }

          // Features 2-3: Constant rates (no process noise — single multiply, then broadcast)
          if (n_time_varying_covar >= 2) {
            real contrib = time_varying_coef_02[2] * patient_log_decrease_rate[j, 1];
            for (t in 1:T_max)
              log_cond_surv_02[j, t] += contrib;
          }
          if (n_time_varying_covar >= 3) {
            real contrib = time_varying_coef_02[3] * patient_log_growth_rate[j, 1];
            for (t in 1:T_max)
              log_cond_surv_02[j, t] += contrib;
          }
        }
      }
    }
  }

  // ── 0→1: Sparse visit-gated latent PSA ──────────────────────────────────
  if (enable_ms_pop_time_varying_cov && enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01) {
    profile("tv covariate inline 01") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int v_start; int v_end;
        (v_start, v_end) = get_pos(patient_visit_pos, p);
        int first_visit = t_patient_visits[v_start];
        real rate_d = patient_decrease_rate[j, 1];
        real rate_g = patient_growth_rate[j, 1];
        real log_bpsa = log_baseline_psa[p];
        real base_d = init_log_decrease_patient[j] + rate_d * first_visit;
        real base_g = init_log_growth_patient[j] - rate_g * first_visit;

        for (v in v_start:v_end) {
          int wk = t_patient_visits[v];
          if (wk >= 1 && wk <= max_all_t) {
            real s_d = base_d - rate_d * wk;
            real s_g = base_g + rate_g * wk;
            real psa_covar = standardize_log_psa(
              log_sum_exp(s_d, s_g), log_bpsa, median_log_psa_obs, iqr_log_psa_obs);
            log_cond_surv_01[j, wk] += time_varying_coef_01[1] * psa_covar;
          }
        }
      }
    }
  }
}
