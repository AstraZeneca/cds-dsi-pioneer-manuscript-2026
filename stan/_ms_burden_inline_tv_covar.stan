// ============================================================================
// SHARED INLINE BURDEN TIME-VARYING COVARIATE (no states_full_grid)
// ============================================================================
// When ms_needs_inline_burden = TRUE (process noise OFF), compute the burden
// trajectory analytically per-patient and add the TV covariate contributions
// directly into log_cond_surv_{01,02,03}. This avoids materializing the dense
// states_full_grid matrix, which can be hundreds of thousands of autodiff vars.
//
// Required aliases (declared in each top-level model's transformed parameters):
//   vector[n_patients] log_baseline_burden
//   real median_log_burden_obs
//   real iqr_log_burden_obs
//
// MUST be included AFTER:  modules/multistate/transformed_parameters.stan
//   (needs log_cond_surv_{01,02,03} constructed but NOT yet transformed)
// MUST be included BEFORE: modules/multistate/cond_surv_transform.stan
//   (which applies -exp() to convert log-hazard to log conditional survival)

if (ms_needs_inline_burden) {
  // ── 0→2: Dense inline up to per-patient event time ────────────────────────
  if (enable_ms_02 && enable_ms_02_time_varying_cov &&
      enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
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
          real log_bburden = log_baseline_burden[p];
          real base_d = init_log_decrease_patient[j] + rate_d * first_visit;
          real base_g = init_log_growth_patient[j] - rate_g * first_visit;

          for (t in 1:T_max) {
            real s_d = base_d - rate_d * t;
            real s_g = base_g + rate_g * t;
            real std_burden = standardize_log_burden(
              log_sum_exp(s_d, s_g), log_bburden, median_log_burden_obs, iqr_log_burden_obs);
            log_cond_surv_02[j, t] += time_varying_coef_02[1] * std_burden;
          }

          if (enable_ms_velocity_basis) {
            // Velocity: central difference of capped analytic log-burden at each t.
            for (t in 1:T_max) {
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (t - 1), base_g + rate_g * (t - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (t + 1), base_g + rate_g * (t + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;  // central difference, unit grid
              log_cond_surv_02[j, t] += time_varying_coef_02[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          } else {
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
  }

  // ── 0→3: Dense inline up to per-patient dropout time ─────────────────────
  if (enable_ms_03 && enable_ms_03_time_varying_cov &&
      enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    profile("tv covariate inline 03") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, p);
        int first_visit = t_patient_visits[visit_start];
        int T_max = ms_time_03[p];

        if (T_max > 0) {
          real rate_d = patient_decrease_rate[j, 1];
          real rate_g = patient_growth_rate[j, 1];
          real log_bburden = log_baseline_burden[p];
          real base_d = init_log_decrease_patient[j] + rate_d * first_visit;
          real base_g = init_log_growth_patient[j] - rate_g * first_visit;

          for (t in 1:T_max) {
            real s_d = base_d - rate_d * t;
            real s_g = base_g + rate_g * t;
            real std_burden = standardize_log_burden(
              log_sum_exp(s_d, s_g), log_bburden, median_log_burden_obs, iqr_log_burden_obs);
            log_cond_surv_03[j, t] += time_varying_coef_03[1] * std_burden;
          }

          if (enable_ms_velocity_basis) {
            for (t in 1:T_max) {
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (t - 1), base_g + rate_g * (t - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (t + 1), base_g + rate_g * (t + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;
              log_cond_surv_03[j, t] += time_varying_coef_03[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          } else {
            if (n_time_varying_covar >= 2) {
              real contrib = time_varying_coef_03[2] * patient_log_decrease_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_03[j, t] += contrib;
            }
            if (n_time_varying_covar >= 3) {
              real contrib = time_varying_coef_03[3] * patient_log_growth_rate[j, 1];
              for (t in 1:T_max)
                log_cond_surv_03[j, t] += contrib;
            }
          }
        }
      }
    }
  }

  // ── 0→1: Sparse visit-gated latent burden trajectory ─────────────────────
  if (enable_ms_01 && enable_ms_pop_time_varying_cov &&
      enable_ms_visit_gated_01 && enable_ms_visit_gated_latent_01) {
    profile("tv covariate inline 01") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int v_start, v_end;
        (v_start, v_end) = get_pos(patient_visit_pos, p);
        int first_visit = t_patient_visits[v_start];
        real rate_d = patient_decrease_rate[j, 1];
        real rate_g = patient_growth_rate[j, 1];
        real log_bburden = log_baseline_burden[p];
        real base_d = init_log_decrease_patient[j] + rate_d * first_visit;
        real base_g = init_log_growth_patient[j] - rate_g * first_visit;

        for (v in v_start:v_end) {
          int wk = t_patient_visits[v];
          if (wk >= 1 && wk <= max_all_t) {
            real s_d = base_d - rate_d * wk;
            real s_g = base_g + rate_g * wk;
            real burden_covar = standardize_log_burden(
              log_sum_exp(s_d, s_g), log_bburden, median_log_burden_obs, iqr_log_burden_obs);
            log_cond_surv_01[j, wk] += time_varying_coef_01[1] * burden_covar;

            if (enable_ms_velocity_basis) {
              // Velocity at visit week wk: central difference of capped analytic
              // log-burden over a unit (weekly) grid straddling wk.
              real g_prev = fmin(log_sum_exp(base_d - rate_d * (wk - 1), base_g + rate_g * (wk - 1)), 10.0);
              real g_next = fmin(log_sum_exp(base_d - rate_d * (wk + 1), base_g + rate_g * (wk + 1)), 10.0);
              real vel = (g_next - g_prev) / 2.0;
              log_cond_surv_01[j, wk] += time_varying_coef_01[2]
                * standardize_velocity(vel, median_velocity_burden_obs, iqr_velocity_burden_obs);
            }
          }
        }
      }
    }
  }
}
