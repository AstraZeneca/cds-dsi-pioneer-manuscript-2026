// gr_decay/transformed_parameters.stan
// Per-patient decay rate kappa_i = exp(log-scale linear predictor):
//   gr_decay_log_loc_pop + pop-cov + Σ_lv level intercepts + Σ_lv level slopes
// Mirrors frac/transformed_parameters.stan. Indexed by forecast-local patient j.
// When enable_gr_decay is off, kappa is zeros (never consumed — warp takes t-branch).

vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
vector[n_forecast_patients] gr_decay_kappa = zeros_vector(n_forecast_patients);

if (enable_gr_decay) {
  // Population covariate effects
  vector[n_forecast_patients] gr_decay_linpred_pop = enable_pop_cov_gr_decay ?
    (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop) : zeros_vector(n_forecast_patients);

  // ===== SD EXPANSION (frac-style) =====
  array[n_levels] real gr_decay_sd_level_intercept;
  {
    int sd_idx = 0;
    for (lv in 1:n_levels) {
      if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_FE) {
        gr_decay_sd_level_intercept[lv] = gr_decay_fe_sd_level_intercept[lv];
      } else if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
                 enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
        sd_idx += 1;
        gr_decay_sd_level_intercept[lv] = gr_decay_sd_level_intercept_raw[sd_idx];
      } else {
        gr_decay_sd_level_intercept[lv] = 0.0;
      }
    }
  }

  // ===== INTERCEPT EFFECTS =====
  vector[n_enabled_groups_gr_decay_intercept] gr_decay_scaled_level_intercept;
  for (lv in 1:n_levels) {
    int mode = enable_level_intercept_gr_decay[lv];
    if (mode == LEVEL_MODE_NONE) continue;
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_gr_decay_intercept, lv);
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_gr_decay_intercept[lv];
      int c_hi = cp_level_pos_gr_decay_intercept[lv + 1] - 1;
      gr_decay_scaled_level_intercept[e_lo:e_hi] = gr_decay_cp_level_intercept[c_lo:c_hi];
    } else {
      int r_lo = raw_level_pos_gr_decay_intercept[lv];
      int r_hi = raw_level_pos_gr_decay_intercept[lv + 1] - 1;
      gr_decay_scaled_level_intercept[e_lo:e_hi] =
        gr_decay_sd_level_intercept[lv] * gr_decay_raw_level_intercept[r_lo:r_hi];
    }
  }

  vector[n_forecast_patients] gr_decay_linpred_level_intercepts = zeros_vector(n_forecast_patients);
  for (lv in 1:n_levels) {
    if (enable_level_intercept_gr_decay[lv]) {
      gr_decay_linpred_level_intercepts += gr_decay_scaled_level_intercept[patient_gr_decay_intercept_flat_idx[forecast_patient_idx, lv]];
    }
  }

  // ===== COVARIATE SLOPE EFFECTS =====
  matrix[n_enabled_groups_gr_decay_slope, n_covar] gr_decay_scaled_level_slope;
  if (n_covar > 0 && n_enabled_groups_gr_decay_slope > 0) {
    for (lv in 1:n_levels) {
      if (!enable_level_cov_gr_decay[lv]) continue;
      int mode = enable_level_intercept_gr_decay[lv];
      int e_lo, e_hi;
      (e_lo, e_hi) = get_pos(enabled_level_pos_gr_decay_slope, lv);
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_slope[lv];
        int c_hi = cp_level_pos_gr_decay_slope[lv + 1] - 1;
        gr_decay_scaled_level_slope[e_lo:e_hi, :] = gr_decay_cp_level_slope[c_lo:c_hi, :];
      } else {
        int r_lo = raw_level_pos_gr_decay_slope[lv];
        int r_hi = raw_level_pos_gr_decay_slope[lv + 1] - 1;
        gr_decay_scaled_level_slope[e_lo:e_hi, :] =
          gr_decay_raw_level_slope[r_lo:r_hi, :] .*
          rep_matrix(gr_decay_sd_level_slope[lv]', r_hi - r_lo + 1);
      }
    }
  }

  vector[n_forecast_patients] gr_decay_linpred_level_slopes = zeros_vector(n_forecast_patients);
  if (n_covar > 0) {
    for (lv in 1:n_levels) {
      if (enable_level_cov_gr_decay[lv]) {
        gr_decay_linpred_level_slopes += rows_dot_product(
          Q_covar_design_matrix[forecast_patient_idx, :],
          gr_decay_scaled_level_slope[patient_gr_decay_slope_flat_idx[forecast_patient_idx, lv], :]
        );
      }
    }
  }

  // ===== FINAL LINEAR PREDICTOR =====
  gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients)
    + gr_decay_linpred_pop
    + gr_decay_linpred_level_intercepts
    + gr_decay_linpred_level_slopes;
  gr_decay_kappa = exp(gr_decay_log_loc_patient);
}
