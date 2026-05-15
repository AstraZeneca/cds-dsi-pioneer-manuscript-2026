// frac/transformed_parameters.stan — active fraction module transformed parameters
// Optimized: uses pre-computed flat indices from transformed_data for direct gather

// Population covariate effects
vector[n_forecast_patients] frac_linpred_pop = enable_pop_cov_frac ?
  (Q_covar_design_matrix[forecast_patient_idx, :] * frac_coef_qr_pop) : zeros_vector(n_forecast_patients);

// ===== SD EXPANSION =====
// Assemble full n_levels SD array from RE/RE_CP free params and FE hyperparams.
// FE levels (mode=1): use fixed frac_fe_sd_level_intercept[lv] from data.
// RE/RE_CP levels (mode=2/4): use free parameter frac_sd_level_intercept_raw (sequential counter).
// Disabled levels (mode=0): set to 0.0 (never used in intercept scaling).
array[n_levels] real<lower=0> frac_sd_level_intercept;
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_frac[lv] == LEVEL_MODE_FE) {
      frac_sd_level_intercept[lv] = frac_fe_sd_level_intercept[lv];
    } else if (enable_level_intercept_frac[lv] == LEVEL_MODE_RE ||
               enable_level_intercept_frac[lv] == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      frac_sd_level_intercept[lv] = frac_sd_level_intercept_raw[sd_idx];
    } else {
      frac_sd_level_intercept[lv] = 0.0;
    }
  }
}

// ===== INTERCEPT EFFECTS =====
// Step 1: Route each level to raw-NCP or CP bucket and scale appropriately
vector[n_enabled_groups_frac_intercept] frac_scaled_level_intercept;
for (lv in 1:n_levels) {
  int mode = enable_level_intercept_frac[lv];
  if (mode == LEVEL_MODE_NONE) continue;
  int e_lo, e_hi;
  (e_lo, e_hi) = get_pos(enabled_level_pos_frac_intercept, lv);

  if (mode == LEVEL_MODE_RE_CP) {
    int c_lo = cp_level_pos_frac_intercept[lv];
    int c_hi = cp_level_pos_frac_intercept[lv + 1] - 1;
    frac_scaled_level_intercept[e_lo:e_hi] = frac_cp_level_intercept[c_lo:c_hi];
  } else {
    int r_lo = raw_level_pos_frac_intercept[lv];
    int r_hi = raw_level_pos_frac_intercept[lv + 1] - 1;
    frac_scaled_level_intercept[e_lo:e_hi] =
      frac_sd_level_intercept[lv] * frac_raw_level_intercept[r_lo:r_hi];
  }
}

// Step 2: Gather using pre-computed flat indices
vector[n_forecast_patients] frac_linpred_level_intercepts = zeros_vector(n_forecast_patients);
for (lv in 1:n_levels) {
  if (enable_level_intercept_frac[lv]) {
    frac_linpred_level_intercepts += frac_scaled_level_intercept[patient_frac_intercept_flat_idx[forecast_patient_idx, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
// Step 1: Route each level to raw-NCP or CP bucket and scale appropriately
matrix[n_enabled_groups_frac_slope, n_covar] frac_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_frac_slope > 0) {
  for (lv in 1:n_levels) {
    if (!enable_level_cov_frac[lv]) continue;
    int mode = enable_level_intercept_frac[lv];
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_frac_slope, lv);

    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_frac_slope[lv];
      int c_hi = cp_level_pos_frac_slope[lv + 1] - 1;
      frac_scaled_level_slope[e_lo:e_hi, :] = frac_cp_level_slope[c_lo:c_hi, :];
    } else {
      int r_lo = raw_level_pos_frac_slope[lv];
      int r_hi = raw_level_pos_frac_slope[lv + 1] - 1;
      frac_scaled_level_slope[e_lo:e_hi, :] =
        frac_raw_level_slope[r_lo:r_hi, :] .*
        rep_matrix(frac_sd_level_slope[lv]', r_hi - r_lo + 1);
    }
  }
}

// Step 2: Gather and compute dot products
vector[n_forecast_patients] frac_linpred_level_slopes = zeros_vector(n_forecast_patients);
if (n_covar > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_frac[lv]) {
      frac_linpred_level_slopes += rows_dot_product(
        Q_covar_design_matrix[forecast_patient_idx, :],
        frac_scaled_level_slope[patient_frac_slope_flat_idx[forecast_patient_idx, lv], :]
      );
    }
  }
}

// ===== FINAL LINEAR PREDICTOR =====
vector[n_forecast_patients] frac_logit_loc_patient = frac_logit_loc_pop
  + frac_linpred_pop
  + frac_linpred_level_intercepts
  + frac_linpred_level_slopes;

vector[n_forecast_patients] frac_log_decrease_patient = log_inv_logit(frac_logit_loc_patient);
vector[n_forecast_patients] frac_log_growth_patient   = log1m_inv_logit(frac_logit_loc_patient);
