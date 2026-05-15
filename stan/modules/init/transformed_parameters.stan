// init/transformed_parameters.stan — active initial proportion module
// Optimized: uses pre-computed flat indices from transformed_data for direct gather

// Population covariate effects
vector[n_forecast_patients] init_linpred_pop = enable_pop_cov_init ?
  (Q_covar_design_matrix[forecast_patient_idx, :] * init_coef_qr_pop) : zeros_vector(n_forecast_patients);

// ===== SD EXPANSION =====
// Assemble full n_levels SD array from RE free params and FE hyperparams.
// FE levels (mode=1): use fixed init_fe_sd_level_intercept[lv] from data.
// RE levels (mode=2): use free parameter init_sd_level_intercept_raw (sequential counter).
// Disabled levels (mode=0): set to 0.0 (never used in intercept scaling).
array[n_levels] real<lower=0> init_sd_level_intercept;
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_init[lv] == LEVEL_MODE_FE) {
      init_sd_level_intercept[lv] = init_fe_sd_level_intercept[lv];
    } else if (enable_level_intercept_init[lv] == LEVEL_MODE_RE ||
               enable_level_intercept_init[lv] == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      init_sd_level_intercept[lv] = init_sd_level_intercept_raw[sd_idx];
    } else {
      init_sd_level_intercept[lv] = 0.0;
    }
  }
}

// ===== INTERCEPT EFFECTS =====
// Step 1: Scale all raw effects at once
vector[n_enabled_groups_init_intercept] init_scaled_level_intercept;
for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    int mode = enable_level_intercept_init[lv];
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_init_intercept, lv);
    if (mode == LEVEL_MODE_RE_CP) {
      // CP: scaled = centered (identity — the _cp_ vector is already at natural scale)
      int c_lo = cp_level_pos_init_intercept[lv];
      int c_hi = cp_level_pos_init_intercept[lv + 1] - 1;
      init_scaled_level_intercept[e_lo:e_hi] = init_cp_level_intercept[c_lo:c_hi];
    } else {
      // NCP (FE, RE, RE_GP): scaled = sd[lv] * raw
      int r_lo = raw_level_pos_init_intercept[lv];
      int r_hi = raw_level_pos_init_intercept[lv + 1] - 1;
      init_scaled_level_intercept[e_lo:e_hi] =
        init_sd_level_intercept[lv] * init_raw_level_intercept[r_lo:r_hi];
    }
  }
}

// Step 2: Gather using pre-computed flat indices
vector[n_forecast_patients] init_linpred_level_intercepts = zeros_vector(n_forecast_patients);
for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    init_linpred_level_intercepts += init_scaled_level_intercept[patient_init_intercept_flat_idx[forecast_patient_idx, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
// Step 1: Scale all raw slope effects at once
matrix[n_enabled_groups_init_slope, n_covar] init_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_init_slope > 0) {
  for (lv in 1:n_levels) {
    if (!enable_level_cov_init[lv]) continue;
    int mode = enable_level_intercept_init[lv];
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_init_slope, lv);
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_init_slope[lv];
      int c_hi = cp_level_pos_init_slope[lv + 1] - 1;
      init_scaled_level_slope[e_lo:e_hi, :] = init_cp_level_slope[c_lo:c_hi, :];
    } else {
      int r_lo = raw_level_pos_init_slope[lv];
      int r_hi = raw_level_pos_init_slope[lv + 1] - 1;
      init_scaled_level_slope[e_lo:e_hi, :] =
        init_raw_level_slope[r_lo:r_hi, :] .*
        rep_matrix(init_sd_level_slope[lv]', r_hi - r_lo + 1);
    }
  }
}

// Step 2: Gather and compute dot products
vector[n_forecast_patients] init_linpred_level_slopes = zeros_vector(n_forecast_patients);
if (n_covar > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_init[lv]) {
      init_linpred_level_slopes += rows_dot_product(
        Q_covar_design_matrix[forecast_patient_idx, :],
        init_scaled_level_slope[patient_init_slope_flat_idx[forecast_patient_idx, lv], :]
      );
    }
  }
}

// ===== FINAL LINEAR PREDICTOR =====
vector[n_forecast_patients] init_logit_loc_patient = init_logit_loc_pop
  + init_linpred_pop
  + init_linpred_level_intercepts
  + init_linpred_level_slopes;

vector[n_forecast_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);
vector[n_forecast_patients] init_log_growth_patient   = log1m_inv_logit(init_logit_loc_patient);
