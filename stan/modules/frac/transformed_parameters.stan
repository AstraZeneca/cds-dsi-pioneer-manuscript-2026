// frac/transformed_parameters.stan — active fraction module transformed parameters
// Optimized: uses pre-computed flat indices from transformed_data for direct gather

// Population covariate effects
vector[n_forecast_patients] frac_linpred_pop = enable_pop_cov_frac ?
  (Q_covar_design_matrix[forecast_patient_idx, :] * frac_coef_qr_pop) : zeros_vector(n_forecast_patients);

// ===== INTERCEPT EFFECTS =====
// Step 1: Scale all raw effects at once
vector[n_enabled_groups_frac_intercept] frac_scaled_level_intercept;
for (lv in 1:n_levels) {
  if (enable_level_intercept_frac[lv]) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(enabled_level_pos_frac_intercept, lv);
    frac_scaled_level_intercept[lv_start:lv_end] =
      frac_sd_level_intercept[lv] * frac_raw_level_intercept[lv_start:lv_end];
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
// Step 1: Scale all raw slope effects at once
matrix[n_enabled_groups_frac_slope, n_covar] frac_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_frac_slope > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_frac[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_frac_slope, lv);
      frac_scaled_level_slope[lv_start:lv_end, :] =
        frac_raw_level_slope[lv_start:lv_end, :] .*
        rep_matrix(frac_sd_level_slope[lv]', lv_end - lv_start + 1);
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
