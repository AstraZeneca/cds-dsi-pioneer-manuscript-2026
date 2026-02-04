// init/transformed_parameters.stan — active initial proportion module
// Optimized: uses pre-computed flat indices from transformed_data for direct gather

// Population covariate effects
vector[n_patients] init_linpred_pop = enable_pop_cov_init ?
  (Q_covar_design_matrix * init_coef_qr_pop) : zeros_vector(n_patients);

// ===== INTERCEPT EFFECTS =====
// Step 1: Scale all raw effects at once
vector[n_enabled_groups_init_intercept] init_scaled_level_intercept;
for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(enabled_level_pos_init_intercept, lv);
    init_scaled_level_intercept[lv_start:lv_end] =
      init_sd_level_intercept[lv] * init_raw_level_intercept[lv_start:lv_end];
  }
}

// Step 2: Gather using pre-computed flat indices
vector[n_patients] init_linpred_level_intercepts = zeros_vector(n_patients);
for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    init_linpred_level_intercepts += init_scaled_level_intercept[patient_init_intercept_flat_idx[, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
// Step 1: Scale all raw slope effects at once
matrix[n_enabled_groups_init_slope, n_covar] init_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_init_slope > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_init[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_init_slope, lv);
      init_scaled_level_slope[lv_start:lv_end, :] =
        init_raw_level_slope[lv_start:lv_end, :] .*
        rep_matrix(init_sd_level_slope[lv]', lv_end - lv_start + 1);
    }
  }
}

// Step 2: Gather and compute dot products
vector[n_patients] init_linpred_level_slopes = zeros_vector(n_patients);
if (n_covar > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_init[lv]) {
      init_linpred_level_slopes += rows_dot_product(
        Q_covar_design_matrix,
        init_scaled_level_slope[patient_init_slope_flat_idx[, lv], :]
      );
    }
  }
}

// ===== FINAL LINEAR PREDICTOR =====
vector[n_patients] init_logit_loc_patient = init_logit_loc_pop
  + init_linpred_pop
  + init_linpred_level_intercepts
  + init_linpred_level_slopes;

vector[n_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);
vector[n_patients] init_log_growth_patient   = log1m_inv_logit(init_logit_loc_patient);
