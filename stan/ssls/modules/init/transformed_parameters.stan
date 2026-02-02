// init/transformed_parameters.stan — active initial proportion module
// Multi-level hierarchy: loop over all levels to accumulate effects

// Population covariate effects
vector[n_patients] init_linpred_pop = enable_pop_cov_init ?
  (Q_covar_design_matrix * init_coef_qr_pop) : rep_vector(0, n_patients);

// ===== INTERCEPT EFFECTS =====
vector[n_patients] init_linpred_level_intercepts = rep_vector(0, n_patients);

for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    int lv_start = level_pos[lv];
    int lv_end = level_pos[lv + 1] - 1;

    vector[n_groups_per_level[lv]] level_effects =
      init_sd_level_intercept[lv] * init_raw_level_intercept[lv_start:lv_end];

    init_linpred_level_intercepts += level_effects[patient_level_groups[, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
vector[n_patients] init_linpred_level_slopes = rep_vector(0, n_patients);

for (lv in 1:n_levels) {
  if (enable_level_cov_init[lv] && n_covar > 0) {
    int lv_start = level_pos[lv];
    int lv_end = level_pos[lv + 1] - 1;

    matrix[n_groups_per_level[lv], n_covar] level_slopes_qr =
      init_raw_level_slope[lv_start:lv_end, :] .*
      rep_matrix(init_sd_level_slope[lv]', n_groups_per_level[lv]);

    init_linpred_level_slopes += rows_dot_product(
      Q_covar_design_matrix,
      level_slopes_qr[patient_level_groups[, lv], :]
    );
  }
}

// ===== FINAL LINEAR PREDICTOR =====
vector[n_patients] init_logit_loc_patient = init_logit_loc_pop
  + init_linpred_pop
  + init_linpred_level_intercepts
  + init_linpred_level_slopes;

vector[n_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);
vector[n_patients] init_log_growth_patient   = log1m_inv_logit(init_logit_loc_patient);
