// init/priors.stan — active priors for initial proportion module
// Multi-level hierarchy: loop over all levels for unified prior structure
// Uses enabled position arrays for efficient indexing into compacted parameter arrays

// Population intercept
init_logit_loc_pop ~ normal(init_logit_loc_pop_mean, init_logit_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_init) {
  init_coef_qr_pop ~ normal(init_coef_qr_pop_mean, init_coef_qr_pop_sd);
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
for (lv in 1:n_levels) {
  // SD hyperpriors are always applied (arrays are always n_levels)
  init_sd_level_intercept[lv] ~ normal(0, init_sd_level_intercept_sd[lv]);
  if (n_covar > 0) {
    init_sd_level_slope[lv] ~ normal(0, init_sd_level_slope_sd[lv]);
  }

  // Intercept raw effects - only apply prior to enabled levels
  // (parameter array is sized by enabled groups only)
  if (enable_level_intercept_init[lv]) {
    int lv_start = enabled_level_pos_init_intercept[lv];
    int lv_end = enabled_level_pos_init_intercept[lv + 1] - 1;
    if (enable_student_t_hierarchy)
      init_raw_level_intercept[lv_start:lv_end] ~ student_t(init_nu_level[lv], 0, 1);
    else
      init_raw_level_intercept[lv_start:lv_end] ~ std_normal();
  }

  // Slope raw effects - only apply prior to enabled levels
  if (enable_level_cov_init[lv] && n_covar > 0) {
    int lv_start = enabled_level_pos_init_slope[lv];
    int lv_end = enabled_level_pos_init_slope[lv + 1] - 1;
    if (enable_student_t_hierarchy)
      to_vector(init_raw_level_slope[lv_start:lv_end, :]) ~ student_t(init_nu_level[lv], 0, 1);
    else
      to_vector(init_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
  }
}

// Student-t nu priors (only when enabled)
if (enable_student_t_hierarchy) {
  for (lv in 1:n_levels)
    init_nu_level[lv] ~ gamma(init_nu_level_prior_alpha[lv], init_nu_level_prior_beta[lv]);
}
