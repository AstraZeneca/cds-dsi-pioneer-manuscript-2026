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
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    // SD prior only for RE levels (mode=2); FE levels use fixed hyperparameter
    if (enable_level_intercept_init[lv] == LEVEL_MODE_RE) {
      sd_idx += 1;
      init_sd_level_intercept_raw[sd_idx] ~ normal(0, init_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      init_sd_level_slope[lv] ~ normal(0, init_sd_level_slope_sd[lv]);
    }

    // Intercept raw effects - only apply prior to enabled levels (FE and RE both)
    // (parameter array is sized by enabled groups only)
    if (enable_level_intercept_init[lv]) {
      int lv_start = enabled_level_pos_init_intercept[lv];
      int lv_end = enabled_level_pos_init_intercept[lv + 1] - 1;
      init_raw_level_intercept[lv_start:lv_end] ~ std_normal();
    }

    // Slope raw effects - only apply prior to enabled levels
    if (enable_level_cov_init[lv] && n_covar > 0) {
      int lv_start = enabled_level_pos_init_slope[lv];
      int lv_end = enabled_level_pos_init_slope[lv + 1] - 1;
      to_vector(init_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}
