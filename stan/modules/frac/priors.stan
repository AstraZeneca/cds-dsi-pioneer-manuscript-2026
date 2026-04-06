// frac/priors.stan — active priors for fraction (decrease share) module
// Multi-level hierarchy: loop over all levels for unified prior structure
// Uses enabled position arrays for efficient indexing into compacted parameter arrays

// Population intercept
frac_logit_loc_pop ~ normal(frac_logit_loc_pop_mean, frac_logit_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_frac) {
  frac_coef_qr_pop ~ normal(frac_coef_qr_pop_mean, frac_coef_qr_pop_sd);
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    // SD prior only for RE levels (mode=2); FE levels use fixed hyperparameter
    if (enable_level_intercept_frac[lv] == LEVEL_MODE_RE) {
      sd_idx += 1;
      frac_sd_level_intercept_raw[sd_idx] ~ normal(0, frac_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      frac_sd_level_slope[lv] ~ normal(0, frac_sd_level_slope_sd[lv]);
    }

    // Intercept raw effects - only apply prior to enabled levels (FE and RE both)
    // (parameter array is sized by enabled groups only)
    if (enable_level_intercept_frac[lv]) {
      int lv_start = enabled_level_pos_frac_intercept[lv];
      int lv_end = enabled_level_pos_frac_intercept[lv + 1] - 1;
      frac_raw_level_intercept[lv_start:lv_end] ~ std_normal();
    }

    // Slope raw effects - only apply prior to enabled levels
    if (enable_level_cov_frac[lv] && n_covar > 0) {
      int lv_start = enabled_level_pos_frac_slope[lv];
      int lv_end = enabled_level_pos_frac_slope[lv + 1] - 1;
      to_vector(frac_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}
