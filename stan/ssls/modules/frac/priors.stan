// frac/priors.stan — active priors for fraction (decrease share) module
// Multi-level hierarchy: loop over all levels for unified prior structure

// Population intercept
frac_logit_loc_pop ~ normal(frac_logit_loc_pop_mean, frac_logit_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_frac) {
  frac_coef_qr_pop ~ normal(frac_coef_qr_pop_mean, frac_coef_qr_pop_sd);
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
for (lv in 1:n_levels) {
  int lv_start = level_pos[lv];
  int lv_end = level_pos[lv + 1] - 1;

  // Intercept SD hyperprior
  frac_sd_level_intercept[lv] ~ normal(0, frac_sd_level_intercept_sd[lv]);

  // Intercept raw effects - ALWAYS apply prior to keep parameters bounded
  frac_raw_level_intercept[lv_start:lv_end] ~ std_normal();

  // Slope SD hyperprior - ALWAYS apply to keep parameter bounded
  if (n_covar > 0) {
    frac_sd_level_slope[lv] ~ normal(0, frac_sd_level_slope_sd[lv]);
  }
  // Slope raw effects - ALWAYS apply prior
  if (n_covar > 0) {
    to_vector(frac_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
  }
}
