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
    // SD prior for RE and RE_CP levels; FE levels use fixed hyperparameter
    if (enable_level_intercept_frac[lv] == LEVEL_MODE_RE ||
        enable_level_intercept_frac[lv] == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      frac_sd_level_intercept_raw[sd_idx] ~ normal(0, frac_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      frac_sd_level_slope[lv] ~ normal(0, frac_sd_level_slope_sd[lv]);
    }

    // Intercept effects — raw-path (FE/RE) or CP-path (RE_CP)
    if (enable_level_intercept_frac[lv]) {
      if (enable_level_intercept_frac[lv] == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_frac_intercept[lv];
        int c_hi = cp_level_pos_frac_intercept[lv + 1] - 1;
        if (enable_student_t_hierarchy)
          frac_cp_level_intercept[c_lo:c_hi] ~ student_t(frac_nu_level[lv], 0, frac_sd_level_intercept_raw[sd_idx]);
        else
          frac_cp_level_intercept[c_lo:c_hi] ~ normal(0, frac_sd_level_intercept_raw[sd_idx]);
      } else {
        int r_lo = raw_level_pos_frac_intercept[lv];
        int r_hi = raw_level_pos_frac_intercept[lv + 1] - 1;
        if (enable_student_t_hierarchy)
          frac_raw_level_intercept[r_lo:r_hi] ~ student_t(frac_nu_level[lv], 0, 1);
        else
          frac_raw_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // Slope effects — raw-path (FE/RE) or CP-path (RE_CP)
    if (enable_level_cov_frac[lv] && n_covar > 0) {
      int mode = enable_level_intercept_frac[lv];
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_frac_slope[lv];
        int c_hi = cp_level_pos_frac_slope[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_covar) {
            if (enable_student_t_hierarchy)
              frac_cp_level_slope[c_lo:c_hi, k] ~ student_t(frac_nu_level[lv], 0, frac_sd_level_slope[lv, k]);
            else
              frac_cp_level_slope[c_lo:c_hi, k] ~ normal(0, frac_sd_level_slope[lv, k]);
          }
        }
      } else {
        int r_lo = raw_level_pos_frac_slope[lv];
        int r_hi = raw_level_pos_frac_slope[lv + 1] - 1;
        if (enable_student_t_hierarchy)
          to_vector(frac_raw_level_slope[r_lo:r_hi, :]) ~ student_t(frac_nu_level[lv], 0, 1);
        else
          to_vector(frac_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
      }
    }
  }
}

// Student-t nu priors (only when enabled)
if (enable_student_t_hierarchy) {
  for (lv in 1:n_levels)
    frac_nu_level[lv] ~ gamma(frac_nu_level_prior_alpha[lv], frac_nu_level_prior_beta[lv]);
}
