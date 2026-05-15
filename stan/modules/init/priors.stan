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
    int mode = enable_level_intercept_init[lv];

    // SD prior for any level that samples its SD as a free parameter (RE or RE_CP)
    if (mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      init_sd_level_intercept_raw[sd_idx] ~ normal(0, init_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      init_sd_level_slope[lv] ~ normal(0, init_sd_level_slope_sd[lv]);
    }

    // RAW path: FE, RE, RE_GP all sample from std_normal / student_t(nu, 0, 1)
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_init_intercept[lv];
      int r_hi = raw_level_pos_init_intercept[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          init_raw_level_intercept[r_lo:r_hi] ~ student_t(init_nu_level[lv], 0, 1);
        else
          init_raw_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // CP path: RE_CP samples directly at scale sd
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_init_intercept[lv];
      int c_hi = cp_level_pos_init_intercept[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          init_cp_level_intercept[c_lo:c_hi]
            ~ student_t(init_nu_level[lv], 0, init_sd_level_intercept[lv]);
        else
          init_cp_level_intercept[c_lo:c_hi]
            ~ normal(0, init_sd_level_intercept[lv]);
      }
    }

    // Slope effects — parameterization follows the level's intercept mode
    if (enable_level_cov_init[lv] && n_covar > 0) {
      // RAW path (FE, RE, RE_GP)
      if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
        int r_lo = raw_level_pos_init_slope[lv];
        int r_hi = raw_level_pos_init_slope[lv + 1] - 1;
        if (r_hi >= r_lo) {
          if (enable_student_t_hierarchy)
            to_vector(init_raw_level_slope[r_lo:r_hi, :]) ~ student_t(init_nu_level[lv], 0, 1);
          else
            to_vector(init_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
        }
      }

      // CP path (RE_CP): sample at scale sd_level_slope[lv]
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_init_slope[lv];
        int c_hi = cp_level_pos_init_slope[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_covar) {
            if (enable_student_t_hierarchy)
              init_cp_level_slope[c_lo:c_hi, k]
                ~ student_t(init_nu_level[lv], 0, init_sd_level_slope[lv, k]);
            else
              init_cp_level_slope[c_lo:c_hi, k]
                ~ normal(0, init_sd_level_slope[lv, k]);
          }
        }
      }
    }
  }
}

// Student-t nu priors (only when enabled)
if (enable_student_t_hierarchy) {
  for (lv in 1:n_levels)
    init_nu_level[lv] ~ gamma(init_nu_level_prior_alpha[lv], init_nu_level_prior_beta[lv]);
}
