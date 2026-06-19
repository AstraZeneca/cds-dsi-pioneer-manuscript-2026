// gr_decay/priors.stan — flag-gated priors on log(kappa) with level hierarchy.
// Mirrors frac/priors.stan (minus Student-t). Entire block omitted when off.

if (enable_gr_decay) {
  gr_decay_log_loc_pop[1] ~ normal(gr_decay_log_loc_pop_mean, gr_decay_log_loc_pop_sd);
  if (enable_pop_cov_gr_decay) {
    gr_decay_coef_qr_pop ~ normal(gr_decay_coef_qr_pop_mean, gr_decay_coef_qr_pop_sd);
  }

  int sd_idx = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
        enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
      sd_idx += 1;
      gr_decay_sd_level_intercept_raw[sd_idx] ~ normal(0, gr_decay_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      gr_decay_sd_level_slope[lv] ~ normal(0, gr_decay_sd_level_slope_sd[lv]);
    }

    // Intercept effects
    if (enable_level_intercept_gr_decay[lv]) {
      if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_intercept[lv];
        int c_hi = cp_level_pos_gr_decay_intercept[lv + 1] - 1;
        gr_decay_cp_level_intercept[c_lo:c_hi] ~ normal(0, gr_decay_sd_level_intercept_raw[sd_idx]);
      } else {
        int r_lo = raw_level_pos_gr_decay_intercept[lv];
        int r_hi = raw_level_pos_gr_decay_intercept[lv + 1] - 1;
        gr_decay_raw_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // Slope effects
    if (enable_level_cov_gr_decay[lv] && n_covar > 0) {
      int mode = enable_level_intercept_gr_decay[lv];
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_gr_decay_slope[lv];
        int c_hi = cp_level_pos_gr_decay_slope[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_covar) {
            gr_decay_cp_level_slope[c_lo:c_hi, k] ~ normal(0, gr_decay_sd_level_slope[lv, k]);
          }
        }
      } else {
        int r_lo = raw_level_pos_gr_decay_slope[lv];
        int r_hi = raw_level_pos_gr_decay_slope[lv + 1] - 1;
        to_vector(gr_decay_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
      }
    }
  }
}
