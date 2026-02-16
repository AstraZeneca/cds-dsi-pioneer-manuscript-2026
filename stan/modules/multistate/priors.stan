// ============================================================================
// Multistate Hazard Model Priors
// ============================================================================

// ============================================================================
// 0→1 TRANSITION PRIORS
// ============================================================================

if (enable_ms_01) {
  // Population-level baseline hazard GP
  ms_log_lambda_gp_01_pop_alpha[1] ~ inv_gamma(
    ms_log_lambda_gp_01_pop_alpha_alpha, ms_log_lambda_gp_01_pop_alpha_beta
  );
  ms_log_lambda_gp_01_pop_rho[1] ~ inv_gamma(
    ms_log_lambda_gp_01_pop_rho_alpha, ms_log_lambda_gp_01_pop_rho_beta
  );
  ms_log_lambda_gp_01_pop_intercept[1] ~ normal(
    ms_log_lambda_gp_01_pop_intercept_mean, ms_log_lambda_gp_01_pop_intercept_sd
  );
  to_vector(ms_log_lambda_gp_01_pop_eta) ~ std_normal();

  // Level-level baseline hazard GP
  // Always set priors for all levels (even disabled ones) to avoid improper posteriors
  for (lv in 1:n_levels) {
    ms_log_lambda_gp_01_level_alpha[lv] ~ inv_gamma(
      ms_log_lambda_gp_01_level_alpha_alpha[lv], ms_log_lambda_gp_01_level_alpha_beta[lv]
    );
    ms_log_lambda_gp_01_level_rho[lv] ~ inv_gamma(
      ms_log_lambda_gp_01_level_rho_alpha[lv], ms_log_lambda_gp_01_level_rho_beta[lv]
    );
    ms_log_lambda_gp_01_level_intercept_sd[lv] ~ normal(
      0, ms_log_lambda_gp_01_level_intercept_sd_sd[lv]
    );
  }

  // Group-level effects (only for enabled levels)
  if (n_enabled_groups_ms_baseline_01 > 0) {
    ms_raw_log_lambda_gp_01_level_intercept ~ std_normal();
    to_vector(ms_log_lambda_gp_01_level_eta) ~ std_normal();
  }

  // Time-varying covariate coefficients
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    ms_time_varying_coef_01 ~ normal(ms_time_varying_coef_01_mean, ms_time_varying_coef_01_sd);
  }

  // Time-invariant covariate coefficients
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    ms_time_invariant_coef_qr_01 ~ normal(ms_time_invariant_coef_01_mean, ms_time_invariant_coef_01_sd);
  }

  // Multi-level random slopes
  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      ms_sd_level_slope_01[lv] ~ normal(0, ms_sd_level_slope_01_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int lv_start = enabled_level_pos_ms_slope[lv];
      int lv_end = enabled_level_pos_ms_slope[lv + 1] - 1;
      to_vector(ms_raw_level_slope_01[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}

// ============================================================================
// 0→2 TRANSITION PRIORS
// ============================================================================

if (enable_ms_02) {
  // Population-level baseline hazard GP
  ms_log_lambda_gp_02_pop_alpha[1] ~ inv_gamma(
    ms_log_lambda_gp_02_pop_alpha_alpha, ms_log_lambda_gp_02_pop_alpha_beta
  );
  ms_log_lambda_gp_02_pop_rho[1] ~ inv_gamma(
    ms_log_lambda_gp_02_pop_rho_alpha, ms_log_lambda_gp_02_pop_rho_beta
  );
  ms_log_lambda_gp_02_pop_intercept[1] ~ normal(
    ms_log_lambda_gp_02_pop_intercept_mean, ms_log_lambda_gp_02_pop_intercept_sd
  );
  to_vector(ms_log_lambda_gp_02_pop_eta) ~ std_normal();

  // Level-level baseline hazard GP
  for (lv in 1:n_levels) {
    ms_log_lambda_gp_02_level_alpha[lv] ~ inv_gamma(
      ms_log_lambda_gp_02_level_alpha_alpha[lv], ms_log_lambda_gp_02_level_alpha_beta[lv]
    );
    ms_log_lambda_gp_02_level_rho[lv] ~ inv_gamma(
      ms_log_lambda_gp_02_level_rho_alpha[lv], ms_log_lambda_gp_02_level_rho_beta[lv]
    );
    ms_log_lambda_gp_02_level_intercept_sd[lv] ~ normal(
      0, ms_log_lambda_gp_02_level_intercept_sd_sd[lv]
    );
  }
  if (n_enabled_groups_ms_baseline_02 > 0) {
    ms_raw_log_lambda_gp_02_level_intercept ~ std_normal();
    to_vector(ms_log_lambda_gp_02_level_eta) ~ std_normal();
  }

  // Time-varying covariate coefficients
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    ms_time_varying_coef_02 ~ normal(ms_time_varying_coef_02_mean, ms_time_varying_coef_02_sd);
  }

  // Time-invariant covariate coefficients
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    ms_time_invariant_coef_qr_02 ~ normal(ms_time_invariant_coef_02_mean, ms_time_invariant_coef_02_sd);
  }

  // Multi-level random slopes
  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      ms_sd_level_slope_02[lv] ~ normal(0, ms_sd_level_slope_02_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int lv_start = enabled_level_pos_ms_slope[lv];
      int lv_end = enabled_level_pos_ms_slope[lv + 1] - 1;
      to_vector(ms_raw_level_slope_02[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}

// ============================================================================
// 1→2 TRANSITION PRIORS (Sojourn Time GP)
// ============================================================================

if (need_12_s_gp) {
  ms_log_lambda_gp_12_s_pop_alpha[1] ~ inv_gamma(
    ms_log_lambda_gp_12_s_pop_alpha_alpha, ms_log_lambda_gp_12_s_pop_alpha_beta
  );
  ms_log_lambda_gp_12_s_pop_rho[1] ~ inv_gamma(
    ms_log_lambda_gp_12_s_pop_rho_alpha, ms_log_lambda_gp_12_s_pop_rho_beta
  );
  ms_log_lambda_gp_12_s_pop_intercept[1] ~ normal(
    ms_log_lambda_gp_12_s_pop_intercept_mean, ms_log_lambda_gp_12_s_pop_intercept_sd
  );
  to_vector(ms_log_lambda_gp_12_s_pop_eta) ~ std_normal();

  for (lv in 1:n_levels) {
    ms_log_lambda_gp_12_s_level_alpha[lv] ~ inv_gamma(
      ms_log_lambda_gp_12_s_level_alpha_alpha[lv], ms_log_lambda_gp_12_s_level_alpha_beta[lv]
    );
    ms_log_lambda_gp_12_s_level_rho[lv] ~ inv_gamma(
      ms_log_lambda_gp_12_s_level_rho_alpha[lv], ms_log_lambda_gp_12_s_level_rho_beta[lv]
    );
    ms_log_lambda_gp_12_s_level_intercept_sd[lv] ~ normal(
      0, ms_log_lambda_gp_12_s_level_intercept_sd_sd[lv]
    );
  }
  if (n_enabled_groups_ms_baseline_12_s > 0) {
    ms_raw_log_lambda_gp_12_s_level_intercept ~ std_normal();
    to_vector(ms_log_lambda_gp_12_s_level_eta) ~ std_normal();
  }
}

// ============================================================================
// 1→2 TRANSITION PRIORS (Clock-forward Time GP)
// ============================================================================

if (need_12_t_gp) {
  ms_log_lambda_gp_12_t_pop_alpha[1] ~ inv_gamma(
    ms_log_lambda_gp_12_t_pop_alpha_alpha, ms_log_lambda_gp_12_t_pop_alpha_beta
  );
  ms_log_lambda_gp_12_t_pop_rho[1] ~ inv_gamma(
    ms_log_lambda_gp_12_t_pop_rho_alpha, ms_log_lambda_gp_12_t_pop_rho_beta
  );
  ms_log_lambda_gp_12_t_pop_intercept[1] ~ normal(
    ms_log_lambda_gp_12_t_pop_intercept_mean, ms_log_lambda_gp_12_t_pop_intercept_sd
  );
  to_vector(ms_log_lambda_gp_12_t_pop_eta) ~ std_normal();

  for (lv in 1:n_levels) {
    ms_log_lambda_gp_12_t_level_alpha[lv] ~ inv_gamma(
      ms_log_lambda_gp_12_t_level_alpha_alpha[lv], ms_log_lambda_gp_12_t_level_alpha_beta[lv]
    );
    ms_log_lambda_gp_12_t_level_rho[lv] ~ inv_gamma(
      ms_log_lambda_gp_12_t_level_rho_alpha[lv], ms_log_lambda_gp_12_t_level_rho_beta[lv]
    );
    ms_log_lambda_gp_12_t_level_intercept_sd[lv] ~ normal(
      0, ms_log_lambda_gp_12_t_level_intercept_sd_sd[lv]
    );
  }
  if (n_enabled_groups_ms_baseline_12_t > 0) {
    ms_raw_log_lambda_gp_12_t_level_intercept ~ std_normal();
    to_vector(ms_log_lambda_gp_12_t_level_eta) ~ std_normal();
  }
}

// 1→2 covariate priors
if (enable_ms_12) {
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    ms_time_varying_coef_12 ~ normal(ms_time_varying_coef_12_mean, ms_time_varying_coef_12_sd);
  }

  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    ms_time_invariant_coef_qr_12 ~ normal(ms_time_invariant_coef_12_mean, ms_time_invariant_coef_12_sd);
  }

  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      ms_sd_level_slope_12[lv] ~ normal(0, ms_sd_level_slope_12_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int lv_start = enabled_level_pos_ms_slope[lv];
      int lv_end = enabled_level_pos_ms_slope[lv + 1] - 1;
      to_vector(ms_raw_level_slope_12[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}
