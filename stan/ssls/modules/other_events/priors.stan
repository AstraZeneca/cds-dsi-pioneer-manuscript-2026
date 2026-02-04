// ============================================================================
// Other Events Model Priors
// ============================================================================

if (n_causes > 0) {
  // --- Baseline Hazard GP: Population Level ---
  for (k in 1:n_causes) {
    // Use inv_gamma prior to keep alpha away from 0 (avoids GP numerical instability)
    log_lambda_gp_pop_alpha[k] ~ inv_gamma(oe_log_lambda_gp_pop_alpha_alpha[k], oe_log_lambda_gp_pop_alpha_beta[k]);
    log_lambda_gp_pop_rho[k] ~ inv_gamma(oe_log_lambda_gp_pop_rho_alpha[k], oe_log_lambda_gp_pop_rho_beta[k]);
    log_lambda_gp_pop_intercept[k] ~ normal(oe_log_lambda_gp_pop_intercept_mean[k], oe_log_lambda_gp_pop_intercept_sd[k]);
    to_vector(log_lambda_gp_pop_eta[k]) ~ std_normal();
    
    // Baseline Hazard GP: Trial Level
    if (oe_enable_trial_baseline_hazard) {
      // Use inv_gamma prior to keep alpha away from 0 (avoids GP numerical instability)
      log_lambda_gp_trial_alpha[k] ~ inv_gamma(oe_log_lambda_gp_trial_alpha_alpha[k], oe_log_lambda_gp_trial_alpha_beta[k]);
      log_lambda_gp_trial_rho[k] ~ inv_gamma(oe_log_lambda_gp_trial_rho_alpha[k], oe_log_lambda_gp_trial_rho_beta[k]);
      log_lambda_gp_trial_intercept_sd[k] ~ normal(0, oe_log_lambda_gp_trial_intercept_sd_sd[k]);
      raw_log_lambda_gp_trial_intercept[k] ~ std_normal();
      to_vector(log_lambda_gp_trial_eta[k]) ~ std_normal();
    }
  }

  // --- Proportional Hazard: Covariate Effects ---
  for (k in 1:n_causes) {
    // Population-level time-varying tumor coefficients (vector)
    if (oe_enable_pop_tumor_cov) {
      oe_tumor_coef_pop[k] ~ normal(oe_tumor_coef_pop_mean[k],
                                     oe_tumor_coef_pop_sd[k]);
    }

    // Population-level non-tumor covariates (QR space, vector)
    if (oe_enable_pop_cov) {
      oe_covar_coef_qr_pop[k] ~ normal(oe_covar_coef_qr_pop_mean[k],
                                        oe_covar_coef_qr_pop_sd[k]);
    }

    // Multi-level hierarchical priors (non-tumor covariates)
    // Uses enabled position arrays for efficient indexing into compacted parameter arrays
    for (lv in 1:n_levels) {
      // SD hyperpriors are always applied (arrays are always n_levels)
      if (n_covar > 0) {
        oe_sd_level_slope[k, lv] ~ normal(0, oe_sd_level_slope_sd[k, lv]);
      }

      // Slope raw effects - only apply prior to enabled levels
      // (parameter array is sized by enabled groups only)
      if (oe_enable_level_cov[lv] && n_covar > 0) {
        int lv_start = enabled_level_pos_oe_slope[lv];
        int lv_end = enabled_level_pos_oe_slope[lv + 1] - 1;
        to_vector(oe_raw_level_slope[k, lv_start:lv_end, :]) ~ std_normal();
      }
    }
  }
}
