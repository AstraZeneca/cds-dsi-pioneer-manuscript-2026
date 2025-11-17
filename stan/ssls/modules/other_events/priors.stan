// ============================================================================
// Other Events Model Priors
// ============================================================================

if (n_causes > 0) {
  // --- Baseline Hazard GP: Population Level ---
  for (k in 1:n_causes) {
    log_lambda_gp_pop_alpha[k] ~ normal(0, oe_log_lambda_gp_pop_alpha_sd[k]);
    log_lambda_gp_pop_rho[k] ~ inv_gamma(oe_log_lambda_gp_pop_rho_alpha[k], oe_log_lambda_gp_pop_rho_beta[k]);
    log_lambda_gp_pop_intercept[k] ~ normal(oe_log_lambda_gp_pop_intercept_mean[k], oe_log_lambda_gp_pop_intercept_sd[k]);
    to_vector(log_lambda_gp_pop_eta[k]) ~ std_normal();
    
    // Baseline Hazard GP: Trial Level
    if (oe_enable_trial_baseline_hazard) {
      log_lambda_gp_trial_alpha[k] ~ normal(0, oe_log_lambda_gp_trial_alpha_sd[k]);
      log_lambda_gp_trial_rho[k] ~ inv_gamma(oe_log_lambda_gp_trial_rho_alpha[k], oe_log_lambda_gp_trial_rho_beta[k]);
      log_lambda_gp_trial_intercept_sd[k] ~ normal(0, oe_log_lambda_gp_trial_intercept_sd_sd[k]);
      raw_log_lambda_gp_trial_intercept[k] ~ std_normal();
      to_vector(log_lambda_gp_trial_eta[k]) ~ std_normal();
    }
  }

  // --- Proportional Hazard: Covariate Effects ---
  if (oe_enable_pop_cov || oe_enable_pop_tumor_cov || oe_enable_trial_cov) {
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
      
      // Trial-level hierarchical priors (non-tumor covariates only)
      if (oe_enable_trial_cov) {
        oe_sd_trial_slope[k] ~ normal(0, oe_sd_trial_slope_sd[k]);
        to_vector(oe_raw_trial_slope[k]) ~ std_normal();
      }
    }
  }
}
