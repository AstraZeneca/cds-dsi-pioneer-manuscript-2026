log_crcr_lambda_gp_alpha ~ normal(0, log_crcr_lambda_gp_alpha_sd);
log_crcr_lambda_gp_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta);
to_vector(log_crcr_lambda_gp_eta) ~ std_normal();
log_crcr_lambda_gp_intercept ~ normal(log_crcr_lambda_gp_intercept_mean, log_crcr_lambda_gp_intercept_sd);

for (k in 1:n_causes) {
  crcr_tumor_stim_pop_coef[, k] ~ normal(0, crcr_tumor_stim_pop_coef_sd);
  crcr_covar_effect[, k] ~ normal(0, crcr_covar_effect_sd);
}

if (add_trial_level) { 
  // log_crcr_lambda_gp_trial_alpha ~ normal(0, log_crcr_lambda_gp_trial_alpha_sd);
  // log_crcr_lambda_gp_trial_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta);
  
  log_crcr_lambda_gp_trial_intercept_sd ~ normal(0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  to_vector(raw_log_crcr_lambda_gp_trial_intercept) ~ std_normal(); 
  
  if (add_trial_level_glm) {
    crcr_covar_trial_sd ~ normal(0, crcr_covar_trial_sd_sd);
    // crcr_covar_trial_corr ~ lkj_corr_cholesky(crcr_covar_trial_corr_eta);
  }
  
  for (s in 1:n_trials) {
    if (add_trial_level_glm) {
      to_vector(raw_crcr_covar_trial_coef[s]) ~ std_normal();
    }
      
    // to_vector(log_crcr_lambda_gp_trial_eta[s]) ~ std_normal();
  }
}
