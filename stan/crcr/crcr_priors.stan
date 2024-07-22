log_crcr_lambda_gp_alpha ~ normal(0, log_crcr_lambda_gp_alpha_sd);
log_crcr_lambda_gp_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta);
to_vector(log_crcr_lambda_gp_eta) ~ std_normal();
log_crcr_lambda_gp_intercept ~ normal(log_crcr_lambda_gp_intercept_mean, log_crcr_lambda_gp_intercept_sd);

if (add_trial_level) { 
  log_crcr_lambda_gp_trial_intercept_sd ~ normal(0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  to_vector(raw_log_crcr_lambda_gp_trial_intercept) ~ std_normal(); 
  
  log_crcr_lambda_gp_trial_alpha ~ normal(0, log_crcr_lambda_gp_trial_alpha_sd);
  log_crcr_lambda_gp_trial_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta);
  
  for (s in 1:n_trials) {
    for (k in 1:n_causes) {
      to_vector(log_crcr_lambda_gp_trial_eta[s, k]) ~ std_normal();
    }
  }
}

for (k in 1:n_causes) {
  crcr_tumor_stim_pop_coef[, k] ~ normal(0, crcr_tumor_stim_pop_coef_sd);
 
  crcr_covar_effect[, k] ~ normal(0, crcr_covar_effect_sd);
}