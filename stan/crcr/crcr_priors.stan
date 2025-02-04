for (s in 1:n_base_separate_trials) {
  log_crcr_lambda_gp_alpha[s] ~ normal(0, log_crcr_lambda_gp_alpha_sd[s]);
  to_vector(log_crcr_lambda_gp_rho[s]) ~ inv_gamma(log_crcr_lambda_gp_rho_alpha[s], log_crcr_lambda_gp_rho_beta[s]);
  
  log_crcr_lambda_gp_intercept[s] ~ normal(log_crcr_lambda_gp_intercept_mean[s], log_crcr_lambda_gp_intercept_sd[s]);
  to_vector(log_crcr_lambda_gp_eta[s]) ~ std_normal();
}

for (s in 1:n_prop_separate_trials) {
  for (k in 1:n_causes) {
    crcr_tumor_stim_pop_coef[s, , k] ~ normal(0, crcr_tumor_stim_pop_coef_sd[s]);
    crcr_covar_effect[s, , k] ~ normal(crcr_covar_effect_mean[s], crcr_covar_effect_sd[s]);
  }
}

if (add_trial_level_baseline_hazard) { 
  // TODO Look further into how to make this fully hierarchical
  log_crcr_lambda_gp_trial_alpha ~ normal(0, log_crcr_lambda_gp_trial_alpha_sd);
  log_crcr_lambda_gp_trial_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha[1], log_crcr_lambda_gp_rho_beta[1]);
  
  log_crcr_lambda_gp_trial_intercept_sd ~ normal(0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  to_vector(raw_log_crcr_lambda_gp_trial_intercept) ~ std_normal(); 
  
  for (s in 1:n_trials) {
    to_vector(log_crcr_lambda_gp_trial_eta[s]) ~ std_normal();
  }
}

if (add_trial_level_prop_hazard && !no_prop_hazard) {
  crcr_covar_trial_sd ~ normal(0, crcr_covar_trial_sd_sd);
  L_crcr_covar_trial_corr ~ lkj_corr_cholesky(crcr_covar_trial_corr_eta);
  
  for (s in 1:n_trials) {
    to_vector(raw_crcr_covar_trial_coef[s]) ~ std_normal();
  }
}
