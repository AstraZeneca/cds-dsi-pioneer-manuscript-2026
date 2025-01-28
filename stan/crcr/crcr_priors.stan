for (k in 1:n_causes) {
  log_crcr_lambda_gp_alpha[k] ~ normal(0, log_crcr_lambda_gp_alpha_sd);
  log_crcr_lambda_gp_rho[k] ~ inv_gamma(log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta);
  
  for (s in 1:n_base_separate_trials) {
    log_crcr_lambda_gp_intercept[k] ~ normal(log_crcr_lambda_gp_intercept_mean[s, k], log_crcr_lambda_gp_intercept_sd[s]);
  }
  
  to_vector(log_crcr_lambda_gp_eta[k]) ~ std_normal();
  
  for (s in 1:n_prop_separate_trials) {
    crcr_tumor_stim_pop_coef[k, s] ~ normal(0, crcr_tumor_stim_pop_coef_sd[s]);
    crcr_covar_effect[k, s] ~ normal(crcr_covar_effect_mean[s], crcr_covar_effect_sd[s]);
  }
}

if (add_trial_level_baseline_hazard) { 
  // TODO Look further into how to make this fully hierarchical
  log_crcr_lambda_gp_trial_alpha ~ normal(0, log_crcr_lambda_gp_trial_alpha_sd);
  log_crcr_lambda_gp_trial_rho ~ inv_gamma(log_crcr_lambda_gp_rho_alpha[1], log_crcr_lambda_gp_rho_beta[1]);
  
  log_crcr_lambda_gp_trial_intercept_sd ~ normal(0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  
  for (k in 1:n_causes) {
    raw_log_crcr_lambda_gp_trial_intercept[k] ~ std_normal(); 
    to_vector(log_crcr_lambda_gp_trial_eta[k]) ~ std_normal();
  }
}

if (add_trial_level_prop_hazard && !no_prop_hazard) {
  crcr_covar_trial_sd ~ normal(0, crcr_covar_trial_sd_sd);
  L_crcr_covar_trial_corr ~ lkj_corr_cholesky(crcr_covar_trial_corr_eta);
 
  for (k in 1:n_causes) { 
    for (s in 1:n_trials) {
      raw_crcr_covar_trial_coef[k, s] ~ std_normal();
    }
  }
}
