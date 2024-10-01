log_lambda_gp_alpha ~ normal(0, log_lambda_gp_alpha_sd);
log_lambda_gp_rho ~ inv_gamma(log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
log_lambda_gp_eta ~ std_normal();
log_lambda_gp_intercept ~ normal(log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd);

log_lambda_gp_trial_intercept_sd ~ normal(0, log_lambda_gp_trial_intercept_sd_sd);
raw_log_lambda_gp_trial_intercept ~ std_normal(); 

// TODO separate hyperparam for these parameters 
log_lambda_gp_trial_alpha ~ normal(0, log_lambda_gp_trial_alpha_sd);
log_lambda_gp_trial_rho ~ inv_gamma(log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);

if (add_trial_level) { 
  for (s in 1:n_trials) {
    log_lambda_gp_trial_eta[s] ~ std_normal();
  }
}