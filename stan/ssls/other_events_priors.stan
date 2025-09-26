// Relocated other_events_priors.stan
for (k in 1:n_causes) {
  log_lambda_gp_pop_alpha[k] ~ normal(0, log_lambda_gp_pop_alpha_sd[k]);
  log_lambda_gp_pop_rho[k] ~ inv_gamma(log_lambda_gp_pop_rho_alpha[k], log_lambda_gp_pop_rho_beta[k]);
  log_lambda_gp_pop_intercept[k] ~ normal(log_lambda_gp_pop_intercept_mean[k], log_lambda_gp_pop_intercept_sd[k]);
  to_vector(log_lambda_gp_pop_eta[k]) ~ std_normal();
  if (add_trial_level_baseline_hazard) {
    log_lambda_gp_trial_alpha[k] ~ normal(0, log_lambda_gp_trial_alpha_sd[k]);
    log_lambda_gp_trial_rho[k] ~ inv_gamma(log_lambda_gp_trial_rho_alpha[k], log_lambda_gp_trial_rho_beta[k]);
    log_lambda_gp_trial_intercept_sd[k] ~ normal(0, log_lambda_gp_trial_intercept_sd_sd[k]);
    raw_log_lambda_gp_trial_intercept[k] ~ std_normal();
    to_vector(log_lambda_gp_trial_eta[k]) ~ std_normal();
  }
}
