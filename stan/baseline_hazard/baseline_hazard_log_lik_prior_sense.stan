if (prior_sense) {
  for (s in 1:(separate_baseline_hazard ? n_trials : 1)) {
    lprior += normal_lpdf(log_lambda_gp_alpha[s] | 0, log_lambda_gp_alpha_sd[s]) + inv_gamma_lpdf(log_lambda_gp_rho[s] | log_lambda_gp_rho_alpha[s], log_lambda_gp_rho_beta[s]) +
      normal_lpdf(log_lambda_gp_intercept[s] | log_lambda_gp_intercept_mean[s], log_lambda_gp_intercept_sd[s]);
  }
  
  if (add_trial_level_prop_hazard) {
    lprior += normal_lpdf(log_lambda_gp_trial_intercept_sd | 0, log_lambda_gp_trial_intercept_sd_sd) + normal_lpdf(log_lambda_gp_trial_alpha | 0, log_lambda_gp_trial_alpha_sd) +
      inv_gamma_lpdf(log_lambda_gp_trial_rho | log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
  }
}