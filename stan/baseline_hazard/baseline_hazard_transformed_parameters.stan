// Log baseline hazard
vector[max_all_t] log_lambda = calc_gp_pred(pfs_range, log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, delta, log_lambda_gp_eta);

// Trial-level variation in log baseline hazard 
vector[add_trial_level ? n_trials : 0] log_lambda_gp_trial_intercept;
array[add_trial_level ? n_trials : 0] vector[max_all_t] log_trial_lambda_residual;
array[n_trials] vector[max_all_t] log_trial_lambda = rep_array(log_lambda, n_trials); // Need to initialize 

if (add_trial_level) {
  // We're using uncentered hierarchical effects here to reduce divergent transitions
  log_lambda_gp_trial_intercept = raw_log_lambda_gp_trial_intercept * log_lambda_gp_trial_intercept_sd;
  
  for (s in 1:n_trials) {
    log_trial_lambda_residual[s] = 
      calc_gp_pred(pfs_range, log_lambda_gp_trial_intercept[s], log_lambda_gp_trial_alpha, log_lambda_gp_trial_rho, delta, log_lambda_gp_trial_eta[s]); 
    log_trial_lambda[s] = log_lambda + log_trial_lambda_residual[s];
  }
}