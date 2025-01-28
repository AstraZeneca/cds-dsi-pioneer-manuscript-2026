// Log baseline hazard
matrix[n_base_separate_trials, max_all_t] log_lambda;

for (s in 1:n_base_separate_trials) {
  log_lambda[s] = calc_gp_pred(pfs_range, log_lambda_gp_intercept[s], log_lambda_gp_alpha[s], log_lambda_gp_rho[s], delta, log_lambda_gp_eta[s]);
}

matrix[n_trials, max_all_t] log_trial_lambda;

if (separate_baseline_hazard) {
  log_trial_lambda = log_lambda; 
} else {
  log_trial_lambda = rep_matrix(log_lambda[1], n_trials);
}

// Trial-level variation in log baseline hazard 
vector[add_trial_level_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_intercept;
matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_all_t] log_trial_lambda_residual;

if (add_trial_level_baseline_hazard) {
  // We're using uncentered hierarchical effects here to reduce divergent transitions
  log_lambda_gp_trial_intercept = raw_log_lambda_gp_trial_intercept * log_lambda_gp_trial_intercept_sd;
  
  for (s in 1:n_trials) {
    log_trial_lambda_residual[s] = calc_gp_pred(pfs_range, log_lambda_gp_trial_intercept[s], log_lambda_gp_trial_alpha, log_lambda_gp_trial_rho, delta, log_lambda_gp_trial_eta[s]); 
    log_trial_lambda[s] += log_trial_lambda_residual[s];
  }
}