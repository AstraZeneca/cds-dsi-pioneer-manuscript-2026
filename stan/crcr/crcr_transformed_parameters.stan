// Log baseline hazard
matrix[max_confresp_week, n_causes] log_crcr_lambda; 

// Trial-level variation in log baseline hazard 
matrix[add_trial_level ? n_trials : 0, n_causes] log_crcr_lambda_gp_trial_intercept;

for (k in 1:n_causes) {
  log_crcr_lambda[, k] = 
    calc_gp_pred(confresp_range, log_crcr_lambda_gp_intercept[k], log_crcr_lambda_gp_alpha[k], log_crcr_lambda_gp_rho[k], delta, log_crcr_lambda_gp_eta[, k]);
}

array[n_trials] matrix[max_confresp_week, n_causes] log_crcr_trial_lambda = rep_array(log_crcr_lambda, n_trials); // Need to initialize 

if (add_trial_level) {
  for (k in 1:n_causes) {
    // We're using uncentered hierarchical effects here to reduce divergent transitions
    log_crcr_lambda_gp_trial_intercept[, k] = raw_log_crcr_lambda_gp_trial_intercept[, k] * log_crcr_lambda_gp_trial_intercept_sd[k];
    
    for (s in 1:n_trials) {
      log_crcr_trial_lambda[s, , k] = log_crcr_lambda[, k] +  
        calc_gp_pred(
          confresp_range, 
          log_crcr_lambda_gp_trial_intercept[s, k], log_crcr_lambda_gp_trial_alpha[k], log_crcr_lambda_gp_trial_rho[k], delta, log_crcr_lambda_gp_trial_eta[s, , k]
        ); 
    }
  }
}

matrix[n_patients, n_causes] patient_log_crcr_hazard_ratio = tumor_sum_covar * crcr_tumor_stim_pop_coef + covar_design_matrix * crcr_covar_effect;  
matrix<upper = 0>[n_crcr_time_periods, n_causes] log_crcr_cond_prob_surv;

{ // Calculate patient-interval conditional probability of survival (without a classified response) 
  int confresp_interval_pos = 1;
  
  for (i in 1:n_patients) {
    int n_intervals = max_confresp_week; 
    int confresp_interval_end = confresp_interval_pos + n_intervals - 1;
    
    for (k in 1:n_causes) {
      log_crcr_cond_prob_surv[confresp_interval_pos:confresp_interval_end, k] = 
        - exp(log_crcr_trial_lambda[patient_trial[i], 1:n_intervals, k] + patient_log_crcr_hazard_ratio[i, k]);
    }

    confresp_interval_pos = confresp_interval_end + 1;
  }
}