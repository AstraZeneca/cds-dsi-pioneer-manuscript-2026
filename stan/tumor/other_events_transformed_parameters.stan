array[n_causes] row_vector[max_all_t] log_pop_lambda; // log baseline hazard

for (k in 1:n_causes) {
  log_pop_lambda[k] = 
    calc_gp_pred(all_tumor_measure_t, log_lambda_gp_pop_intercept[k], log_lambda_gp_pop_alpha[k], log_lambda_gp_pop_rho[k], delta, log_lambda_gp_pop_eta[k]);
}

// Trial level variation in log baseline hazard
array[n_causes] matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_all_t] log_trial_lambda_residual; 
array[n_causes] matrix[n_trials, max_all_t] log_trial_lambda; 

for (k in 1:n_causes) {
  log_trial_lambda[k] = rep_matrix(log_pop_lambda[k], n_trials); 
}

array[n_causes] vector[add_trial_level_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_intercept;

if (add_trial_level_baseline_hazard) {
  for (k in 1:n_causes) {
    log_lambda_gp_trial_intercept[k] = raw_log_lambda_gp_trial_intercept[k] * log_lambda_gp_trial_intercept_sd[k];
    
    for (s in 1:n_trials) {
      log_trial_lambda_residual[k, s] = calc_gp_pred(
        all_tumor_measure_t,
        log_lambda_gp_trial_intercept[k, s], log_lambda_gp_trial_alpha[k], log_lambda_gp_trial_rho[k], delta, log_lambda_gp_trial_eta[k, s]
      );
      
      log_trial_lambda[k, s] += log_trial_lambda_residual[k, s];
    }
  }
}

// // GLM parameters
// array[n_causes] vector[no_prop_hazard ? 0 : n_surv_covar] surv_covar_effect;  
// 
// // GLM hierarchical parameters
// vector<lower = 0>[add_trial_level_prop_hazard && !no_prop_hazard ? n_surv_covar : 0] surv_covar_trial_sd;
// array[n_causes, add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_surv_covar] raw_surv_covar_trial_coef;

// array[n_causes] vector[no_prop_hazard ? 0 : n_patients] patient_log_crcr_hazard_ratio; // Log proportional hazard
array[n_causes] matrix<upper = 0>[n_train_patients, max_all_t] log_cond_prob_surv;

for (s in 1:n_trials) {
  int patient_start, patient_end;
  (patient_start, patient_end) = get_pos(train_trial_patient_pos, s);
  // int patient_end = patient_pos + n_trial_patients[s] - 1;

  for (k in 1:n_causes) {
    log_cond_prob_surv[k, patient_start:patient_end] = rep_matrix(log_trial_lambda[k, s], get_pos_size(train_trial_patient_pos, s));
      
    // if (!no_prop_hazard) {
    //   patient_log_crcr_hazard_ratio[k, patient_pos:patient_end] =
    //     tumor_sum_covar[patient_pos:patient_end] * crcr_covar_trial_coef[k, s, :n_tumor_covar] + 
    //     covar_design_matrix[patient_pos:patient_end] * crcr_covar_trial_coef[k, s, (n_tumor_covar + 1):];
    //   
    //   log_crcr_cond_prob_surv[k, patient_pos:patient_end] += rep_matrix(patient_log_crcr_hazard_ratio[k, patient_pos:patient_end], max_confresp_week);
    // }
    
    log_cond_prob_surv[k, patient_start:patient_end] = - exp(log_cond_prob_surv[k, patient_start:patient_end]); 
  }
}


matrix[n_train_patients, n_causes] patient_response_lp; 

// Non-target progression
patient_response_lp[, 1] = calc_pch_loglik(
  non_target_pfs[train_patients_pos:train_patients_end], 
  non_target_right_censored[train_patients_pos:train_patients_end], 
  zeros_int_array(n_train_patients),
  0, 
  log_cond_prob_surv[1]
);