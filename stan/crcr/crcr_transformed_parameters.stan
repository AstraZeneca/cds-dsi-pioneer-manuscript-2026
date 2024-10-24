array[separate_baseline_hazard ? n_trials : 1] matrix[max_confresp_week, n_causes] log_crcr_lambda; 

profile("crcr population baseline hazards") {
  for (s in 1:(separate_baseline_hazard ? n_trials : 1)) {
    for (k in 1:n_causes) {
      log_crcr_lambda[s, , k] = 
        calc_gp_pred(confresp_range, log_crcr_lambda_gp_intercept[k, s], log_crcr_lambda_gp_alpha[s, k], log_crcr_lambda_gp_rho[s, k], delta, log_crcr_lambda_gp_eta[s, , k]);
    }
  }
}

array[add_trial_level_baseline_hazard ? n_trials : 0] matrix[max_confresp_week, n_causes] log_crcr_trial_lambda_residual;
array[n_trials] matrix[max_confresp_week, n_causes] log_crcr_trial_lambda = separate_baseline_hazard ? log_crcr_lambda : rep_array(log_crcr_lambda[1], n_trials); 

matrix[add_trial_level_baseline_hazard ? n_trials : 0, n_causes] log_crcr_lambda_gp_trial_intercept;
array[n_trials] matrix[n_tumor_covar + n_covar, n_causes] crcr_covar_trial_coef;
array[add_trial_level_prop_hazard ? n_trials : 0] matrix[n_tumor_covar + n_covar, n_causes] crcr_covar_trial_coef_residual;

if (add_trial_level_baseline_hazard) {
  log_crcr_lambda_gp_trial_intercept = diag_post_multiply(raw_log_crcr_lambda_gp_trial_intercept, log_crcr_lambda_gp_trial_intercept_sd);
    
  for (s in 1:n_trials) {
    for (k in 1:n_causes) {
      log_crcr_trial_lambda_residual[s, , k] = calc_gp_pred(
        confresp_range,
        log_crcr_lambda_gp_trial_intercept[s, k], log_crcr_lambda_gp_trial_alpha[k], log_crcr_lambda_gp_trial_rho[k], delta, log_crcr_lambda_gp_trial_eta[s, , k]
      );
      
      log_crcr_trial_lambda[s, , k] += log_crcr_trial_lambda_residual[s, , k];
    }
  }
}

for (s in 1:n_trials) {
  if (add_trial_level_prop_hazard) {
    crcr_covar_trial_coef_residual[s] = diag_pre_multiply(crcr_covar_trial_sd, L_crcr_covar_trial_corr) * raw_crcr_covar_trial_coef[s];
    crcr_covar_trial_coef[s] = append_row(crcr_tumor_stim_pop_coef[1], crcr_covar_effect[1]) + crcr_covar_trial_coef_residual[s];
  } else if (separate_prop_hazard) {
    crcr_covar_trial_coef[s] = append_row(crcr_tumor_stim_pop_coef[s], crcr_covar_effect[s]);
  } else {
    crcr_covar_trial_coef[s] = append_row(crcr_tumor_stim_pop_coef[1], crcr_covar_effect[1]);
  }
}

matrix[n_patients, n_causes] patient_log_crcr_hazard_ratio;

matrix<upper = 0>[n_crcr_time_periods, n_causes] log_crcr_cond_prob_surv;

profile("log_crcr_cond_prob_surv") { // Calculate patient-interval conditional probability of survival (without a classified response) 
  int patient_pos = 1;
  
  for (s in 1:n_trials) {
    int patient_end = patient_pos + n_trial_patients[s] - 1;
   
    patient_log_crcr_hazard_ratio[patient_pos:patient_end] =
      tumor_sum_covar[patient_pos:patient_end] * crcr_covar_trial_coef[s, :n_tumor_covar] +
      covar_design_matrix[patient_pos:patient_end] * crcr_covar_trial_coef[s, (n_tumor_covar + 1):];
    
    for (i in patient_pos:patient_end) {
      int confresp_interval_pos = patient_conf_resp_interval_pos[i];
      int confresp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
      
      log_crcr_cond_prob_surv[confresp_interval_pos:confresp_interval_end] = 
        - exp(log_crcr_trial_lambda[s, 1:max_confresp_week] + rep_matrix(patient_log_crcr_hazard_ratio[i], max_confresp_week));
    }
    
    patient_pos = patient_end + 1;
  }
}