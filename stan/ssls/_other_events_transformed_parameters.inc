// Include fragment: transformed parameters for other events baseline hazard GP
array[n_causes] row_vector[max_all_t] log_pop_lambda;
for (k in 1:n_causes) {
  log_pop_lambda[k] = calc_gp_pred(all_tumor_measure_t, log_lambda_gp_pop_intercept[k], log_lambda_gp_pop_alpha[k], log_lambda_gp_pop_rho[k], delta, log_lambda_gp_pop_eta[k]);
}
array[n_causes] matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_all_t] log_trial_lambda_residual; 
array[n_causes] matrix[n_trials, max_all_t] log_trial_lambda; 
for (k in 1:n_causes) log_trial_lambda[k] = rep_matrix(log_pop_lambda[k], n_trials); 
array[n_causes] vector[add_trial_level_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_intercept;
if (add_trial_level_baseline_hazard) {
  for (k in 1:n_causes) {
    log_lambda_gp_trial_intercept[k] = raw_log_lambda_gp_trial_intercept[k] * log_lambda_gp_trial_intercept_sd[k];
    for (s in 1:n_trials) {
      log_trial_lambda_residual[k, s] = calc_gp_pred(all_tumor_measure_t, log_lambda_gp_trial_intercept[k, s], log_lambda_gp_trial_alpha[k], log_lambda_gp_trial_rho[k], delta, log_lambda_gp_trial_eta[k, s]);
      log_trial_lambda[k, s] += log_trial_lambda_residual[k, s];
    }
  }
}
array[n_causes] matrix<upper = 0>[n_train_patients, max_all_t] log_cond_prob_surv;
for (s in 1:n_trials) {
  int patient_start, patient_end; (patient_start, patient_end) = get_pos(train_trial_patient_pos, s);
  for (k in 1:n_causes) {
    log_cond_prob_surv[k, patient_start:patient_end] = rep_matrix(log_trial_lambda[k, s], get_pos_size(train_trial_patient_pos, s));
    log_cond_prob_surv[k, patient_start:patient_end] = - exp(log_cond_prob_surv[k, patient_start:patient_end]); 
  }
}
