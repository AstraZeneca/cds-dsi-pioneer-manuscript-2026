// frac/transformed_parameters.stan — active fraction module transformed parameters
vector[n_train_trials] frac_effect_trial_intercept = enable_trial_intercept_frac ? frac_sd_trial_intercept * frac_raw_trial_intercept : rep_vector(0, n_train_trials);
vector[n_train_patients] frac_effect_patient_intercept = enable_patient_intercept_frac ? frac_sd_patient_intercept * frac_raw_patient_intercept : rep_vector(0, n_train_patients);
vector[n_train_patients] frac_linpred_pop = enable_pop_cov_frac ? (Q_covar_design_matrix * frac_coef_qr_pop) : rep_vector(0, n_train_patients);

// (Optional future) slope deviations analogous to tr
vector[n_train_patients] frac_linpred_patient_dev = rep_vector(0, n_train_patients);

if (enable_trial_cov_frac) {
  matrix[n_train_trials, n_covar] frac_trial_slope_qr = frac_raw_trial_slope .* rep_matrix(frac_sd_trial_slope', n_train_trials);
  frac_linpred_patient_dev += rows_dot_product(Q_covar_design_matrix, frac_trial_slope_qr[train_patient_trial]);
}
if (enable_patient_cov_frac) {
  matrix[n_train_patients, n_covar] frac_patient_slope_qr = frac_raw_patient_slope .* rep_matrix(frac_sd_patient_slope', n_train_patients);
  frac_linpred_patient_dev += rows_dot_product(Q_covar_design_matrix, frac_patient_slope_qr);
}

vector[n_train_patients] patient_decrease_frac_logit = pop_decrease_frac_logit
  + frac_effect_trial_intercept[train_patient_trial] 
  + frac_effect_patient_intercept
  + frac_linpred_pop + frac_linpred_patient_dev;

vector[n_train_patients] patient_log_decrease_frac = -log1p_exp(-patient_decrease_frac_logit);
vector[n_train_patients] patient_log_growth_frac   = -log1p_exp(patient_decrease_frac_logit);
