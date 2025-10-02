// frac/transformed_parameters.stan — active fraction module transformed parameters
vector[n_train_trials] frac_effect_trial_intercept = enable_trial_intercept_frac ? frac_sd_trial_intercept * frac_raw_trial_intercept : rep_vector(0, n_train_trials);
vector[n_patients] frac_effect_patient_intercept = enable_patient_intercept_frac ? frac_sd_patient_intercept * frac_raw_patient_intercept : rep_vector(0, n_patients);
vector[n_patients] frac_linpred_pop = enable_pop_cov_frac ? (Q_covar_design_matrix * frac_coef_qr_pop) : rep_vector(0, n_patients);

// Trial slope deviations: construct QR-space linear predictor contributions per patient
vector[n_patients] frac_linpred_trial = rep_vector(0, n_patients);

if (enable_trial_cov_frac) {
  matrix[n_train_trials, n_covar] frac_trial_slope_qr = frac_raw_trial_slope .* rep_matrix(frac_sd_trial_slope', n_train_trials);
  frac_linpred_trial = rows_dot_product(Q_covar_design_matrix, frac_trial_slope_qr[patient_trial]);
}

// Patient slope deviations
vector[n_patients] frac_linpred_patient = rep_vector(0, n_patients);

if (enable_patient_cov_frac) {
  matrix[n_patients, n_covar] frac_patient_slope_qr = frac_raw_patient_slope .* rep_matrix(frac_sd_patient_slope', n_patients);
  frac_linpred_patient = rows_dot_product(Q_covar_design_matrix, frac_patient_slope_qr);
}

vector[n_patients] frac_logit_loc_patient = frac_logit_loc_pop
  + frac_linpred_pop + frac_linpred_trial + frac_linpred_patient
  + frac_effect_trial_intercept[patient_trial] 
  + frac_effect_patient_intercept;

vector[n_patients] frac_log_decrease_patient = log_inv_logit(frac_logit_loc_patient);
vector[n_patients] frac_log_growth_patient   = log1m_inv_logit(frac_logit_loc_patient);
