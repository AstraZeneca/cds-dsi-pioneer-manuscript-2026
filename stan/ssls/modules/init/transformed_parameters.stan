// init/transformed_parameters.stan — active initial proportion module
vector[n_train_trials] init_effect_trial_intercept = enable_trial_intercept_init ? init_sd_trial_intercept * init_raw_trial_intercept : rep_vector(0, n_train_trials);
vector[n_patients] init_effect_patient_intercept = enable_patient_intercept_init ? init_sd_patient_intercept * init_raw_patient_intercept : rep_vector(0, n_patients);
vector[n_patients] init_linpred_pop = enable_pop_cov_init ? (Q_covar_design_matrix * init_coef_qr_pop) : rep_vector(0, n_patients);

// Trial slope deviations: construct QR-space linear predictor contributions per patient
vector[n_patients] init_linpred_trial = rep_vector(0, n_patients);

if (enable_trial_cov_init) {
  matrix[n_train_trials, n_covar] init_trial_slope_qr = (init_raw_trial_slope .* rep_matrix(init_sd_trial_slope', n_train_trials));
  init_linpred_trial = rows_dot_product(Q_covar_design_matrix, init_trial_slope_qr[patient_trial]);
}

// Patient slope deviations
vector[n_patients] init_linpred_patient = rep_vector(0, n_patients);

if (enable_patient_cov_init) {
  matrix[n_patients, n_covar] init_patient_slope_qr = (init_raw_patient_slope .* rep_matrix(init_sd_patient_slope', n_patients));
  init_linpred_patient = rows_dot_product(Q_covar_design_matrix, init_patient_slope_qr);
}

vector[n_patients] init_logit_loc_patient = init_logit_loc_pop
  + init_linpred_pop + init_linpred_trial + init_linpred_patient
  + init_effect_trial_intercept[patient_trial]
  + init_effect_patient_intercept;

vector[n_patients] init_log_decrease_patient = log_inv_logit(init_logit_loc_patient);
vector[n_patients] init_log_growth_patient   = log1m_inv_logit(init_logit_loc_patient);
