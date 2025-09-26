// tr/transformed_parameters.stan — active linear predictor assembly for total rate module
vector[n_train_trials] trial_log_total_rate_effect = enable_trial_intercept_tr ? tr_sd_trial_intercept * tr_raw_trial_intercept : rep_vector(0, n_train_trials);
vector[n_train_patients] patient_log_total_rate_effect = enable_patient_intercept_tr ? tr_sd_patient_intercept * tr_raw_patient_intercept : rep_vector(0, n_train_patients);
vector[n_train_patients] tr_linpred_pop = enable_pop_cov_tr ? (Q_covar_design_matrix * tr_coef_qr_pop) : rep_vector(0, n_train_patients);

// Trial slope deviations: construct QR-space linear predictor contributions per patient
vector[n_train_patients] tr_linpred_trial = rep_vector(0, n_train_patients);

if (enable_trial_cov_tr) {
  matrix[n_train_trials, n_covar] tr_trial_slope_qr = (tr_raw_trial_slope .* rep_matrix(tr_sd_trial_slope', n_train_trials));
  tr_linpred_trial = rows_dot_product(Q_covar_design_matrix, tr_trial_slope_qr[train_patient_trial]);
}

// Patient slope deviations
vector[n_train_patients] tr_linpred_patient_dev = rep_vector(0, n_train_patients);

if (enable_patient_cov_tr) {
  matrix[n_train_patients, n_covar] tr_patient_slope_qr = (tr_raw_patient_slope .* rep_matrix(tr_sd_patient_slope', n_train_patients));
  tr_linpred_patient_dev = rows_dot_product(Q_covar_design_matrix, tr_patient_slope_qr);
}

vector[n_train_patients] tr_loc_patient = pop_log_total_rate
  + tr_linpred_pop + tr_linpred_trial + tr_linpred_patient_dev
  + trial_log_total_rate_effect[train_patient_trial]
  + patient_log_total_rate_effect; 

// vector[n_train_patients] patient_log_total_rate = tr_loc_patient;
