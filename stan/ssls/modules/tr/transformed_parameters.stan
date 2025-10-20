// tr/transformed_parameters.stan — active linear predictor assembly for total rate module
vector[n_trials] tr_effect_trial_intercept = enable_trial_intercept_tr ? tr_sd_trial_intercept * tr_raw_trial_intercept : rep_vector(0, n_trials);
vector[n_patients] tr_effect_patient_intercept = enable_patient_intercept_tr ? tr_sd_patient_intercept * tr_raw_patient_intercept : rep_vector(0, n_patients);
vector[n_patients] tr_linpred_pop = enable_pop_cov_tr ? (Q_covar_design_matrix * tr_coef_qr_pop) : rep_vector(0, n_patients);

// Trial slope deviations: construct QR-space linear predictor contributions per patient
vector[n_patients] tr_linpred_trial = rep_vector(0, n_patients);

if (enable_trial_cov_tr) {
  matrix[n_trials, n_covar] tr_trial_slope_qr = (tr_raw_trial_slope .* rep_matrix(tr_sd_trial_slope', n_trials));
  tr_linpred_trial = rows_dot_product(Q_covar_design_matrix, tr_trial_slope_qr[patient_trial]);
}

// Patient slope deviations
vector[n_patients] tr_linpred_patient = rep_vector(0, n_patients);

if (enable_patient_cov_tr) {
  matrix[n_patients, n_covar] tr_patient_slope_qr = (tr_raw_patient_slope .* rep_matrix(tr_sd_patient_slope', n_patients));
  tr_linpred_patient = rows_dot_product(Q_covar_design_matrix, tr_patient_slope_qr);
}

vector[n_patients] tr_loc_patient = tr_loc_pop
  + tr_linpred_pop + tr_linpred_trial + tr_linpred_patient
  + tr_effect_trial_intercept[patient_trial]
  + tr_effect_patient_intercept; 

// vector[n_patients] patient_log_total_rate = tr_loc_patient;
