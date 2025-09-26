// init/transformed_parameters.stan — active initial proportion module
vector[n_train_trials] init_effect_trial_intercept = enable_trial_intercept_init ? init_sd_trial_intercept * init_raw_trial_intercept : rep_vector(0, n_train_trials);
vector[n_train_patients] init_effect_patient_intercept = enable_patient_intercept_init ? init_sd_patient_intercept * init_raw_patient_intercept : rep_vector(0, n_train_patients);
vector[n_train_patients] init_linpred_pop = enable_pop_cov_init ? (Q_covar_design_matrix * init_coef_qr_pop) : rep_vector(0, n_train_patients);

vector[n_train_patients] init_linpred_patient_dev = rep_vector(0, n_train_patients);

if (enable_trial_cov_init) {
  matrix[n_train_trials, n_covar] init_trial_slope_qr = (init_raw_trial_slope .* rep_matrix(init_sd_trial_slope', n_train_trials));
  init_linpred_patient_dev += rows_dot_product(Q_covar_design_matrix, init_trial_slope_qr[train_patient_trial]);
}
if (enable_patient_cov_init) {
  matrix[n_train_patients, n_covar] init_patient_slope_qr = (init_raw_patient_slope .* rep_matrix(init_sd_patient_slope', n_train_patients));
  init_linpred_patient_dev += rows_dot_product(Q_covar_design_matrix, init_patient_slope_qr);
}

vector[n_train_patients] patient_decrease_prop_logis = pop_decrease_prop_logis
  + init_effect_trial_intercept[train_patient_trial]
  + init_effect_patient_intercept 
  + init_linpred_pop + init_linpred_patient_dev;

vector[n_train_patients] patient_log_decrease_prop = -log1p_exp(- patient_decrease_prop_logis);
vector[n_train_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
