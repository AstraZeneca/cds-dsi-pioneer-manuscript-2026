vector[n_trials] trial_tumor_intercept_effect = zeros_vector(n_trials);
vector[n_patients] patient_tumor_intercept_effect = raw_patient_tumor_intercept_effect * patient_tumor_intercept_sd;

vector[n_trials] log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
vector[n_patients] log_patient_tumor_gp_rho_effect = raw_log_patient_tumor_gp_rho_effect * log_patient_tumor_gp_rho_sd;

if (add_trial_level_tumor_intercept) {
  trial_tumor_intercept_effect = raw_trial_tumor_intercept_effect * trial_tumor_intercept_sd;
}

if (add_trial_level_tumor_gp_param) {
  log_trial_tumor_gp_rho_effect = raw_log_trial_tumor_gp_rho_effect * log_trial_tumor_gp_rho_sd;
}

vector[n_patients] patient_tumor_intercept = pop_tumor_intercept + trial_tumor_intercept_effect[patient_trial] + patient_tumor_intercept_effect;
vector<lower = 0>[n_patients] patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial] + log_patient_tumor_gp_rho_effect);

vector[sum(n_patient_unique_visits)] patient_obs_tumor_gp;

for (i in 1:n_patients) {
  int patient_gp_pos, patient_gp_end;
  (patient_gp_pos, patient_gp_end) = get_pos(patient_unique_visits_pos, i);
  
  patient_obs_tumor_gp[patient_gp_pos:patient_gp_end] = calc_gp_pred(
  // patient_obs_tumor_gp[patient_gp_pos:patient_gp_end] = ncp_gp_matern52(
    all_tumor_measure_t[get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i)], 
    pop_tumor_gp_alpha, patient_tumor_gp_rho[i], delta, get_sub_vector(patient_tumor_gp_eta, patient_unique_visits_pos, i) 
  );
}