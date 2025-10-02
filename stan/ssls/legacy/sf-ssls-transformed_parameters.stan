vector[n_patients] patient_log_decrease_rate = tr_loc_patient + frac_log_decrease_patient;
vector[n_patients] patient_log_growth_rate   = tr_loc_patient + frac_log_growth_patient;

// Legacy transformed parameters reduced: fraction & init logic migrated to modules.
// Keep growth lag, GP, process noise, and state evolution until those modules are refactored.
vector[n_patients] patient_log_growth_lag_effect = zeros_vector(n_patients);
vector[n_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_patients);

if (!pop_growth_lag_param_only) {
  patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
  patient_log_growth_lag += patient_log_growth_lag_effect;
}
vector[n_patients] patient_tumor_gp_rho = independ_long_process_noise ? zeros_vector(n_patients) : rep_vector(exp(log_pop_tumor_gp_rho), n_patients);  
matrix[n_total_train_visits, 2] states; 

profile("states") {
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_patients] log_patient_tumor_gp_rho_effect; 
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
  if (!independ_long_process_noise && !pop_rho_param_only) {
    log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
    log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
    patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial] + log_patient_tumor_gp_rho_effect);
  }
  states = calc_states(
    patient_visit_pos,
    t_patient_visits,
    patient_tumor_gp_rho,
    delta,
    pop_process_sd,
    independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
    // raw_patient_process_noise,
    independ_long_process_noise, independ_cross_process_noise,
    append_col(init_log_decrease_patient, init_log_growth_patient),
    exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
    rep_vector(0.0001, n_patients), // exp(patient_log_growth_lag), 
    0.0001, // exp(pop_log_growth_transition_rate),
    n_shards, // && !debug,
    0 // debug 
  );
}
