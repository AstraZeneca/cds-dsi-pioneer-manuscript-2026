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

vector[n_patients] patient_decrease_rate = exp(patient_log_decrease_rate);
vector[n_patients] patient_growth_rate = exp(patient_log_growth_rate);
    

profile("states") {
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_patients] log_patient_tumor_gp_rho_effect; 
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
  if (!independ_long_process_noise && !pop_rho_param_only) {
    log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
    log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
    patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial] + log_patient_tumor_gp_rho_effect);
  }

  // states = calc_states(
  //   patient_visit_pos,
  //   t_patient_visits,
  //   patient_tumor_gp_rho,
  //   delta,
  //   pop_process_sd,
  //   independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
  //   // raw_patient_process_noise,
  //   independ_long_process_noise, independ_cross_process_noise,
  //   append_col(init_log_decrease_patient, init_log_growth_patient),
  //   exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
  //   rep_vector(0.0001, n_patients), // exp(patient_log_growth_lag), 
  //   0.0001, // exp(pop_log_growth_transition_rate),
  //   n_shards, // && !debug,
  //   0 // debug 
  // );

  // ============================================================================
  // Batched matrix computation: V * D' and V * G'
  // Computes cumulative states for ALL patients in single matrix operation
  // ============================================================================
  
  // Step 2: Build increment matrices [n_patients × time_range]
  // Using unit time steps (delta_t = 1 day)
  // D[i,t] = -d[i], G[i,t] = g[i]
  matrix[n_patients, max_unique_visit] D_increments = rep_matrix(-patient_decrease_rate, max_unique_visit);
  matrix[n_patients, max_unique_visit] G_increments = rep_matrix(patient_growth_rate, max_unique_visit);
  
  // Step 3: Single batched matrix multiply for ALL patients
  // Step 4: Add initial states (vectorized broadcast)
  matrix[n_pop_unique_visits, n_patients] cumsum_d = visit_cumsum_mat * D_increments' + rep_matrix(init_log_decrease_patient', n_pop_unique_visits);
  matrix[n_pop_unique_visits, n_patients] cumsum_g = visit_cumsum_mat * G_increments' + rep_matrix(init_log_growth_patient', n_pop_unique_visits);
  
  profile("extract states") {
  // Step 5: Extract states at actual patient visit times
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    states[visit_start:visit_end, 1] = cumsum_d[patient2pop_unique_visit_idx[visit_start:visit_end], i];
    states[visit_start:visit_end, 2] = cumsum_g[patient2pop_unique_visit_idx[visit_start:visit_end], i];
  }}
}
