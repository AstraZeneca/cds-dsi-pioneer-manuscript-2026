vector[n_train_trials] trial_log_net_rate_effect = zeros_vector(n_train_trials);
vector[n_train_patients] patient_log_net_rate_effect = zeros_vector(n_train_patients);
vector[n_train_patients] patient_log_net_rate = rep_vector(pop_log_net_rate, n_train_patients);
// vector[n_train_patients] patient_log_rate_ratio_effect = zeros_vector(n_train_patients);
vector[n_train_patients] patient_log_rate_ratio = rep_vector(pop_log_rate_ratio, n_train_patients);
  
// Calculate linear predictors for rates
vector[n_train_patients] patient_log_net_rate_linpred = Q_covar_design_matrix * QR_pop_log_net_rate_coef;
  
matrix[n_train_trials, n_covar] trial_log_net_rate_coef = rep_matrix(0, n_train_trials, n_covar);
  
if (!pop_covar_coef_only) {
  // Trial-level covariate effects
  trial_log_net_rate_coef = rep_matrix(trial_log_net_rate_coef_sd, n_train_trials) .* raw_trial_log_net_rate_coef;
  
  patient_log_net_rate_linpred += rows_dot_product(Q_covar_design_matrix,
                                                   trial_log_net_rate_coef[train_patient_trial]);
  
  // Patient-level covariate effects (if you want this level of complexity)
  // matrix[n_train_patients, n_covar] patient_log_net_rate_coef = patient_log_net_rate_coef_sd * raw_patient_log_net_rate_coef;
  // patient_log_net_rate_lp += rows_dot_product(covar_design_matrix[train_patients_pos:train_patients_end], patient_log_net_rate_coef);
}

// Add hierarchical covariate effects if needed
if (!pop_rates_param_only) {
  if (add_trial_level_net_rate) {
    trial_log_net_rate_effect = trial_log_net_rate_sd * raw_trial_log_net_rate;
  }
  
  patient_log_net_rate_effect = patient_log_net_rate_sd * raw_patient_log_net_rate;
}

// Update the rate calculations to include covariate effects
patient_log_net_rate += patient_log_net_rate_linpred + trial_log_net_rate_effect[train_patient_trial] + patient_log_net_rate_effect;
  
vector[n_train_patients] patient_log_growth_rate = patient_log_net_rate - log_diff_exp(patient_log_rate_ratio, zeros_vector(n_train_patients));
vector[n_train_patients] patient_log_decrease_rate = patient_log_growth_rate + patient_log_rate_ratio;
  
vector[n_train_patients] patient_log_growth_lag_effect = zeros_vector(n_train_patients);
vector[n_train_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_train_patients);
  
if (!pop_growth_lag_param_only) {
  patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
  patient_log_growth_lag += patient_log_growth_lag_effect;
}
  
vector[n_train_patients] patient_decrease_prop_logis_linpred = Q_covar_design_matrix * QR_pop_decrease_prop_logis_coef;
vector[n_trials] trial_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_trials);
vector[n_trials] trial_decrease_prop_logis_effect = zeros_vector(n_trials);
vector[n_train_patients] patient_decrease_prop_logis = rep_vector(pop_decrease_prop_logis, n_train_patients) + patient_decrease_prop_logis_linpred; 
vector[n_train_patients] patient_decrease_prop_logis_effect = zeros_vector(n_train_patients);
  
if (!pop_initial_states_param_only ) {
  if (add_trial_level_prop) {
    trial_decrease_prop_logis_effect = trial_decrease_prop_logis_sd * raw_trial_decrease_prop_logis;
    trial_decrease_prop_logis += trial_decrease_prop_logis_effect;
  }
 
  patient_decrease_prop_logis_effect = patient_decrease_prop_logis_sd * raw_patient_decrease_prop_logis; 
  patient_decrease_prop_logis += trial_decrease_prop_logis_effect[patient_trial[train_patients_pos:train_patients_end]] + patient_decrease_prop_logis_effect;
} 
  
vector[n_train_patients] patient_log_decrease_prop = -log1p_exp(- patient_decrease_prop_logis);
vector[n_train_patients] patient_log_growth_prop = patient_log_decrease_prop - patient_decrease_prop_logis;
  
vector[n_train_patients] patient_tumor_gp_rho = independ_long_process_noise ? zeros_vector(n_train_patients) : rep_vector(exp(log_pop_tumor_gp_rho), n_train_patients);  
  
matrix[n_total_train_visits, 2] states; 
  
profile("states") {
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] log_patient_tumor_gp_rho_effect; 
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
  
  if (!independ_long_process_noise && !pop_rho_param_only) {
    log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
    log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
    patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial[train_patients_pos:train_patients_end]] + log_patient_tumor_gp_rho_effect);
  }

  states = calc_states(
    train_patient_visit_pos,
    train_patient_visits,
    patient_tumor_gp_rho,
    delta,
    pop_process_sd,
    independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
    // raw_patient_process_noise,
    independ_long_process_noise, independ_cross_process_noise,
    append_col(patient_log_decrease_prop, patient_log_growth_prop),
    exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
    rep_vector(0.0001, n_train_patients), // exp(patient_log_growth_lag), 
    0.0001, // exp(pop_log_growth_transition_rate),
    n_shards, // && !debug,
    0 // debug 
  );
}

vector[n_covar] pop_log_net_rate_coef = R_inv_covar_design_matrix * QR_pop_log_net_rate_coef;      
vector[n_covar] pop_decrease_prop_logis_coef = R_inv_covar_design_matrix * QR_pop_decrease_prop_logis_coef;
